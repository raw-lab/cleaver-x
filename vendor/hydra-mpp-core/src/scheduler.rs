//! The scheduler core: shared cluster state plus the background dispatch loop.
//!
//! Design notes
//! ------------
//! * All shared state lives behind a single `Mutex<Shared>`; a `Condvar`
//!   (`work_cv` for the scheduler, `done_cv` for waiters) are notified on the
//!   `wait()` blocks efficiently instead of busy-polling like the Python loop.
//! * Two resources are reserved per job: **CPU slots** (a count) and **GPU
//!   devices** (concrete ids). A job that asks for `num_cpus = N` reserves N
//!   slots on its node but still runs on a single worker thread — identical
//!   accounting to the Python original. A job that asks for `num_gpus = K` is
//!   pinned to K specific device ids for its lifetime.
//! * Each admitted local job runs on its own `std::thread`, with its assigned
//!   GPU ids placed in a thread-local so the task (and any subprocess it spawns)
//!   can read `CUDA_VISIBLE_DEVICES`. Panics inside a task are caught in the
//!   registry, so a bad job cannot bring the worker (or the node) down.

use std::collections::{BTreeMap, VecDeque};
use std::net::TcpStream;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use serde::de::DeserializeOwned;

use crate::error::{HydraError, Result};
use crate::job::{Completed, Job, JobId, JobResult, JobState, Payload};
use crate::node::{Node, NodeInfo, NodeLink};
use crate::registry::Registry;
use crate::wire::{JobView, NodeView, Outcome, StatusSnapshot, Wire};
use crate::{gpu, log, net, printlog};

/// Mutable cluster state guarded by [`Core`]'s mutex.
pub(crate) struct Shared {
    /// `nodes[0]` is always the local node.
    pub nodes: Vec<Node>,
    /// id → job. A `BTreeMap` keeps FIFO-ish ordering by id for fair dispatch.
    pub queue: BTreeMap<JobId, Job>,
    /// Ids of jobs still in `Pending`, in submission order. Dispatch pops from
    /// the front, so a pass costs O(jobs dispatched) instead of scanning the
    /// whole (possibly huge) queue every time.
    pub pending: VecDeque<JobId>,
    pub next_id: JobId,
    pub running: bool,
}

impl Shared {
    fn local_host(&self) -> String {
        self.nodes
            .first()
            .map(|n| n.hostname.clone())
            .unwrap_or_else(log::hostname)
    }
}

/// A unit of CPU-bound local work handed to a pool worker. The task handler is
/// resolved at dispatch time so the worker needs no registry lookup on the hot
/// path.
pub(crate) struct Work {
    pub id: JobId,
    pub func: Arc<str>,
    pub payload: Payload,
    pub gpus: Vec<usize>,
}

/// Shared runtime owned (via `Arc`) by every scheduler / worker / reader thread.
pub(crate) struct Core {
    pub state: Mutex<Shared>,
    /// The scheduler thread waits here. Notified when there is potentially new
    /// work to dispatch or freed resources (submit, finish, node join/loss,
    /// shutdown). It also wakes every 5 ms, so a missed notify only delays
    /// dispatch — it can never hang the scheduler.
    pub work_cv: Condvar,
    /// Result/dispatch waiters (`wait`, `get_blocking`, `await_dispatch`) wait
    /// here. Notified only when a job's observable state changes — it is
    /// dispatched (`Pending → Running`), completes, or the runtime stops. A
    /// `submit` no longer wakes these callers, since it completes no job.
    pub done_cv: Condvar,
    pub registry: Registry,
    /// Warm pool for **CPU-bound** local jobs: a shared queue + a condvar the
    /// workers park on. Only populated when `use_pool` is true; `cpus(0)` jobs
    /// never use it (they always get a dedicated thread).
    pub pool: (Mutex<VecDeque<Work>>, Condvar),
    /// Cleared on shutdown to wake idle pool workers so they can exit.
    pub pool_alive: AtomicBool,
    /// Whether CPU-bound jobs are routed through the warm pool (vs thread-per-job).
    pub use_pool: bool,
    /// True iff every node is local (mode == Local): enables the zero-copy
    /// typed payload path, which must never be used when a job could be shipped.
    pub local_only: bool,
    /// Optional stack size (bytes) for worker threads; `None` = OS default.
    pub worker_stack: Option<usize>,
}

pub(crate) type CoreRef = Arc<Core>;

impl Core {
    pub(crate) fn new(
        local: Node,
        registry: Registry,
        use_pool: bool,
        local_only: bool,
        worker_stack: Option<usize>,
    ) -> CoreRef {
        Arc::new(Core {
            state: Mutex::new(Shared {
                nodes: vec![local],
                queue: BTreeMap::new(),
                pending: VecDeque::new(),
                next_id: 0,
                running: true,
            }),
            work_cv: Condvar::new(),
            done_cv: Condvar::new(),
            registry,
            pool: (Mutex::new(VecDeque::new()), Condvar::new()),
            pool_alive: AtomicBool::new(true),
            use_pool,
            local_only,
            worker_stack,
        })
    }
}

// --------------------------------------------------------------------------
// Submission
// --------------------------------------------------------------------------

/// Queue a registered task for execution; returns its job id immediately.
pub(crate) fn submit(
    core: &CoreRef,
    func: &str,
    args: Vec<u8>,
    num_cpus: usize,
    num_gpus: usize,
) -> Result<JobId> {
    let fname = core
        .registry
        .intern(func)
        .ok_or_else(|| HydraError::UnknownTask(func.to_string()))?;
    let id = {
        let mut st = core.state.lock().unwrap();
        if !st.running {
            return Err(HydraError::NotRunning);
        }
        st.next_id += 1;
        let id = st.next_id;
        st.queue.insert(
            id,
            Job {
                id,
                func: fname,
                num_cpus,
                num_gpus,
                state: JobState::Pending {
                    args: Payload::Bytes(args),
                },
            },
        );
        st.pending.push_back(id);
        id
    };
    core.work_cv.notify_all();
    Ok(id)
}

/// Submit many jobs for the same task under a **single** lock acquisition and a
/// **single** wakeup. Calling `submit` in a loop costs O(n) lock cycles and O(n)
/// condvar `notify_all`s — and on a busy core each wakeup bounces the scheduler.
/// Batching collapses that to O(1); the win grows with batch size and powers
/// [`crate::Hydra::map`] / [`crate::Hydra::remote_many`].
pub(crate) fn submit_batch(
    core: &CoreRef,
    func: &str,
    args_list: Vec<Vec<u8>>,
    num_cpus: usize,
    num_gpus: usize,
) -> Result<Vec<JobId>> {
    let payloads = args_list.into_iter().map(Payload::Bytes).collect();
    submit_payloads(core, func, payloads, num_cpus, num_gpus)
}

/// The batch-submit primitive over `Payload`s (serialized bytes or typed
/// in-memory values). One lock acquisition + one wakeup for the whole batch.
pub(crate) fn submit_payloads(
    core: &CoreRef,
    func: &str,
    payloads: Vec<Payload>,
    num_cpus: usize,
    num_gpus: usize,
) -> Result<Vec<JobId>> {
    let fname = core
        .registry
        .intern(func)
        .ok_or_else(|| HydraError::UnknownTask(func.to_string()))?;
    let ids = {
        let mut st = core.state.lock().unwrap();
        if !st.running {
            return Err(HydraError::NotRunning);
        }
        let mut ids = Vec::with_capacity(payloads.len());
        for args in payloads {
            st.next_id += 1;
            let id = st.next_id;
            st.queue.insert(
                id,
                Job {
                    id,
                    func: fname.clone(),
                    num_cpus,
                    num_gpus,
                    state: JobState::Pending { args },
                },
            );
            st.pending.push_back(id);
            ids.push(id);
        }
        ids
    };
    if !ids.is_empty() {
        core.work_cv.notify_all();
    }
    Ok(ids)
}

/// Inject an already-computed value into the queue as a finished job.
pub(crate) fn put(core: &CoreRef, name: &str, value_bytes: Vec<u8>) -> Result<JobId> {
    let id = {
        let mut st = core.state.lock().unwrap();
        st.next_id += 1;
        let id = st.next_id;
        let host = st.local_host();
        st.queue.insert(
            id,
            Job {
                id,
                func: Arc::from(name),
                num_cpus: 0,
                num_gpus: 0,
                state: JobState::Done {
                    node: host,
                    gpus: Vec::new(),
                    completed: Completed::Serialized(Outcome::Ok(value_bytes)),
                    runtime_secs: 0.0,
                },
            },
        );
        id
    };
    core.done_cv.notify_all();
    Ok(id)
}

/// Block until `id` leaves the `Pending` state (used by `.blocking(true)`).
pub(crate) fn await_dispatch(core: &CoreRef, id: JobId) {
    let mut guard = core.state.lock().unwrap();
    loop {
        let dispatched = match guard.queue.get(&id) {
            Some(j) => !matches!(j.state, JobState::Pending { .. }),
            None => true,
        };
        if dispatched {
            return;
        }
        guard = core.done_cv.wait(guard).unwrap();
    }
}

// --------------------------------------------------------------------------
// Inspection: get / wait / nodes / snapshot
// --------------------------------------------------------------------------

fn to_result(job: &Job) -> JobResult {
    match &job.state {
        JobState::Pending { .. } => JobResult {
            finished: false,
            func_name: job.func.to_string(),
            num_cpus: job.num_cpus,
            num_gpus: job.num_gpus,
            gpu_ids: Vec::new(),
            runtime_secs: 0.0,
            hostname: None,
            outcome: None,
        },
        JobState::Running { node, gpus } => JobResult {
            finished: false,
            func_name: job.func.to_string(),
            num_cpus: job.num_cpus,
            num_gpus: job.num_gpus,
            gpu_ids: gpus.clone(),
            runtime_secs: 0.0,
            hostname: Some(node.clone()),
            outcome: None,
        },
        JobState::Done {
            node,
            gpus,
            completed,
            runtime_secs,
        } => JobResult {
            finished: true,
            func_name: job.func.to_string(),
            num_cpus: job.num_cpus,
            num_gpus: job.num_gpus,
            gpu_ids: gpus.clone(),
            runtime_secs: *runtime_secs,
            hostname: Some(node.clone()),
            // The bytes path exposes its `Outcome` here for `value()`. A typed
            // result is only retrievable via `get_typed` (it is not serializable
            // and `map_owned` never hands the caller a job id), so map it to None.
            outcome: match completed {
                Completed::Serialized(o) => Some(o.clone()),
                Completed::Typed(_) => None,
            },
        },
    }
}

/// Fetch a job's record. If finished, it is removed from the queue (matching the
/// Python `get()` which pops completed entries to free memory).
pub(crate) fn get(core: &CoreRef, id: JobId) -> Result<JobResult> {
    let mut st = core.state.lock().unwrap();
    let (res, finished) = {
        let job = st.queue.get(&id).ok_or(HydraError::UnknownJob(id))?;
        (to_result(job), job.is_finished())
    };
    if finished {
        st.queue.remove(&id);
    }
    Ok(res)
}

/// Consume a finished job and return its value as `R` — downcasting a typed
/// (zero-copy) result or deserializing a bytes result. Removes the job from the
/// queue like `get`. This is the retrieval half of the end-to-end zero-copy
/// path: a `map_owned` task keeps its result as an in-memory value, and this
/// hands it back with no deserialization. Errors if the job is unknown, not yet
/// finished, the task failed, or the requested type does not match.
pub(crate) fn take_value<R: DeserializeOwned + 'static>(core: &CoreRef, id: JobId) -> Result<R> {
    let mut st = core.state.lock().unwrap();
    match st.queue.get(&id) {
        Some(j) if j.is_finished() => {}
        Some(j) => {
            return Err(HydraError::TaskFailed {
                func: j.func.to_string(),
                message: "job not finished".into(),
            })
        }
        None => return Err(HydraError::UnknownJob(id)),
    }
    let job = st.queue.remove(&id).expect("present and finished");
    drop(st);
    let fname = job.func;
    match job.state {
        JobState::Done { completed, .. } => match completed {
            Completed::Typed(Ok(boxed)) => boxed.downcast::<R>().map(|b| *b).map_err(|_| {
                HydraError::TaskFailed {
                    func: fname.to_string(),
                    message: "typed result downcast failed (wrong type for get_typed?)".into(),
                }
            }),
            Completed::Typed(Err(msg)) => Err(HydraError::TaskFailed {
                func: fname.to_string(),
                message: msg,
            }),
            Completed::Serialized(Outcome::Ok(bytes)) => {
                bincode::deserialize(&bytes).map_err(HydraError::serde)
            }
            Completed::Serialized(Outcome::Err(msg)) => Err(HydraError::TaskFailed {
                func: fname.to_string(),
                message: msg,
            }),
        },
        _ => unreachable!("is_finished was checked under the lock"),
    }
}

fn collect(st: &Shared, pending: &mut Vec<JobId>, ready: &mut Vec<JobId>, max: usize) {
    let mut i = 0;
    while i < pending.len() && ready.len() < max {
        let id = pending[i];
        // A job that is no longer in the queue was already `get()`-ed → ready.
        let done = st.queue.get(&id).map(|j| j.is_finished()).unwrap_or(true);
        if done {
            ready.push(pending.remove(i));
        } else {
            i += 1;
        }
    }
}

/// Partition `ids` into `(ready, pending)`. Blocks (up to `timeout`) until at
/// least `max` jobs are ready or the pending set empties. Mirrors `wait()`.
pub(crate) fn wait(
    core: &CoreRef,
    ids: &[JobId],
    timeout: Option<Duration>,
    max: usize,
) -> (Vec<JobId>, Vec<JobId>) {
    let max = max.max(1);
    let mut pending = ids.to_vec();
    let mut ready = Vec::new();
    let start = Instant::now();

    let mut guard = core.state.lock().unwrap();
    collect(&guard, &mut pending, &mut ready, max);
    while ready.len() < max && !pending.is_empty() {
        match timeout {
            None => {
                guard = core.done_cv.wait(guard).unwrap();
            }
            Some(to) => {
                let elapsed = start.elapsed();
                if elapsed >= to {
                    break;
                }
                let (g, _) = core.done_cv.wait_timeout(guard, to - elapsed).unwrap();
                guard = g;
            }
        }
        collect(&guard, &mut pending, &mut ready, max);
        if let Some(to) = timeout {
            if start.elapsed() >= to {
                break;
            }
        }
    }
    (ready, pending)
}

/// Snapshot of all nodes (returned by `Hydra::nodes()`).
pub(crate) fn nodes(core: &CoreRef) -> Vec<NodeInfo> {
    let st = core.state.lock().unwrap();
    st.nodes
        .iter()
        .enumerate()
        .map(|(i, n)| NodeInfo {
            hostname: n.hostname.clone(),
            address: n.address.clone(),
            total_cpus: n.total_cpus,
            free_cpus: n.free_cpus,
            total_gpus: n.total_gpus(),
            free_gpus: n.free_gpus(),
            gpu_ids: n.gpu_ids.clone(),
            is_local: i == 0,
        })
        .collect()
}

/// Build a [`StatusSnapshot`] for the UDP status monitor.
///
/// Counts cover the whole queue; the per-job `queue` field is a bounded sample
/// (running jobs first) so the reply always fits in one UDP datagram.
pub(crate) fn snapshot(core: &CoreRef) -> StatusSnapshot {
    use crate::wire::STATUS_QUEUE_SAMPLE;

    let st = core.state.lock().unwrap();
    let epoch_secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let nodes = st
        .nodes
        .iter()
        .map(|n| NodeView {
            hostname: n.hostname.clone(),
            address: n.address.clone(),
            num_cpus: n.total_cpus,
            free_cpus: n.free_cpus,
            num_gpus: n.total_gpus(),
            free_gpus: n.free_gpus(),
        })
        .collect();

    // Whole-queue counts.
    let mut pending = 0usize;
    let mut running = 0usize;
    let mut done = 0usize;
    for j in st.queue.values() {
        match &j.state {
            JobState::Pending { .. } => pending += 1,
            JobState::Running { .. } => running += 1,
            JobState::Done { .. } => done += 1,
        }
    }
    let total = st.queue.len();

    // Bounded sample: running jobs first (most interesting), then the rest.
    let view_of = |j: &Job| {
        let (finished, hostname, runtime) = match &j.state {
            JobState::Pending { .. } => (false, None, 0.0),
            JobState::Running { node, .. } => (false, Some(node.clone()), 0.0),
            JobState::Done {
                node, runtime_secs, ..
            } => (true, Some(node.clone()), *runtime_secs),
        };
        JobView {
            id: j.id,
            hostname,
            finished,
            num_cpus: j.num_cpus,
            num_gpus: j.num_gpus,
            runtime_secs: runtime,
            func: j.func.to_string(),
        }
    };
    let mut queue: Vec<JobView> = Vec::new();
    for j in st.queue.values() {
        if matches!(j.state, JobState::Running { .. }) {
            queue.push(view_of(j));
            if queue.len() >= STATUS_QUEUE_SAMPLE {
                break;
            }
        }
    }
    if queue.len() < STATUS_QUEUE_SAMPLE {
        for j in st.queue.values() {
            if !matches!(j.state, JobState::Running { .. }) {
                queue.push(view_of(j));
                if queue.len() >= STATUS_QUEUE_SAMPLE {
                    break;
                }
            }
        }
    }
    let queue_truncated = total > queue.len();

    StatusSnapshot {
        epoch_secs,
        nodes,
        pending,
        running,
        done,
        total,
        queue,
        queue_truncated,
    }
}

// --------------------------------------------------------------------------
// Completion (called by local workers and remote reader threads)
// --------------------------------------------------------------------------

fn job_node_and_gpus(job: &Job, fallback: &str) -> (String, Vec<usize>) {
    match &job.state {
        JobState::Running { node, gpus } | JobState::Done { node, gpus, .. } => {
            (node.clone(), gpus.clone())
        }
        JobState::Pending { .. } => (fallback.to_string(), Vec::new()),
    }
}

/// Mark job `id` finished, return its reserved CPUs and GPUs to its node, and
/// wake waiters.
pub(crate) fn finish_job(core: &CoreRef, id: JobId, completed: Completed, runtime_secs: f64) {
    {
        let mut st = core.state.lock().unwrap();
        let fallback = st.local_host();
        let (node, ncpus, gpus) = match st.queue.get(&id) {
            Some(j) => {
                let (node, gpus) = job_node_and_gpus(j, &fallback);
                (node, j.num_cpus, gpus)
            }
            None => return, // already removed via get()
        };
        if let Some(n) = st.nodes.iter_mut().find(|n| n.hostname == node) {
            n.free_cpus = (n.free_cpus + ncpus).min(n.total_cpus);
            n.release_gpus(&gpus);
        }
        if let Some(j) = st.queue.get_mut(&id) {
            j.state = JobState::Done {
                node,
                gpus,
                completed,
                runtime_secs,
            };
        }
    }
    // A completion both frees resources (wake the scheduler) and finishes a job
    // (wake result/dispatch waiters).
    core.work_cv.notify_all();
    core.done_cv.notify_all();
}

// --------------------------------------------------------------------------
// Dispatch loop
// --------------------------------------------------------------------------

/// One dispatch pass: place pending jobs (from the front of the FIFO) onto any
/// node with free CPU **and** GPU capacity, stopping as soon as the head job
/// cannot be placed.
///
/// * A job whose `num_cpus`/`num_gpus` exceeds *every* node's total capacity can
///   never run; it is failed immediately with a clear error rather than wedging
///   the queue.
/// * Otherwise, if the head job does not currently fit anywhere, the pass stops
///   and waits for resources to free (strict FIFO — predictable, deadlock-free).
fn dispatch_pass(core: &CoreRef) {
    let mut local_jobs: Vec<(JobId, Arc<str>, Payload, Vec<usize>, usize)> = Vec::new();
    let mut remote_jobs: Vec<(Arc<Mutex<TcpStream>>, Wire)> = Vec::new();
    let mut unplaceable: Vec<(JobId, usize, usize)> = Vec::new();
    let mut mismatched: Vec<JobId> = Vec::new();

    {
        let mut st = core.state.lock().unwrap();
        if !st.running {
            return;
        }

        loop {
            let Some(&id) = st.pending.front() else { break };

            // Look up the head job's requirements (skip stale ids defensively).
            // The payload is *not* read here — it is moved out below, only once
            // we have committed to dispatching, so a failed fit leaves it intact.
            let (func, num_cpus, num_gpus) = match st.queue.get(&id) {
                Some(Job {
                    func,
                    num_cpus,
                    num_gpus,
                    state: JobState::Pending { .. },
                    ..
                }) => (func.clone(), *num_cpus, *num_gpus),
                _ => {
                    st.pending.pop_front();
                    continue;
                }
            };

            // Can any node ever run this job (enough *total* CPUs and GPUs)?
            let ever_fits = st
                .nodes
                .iter()
                .any(|n| n.total_cpus >= num_cpus && n.total_gpus() >= num_gpus);
            if !ever_fits {
                st.pending.pop_front();
                unplaceable.push((id, num_cpus, num_gpus));
                continue;
            }

            // Does it fit on a node *right now*?
            let chosen = st.nodes.iter().position(|n| {
                n.total_cpus > 0 && n.free_cpus >= num_cpus && n.free_gpus() >= num_gpus
            });
            let Some(i) = chosen else { break }; // head can't fit yet → wait

            // Reserve CPUs + concrete GPU ids on the chosen node.
            let gpu_ids = match st.nodes[i].reserve_gpus(num_gpus) {
                Some(ids) => ids,
                None => break, // shouldn't happen given the fit check
            };
            st.nodes[i].free_cpus -= num_cpus;
            let host = st.nodes[i].hostname.clone();
            let remote_write = match &st.nodes[i].link {
                NodeLink::Local => None,
                NodeLink::Remote { write } => Some(write.clone()),
            };
            st.pending.pop_front();

            // Move the payload out of Pending as we flip the job to Running.
            let payload = match st.queue.get_mut(&id) {
                Some(j) => match std::mem::replace(
                    &mut j.state,
                    JobState::Running {
                        node: host,
                        gpus: gpu_ids.clone(),
                    },
                ) {
                    JobState::Pending { args } => args,
                    other => {
                        j.state = other; // not actually pending; restore
                        continue;
                    }
                },
                None => continue,
            };

            match remote_write {
                None => local_jobs.push((id, func, payload, gpu_ids, num_cpus)),
                Some(write) => match payload {
                    Payload::Bytes(args) => remote_jobs.push((
                        write,
                        Wire::Task {
                            id,
                            func: func.to_string(),
                            args,
                            num_cpus,
                            gpu_ids,
                        },
                    )),
                    // Typed payloads only exist in local-only runtimes, so this
                    // is unreachable in practice; fail safe rather than panic.
                    Payload::Typed(_) => mismatched.push(id),
                },
            }
        }
    } // release state lock before doing side effects

    for (id, ncpus, ngpus) in unplaceable {
        finish_job(
            core,
            id,
            Completed::err(format!(
                "job reserves {ncpus} CPU(s) and {ngpus} GPU(s) but no node provides that many"
            )),
            0.0,
        );
    }
    for id in mismatched {
        finish_job(
            core,
            id,
            Completed::err("internal: typed payload cannot be sent to a remote node".into()),
            0.0,
        );
    }
    let dispatched = !local_jobs.is_empty() || !remote_jobs.is_empty();
    for (id, func, payload, gpus, ncpus) in local_jobs {
        // cpus(0) work always gets its own thread (unbounded I/O- / GPU-bound
        // concurrency); CPU-bound work uses the warm pool when it is enabled.
        if ncpus >= 1 && core.use_pool {
            enqueue_local(core, id, func, payload, gpus);
        } else {
            spawn_local(core, id, func, payload, gpus);
        }
    }
    for (write, task) in remote_jobs {
        if let Ok(mut s) = write.lock() {
            if let Err(e) = net::send_msg(&mut s, &task) {
                printlog!("WARN: failed to ship task to client: {e}");
            }
        }
    }
    // Jobs just moved Pending → Running: wake `await_dispatch` (`.blocking(true)`)
    // promptly instead of relying on the next unrelated notification.
    if dispatched {
        core.done_cv.notify_all();
    }
}

/// Resolve and run a job's payload. The `Typed` path is end-to-end zero-copy
/// (no argument decode, no result encode); `Bytes` is the serialized wire path.
fn run_payload(core: &CoreRef, func: &str, payload: Payload) -> Completed {
    match payload {
        Payload::Bytes(b) => match core.registry.get(func) {
            Some(t) => Completed::Serialized(t(&b)),
            None => Completed::err(format!("unknown task '{func}'")),
        },
        Payload::Typed(boxed) => match core.registry.get_typed(func) {
            Some(t) => Completed::Typed(t(boxed)),
            None => Completed::Typed(Err(format!("unknown task '{func}'"))),
        },
    }
}

/// Run a single local job on its own thread, pinned to `gpus`.
fn spawn_local(core: &CoreRef, id: JobId, func: Arc<str>, payload: Payload, gpus: Vec<usize>) {
    let core2 = core.clone();
    let mut builder = thread::Builder::new().name(format!("hydra-worker-{id}"));
    if let Some(sz) = core.worker_stack {
        builder = builder.stack_size(sz);
    }
    let spawned = builder.spawn(move || {
        gpu::set_assigned(&gpus);
        let start = Instant::now();
        let outcome = run_payload(&core2, &func, payload);
        finish_job(&core2, id, outcome, start.elapsed().as_secs_f64());
    });
    if let Err(e) = spawned {
        // Could not create the worker thread (e.g. resource limits). Fail the
        // job instead of leaving it stuck in Running forever — this releases its
        // reserved CPUs/GPUs and unblocks any wait()/get_blocking() callers.
        finish_job(
            core,
            id,
            Completed::err(format!("failed to spawn worker thread: {e}")),
            0.0,
        );
    }
}

/// Hand a CPU-bound local job to the warm pool (no per-job thread creation).
/// Used only when `core.use_pool` is true and `num_cpus >= 1`.
fn enqueue_local(core: &CoreRef, id: JobId, func: Arc<str>, payload: Payload, gpus: Vec<usize>) {
    let (lock, cvar) = &core.pool;
    lock.lock().unwrap().push_back(Work {
        id,
        func,
        payload,
        gpus,
    });
    cvar.notify_one();
}

/// One persistent pool worker: pull CPU-bound work, run it (pinned to its GPUs),
/// report. Exits when the pool is drained and the runtime is shutting down.
fn pool_worker(core: CoreRef) {
    let (lock, cvar) = &core.pool;
    loop {
        let work = {
            let mut q = lock.lock().unwrap();
            while q.is_empty() {
                if !core.pool_alive.load(Ordering::Acquire) {
                    return;
                }
                q = cvar.wait(q).unwrap();
            }
            q.pop_front()
        };
        if let Some(w) = work {
            gpu::set_assigned(&w.gpus);
            let start = Instant::now();
            let outcome = run_payload(&core, &w.func, w.payload);
            finish_job(&core, w.id, outcome, start.elapsed().as_secs_f64());
        }
    }
}

/// Start `size` persistent pool workers (called once from `init` when the warm
/// pool is enabled). Admission control caps concurrent CPU-reserving jobs at the
/// local CPU count, so a pool of that size never backs up for leaf tasks.
///
/// Note: jobs routed through the pool (those reserving >= 1 CPU) should be
/// leaves — a pooled task that submits work and blocks on it can starve a fixed
/// pool. Submit such work from the driver, or give it `cpus(0)` (which always
/// runs on its own thread).
pub(crate) fn start_pool(core: &CoreRef, size: usize) {
    for w in 0..size.max(1) {
        let core = core.clone();
        let mut builder = thread::Builder::new().name(format!("hydra-pool-{w}"));
        if let Some(sz) = core.worker_stack {
            builder = builder.stack_size(sz);
        }
        let _ = builder.spawn(move || pool_worker(core));
    }
}

/// The background scheduler thread. Returns when `running` is cleared.
pub(crate) fn run_scheduler(core: CoreRef) {
    printlog!("scheduler started");
    loop {
        {
            let st = core.state.lock().unwrap();
            if !st.running {
                break;
            }
        }
        dispatch_pass(&core);

        let guard = core.state.lock().unwrap();
        if !guard.running {
            break;
        }
        // Wake on the next submission/completion, or retry in 5 ms in case a
        // slot freed up while we were outside the lock.
        let _ = core
            .work_cv
            .wait_timeout(guard, Duration::from_millis(5))
            .unwrap();
    }
    printlog!("scheduler stopped");
}

/// Spawn a reader thread for a connected client that feeds results back.
pub(crate) fn spawn_peer_reader(core: CoreRef, mut read_stream: TcpStream, peer_host: String) {
    let _ = thread::Builder::new()
        .name(format!("hydra-reader-{peer_host}"))
        .spawn(move || {
            loop {
                match net::recv_msg(&mut read_stream) {
                    Ok(Some(Wire::Result {
                        id,
                        func: _,
                        outcome,
                        runtime_secs,
                    })) => finish_job(&core, id, Completed::Serialized(outcome), runtime_secs),
                    Ok(Some(_)) => {} // ignore unexpected frames
                    Ok(None) => {
                        printlog!("INFO: client {peer_host} disconnected");
                        break;
                    }
                    Err(e) => {
                        printlog!("WARN: read error from {peer_host}: {e}");
                        break;
                    }
                }
            }
            // The client is gone. Stop scheduling onto it, and — unlike the
            // Python original, which would leave its in-flight jobs hung forever
            // — fail those jobs explicitly so wait()/get() callers unblock with a
            // clear error instead of deadlocking.
            {
                let mut st = core.state.lock().unwrap();
                let dead: Vec<JobId> = st
                    .queue
                    .values()
                    .filter_map(|j| match &j.state {
                        JobState::Running { node, .. } if *node == peer_host => Some(j.id),
                        _ => None,
                    })
                    .collect();
                for id in dead {
                    if let Some(j) = st.queue.get_mut(&id) {
                        let gpus = match &j.state {
                            JobState::Running { gpus, .. } => gpus.clone(),
                            _ => Vec::new(),
                        };
                        j.state = JobState::Done {
                            node: peer_host.clone(),
                            gpus,
                            completed: Completed::err(format!(
                                "node '{peer_host}' disconnected before returning this job"
                            )),
                            runtime_secs: 0.0,
                        };
                    }
                }
                if let Some(n) = st.nodes.iter_mut().find(|n| n.hostname == peer_host) {
                    n.total_cpus = 0;
                    n.free_cpus = 0;
                    n.gpu_ids.clear();
                    n.gpu_free.clear();
                }
            }
            // The node is gone: capacity changed (scheduler) and its in-flight
            // jobs were just failed (result/dispatch waiters).
            core.work_cv.notify_all();
            core.done_cv.notify_all();
        });
}
