//! # HydraMPP — many-headed parallel processing, in Rust
//!
//! A featherweight, Ray-style distributed task engine. Register plain Rust
//! functions, fire them with [`Hydra::remote`], and collect results with
//! [`Hydra::wait`] / [`Hydra::get`]. The **same binary** runs on a laptop, a
//! workstation, or an entire SLURM allocation — only the [`Config`] changes.
//!
//! ```no_run
//! use hydra_mpp_core::prelude::*;
//!
//! let hydra = Hydra::new();
//! hydra.register("square", |x: u64| x * x);
//! hydra.init(Config::local()).unwrap();
//!
//! let ids: Vec<_> = (0..8).map(|i| hydra.remote("square", &i).unwrap()).collect();
//! let (ready, _) = hydra.wait(&ids, None, ids.len()).unwrap();
//! for id in ready {
//!     let r = hydra.get(id).unwrap();
//!     println!("square -> {}", r.value::<u64>().unwrap());
//! }
//! hydra.shutdown();
//! ```
//!
//! ## How it maps to the Python original
//! * `@hydra.remote` decorator  → [`Hydra::register`] + [`Hydra::remote`]
//! * dispatch *by function name* → preserved (every node registers the same
//!   names), but now type-checked at the boundary
//! * `init(address=…)`           → [`Hydra::init`] with [`Config`]
//! * `--hydraMPP-*` flags         → [`Config::from_args`] (`--hydra-*`)
//! * `get` / `wait` / `put`       → same names, same semantics
//! * UDP status monitor           → preserved (see the `hydra-status` binary)

#![forbid(unsafe_code)]

mod banner;
mod cluster;
mod config;
mod error;
mod gpu;
mod job;
mod log;
mod net;
mod node;
mod registry;
mod scheduler;
mod slurm;
mod wire;

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Duration;

use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::node::Node;
use crate::registry::Registry;
use crate::scheduler::{Core, CoreRef};

pub use crate::config::{Config, Mode, PoolMode, DEFAULT_PORT};
pub use crate::error::{HydraError, Result};
pub use crate::gpu::{cuda_visible_devices, current_gpus};
pub use crate::job::{JobId, JobResult};
pub use crate::node::NodeInfo;
pub use crate::slurm::SlurmClients;
pub use crate::wire::StatusSnapshot;

/// Crate version, surfaced in the banner and the status monitor.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Query a running host's UDP status monitor and return its [`StatusSnapshot`].
///
/// This is what the `hydra-status` binary uses; it sends a single `STATUS`
/// datagram to `address:port` and decodes the reply.
pub fn query_status(address: &str, port: u16, timeout: Duration) -> Result<StatusSnapshot> {
    use std::net::UdpSocket;

    let socket = UdpSocket::bind(("0.0.0.0", 0))?;
    socket.set_read_timeout(Some(timeout))?;
    socket.send_to(wire::STATUS_REQUEST, (address, port))?;

    let mut buf = vec![0u8; 64 * 1024];
    let (n, _) = socket.recv_from(&mut buf)?;
    bincode::deserialize::<StatusSnapshot>(&buf[..n]).map_err(HydraError::serde)
}

/// Everything you need for everyday use.
///
/// ```
/// use hydra_mpp_core::prelude::*;
/// ```
pub mod prelude {
    pub use crate::{Config, Hydra, JobId, JobResult, NodeInfo, PoolMode};
}

/// The HydraMPP runtime handle.
///
/// Create one with [`Hydra::new`], register tasks, then [`Hydra::init`]. The
/// handle is `Send + Sync`, so you can wrap it in an `Arc` and submit work from
/// many threads.
pub struct Hydra {
    registry: Registry,
    inner: Mutex<Option<CoreRef>>,
    started: AtomicBool,
    slurm: Mutex<Option<SlurmClients>>,
}

impl Default for Hydra {
    fn default() -> Self {
        Self::new()
    }
}

impl Hydra {
    /// Create a new, un-initialized runtime. Register tasks before calling
    /// [`Hydra::init`].
    pub fn new() -> Self {
        Hydra {
            registry: Registry::new(),
            inner: Mutex::new(None),
            started: AtomicBool::new(false),
            slurm: Mutex::new(None),
        }
    }

    /// Register a task under `name`.
    ///
    /// The handler takes one argument that implements [`serde::de::DeserializeOwned`]
    /// and returns a value that implements [`serde::Serialize`]. For several
    /// arguments, take a tuple: `|(a, b): (u64, u64)| a + b`.
    ///
    /// Register the same names on every node so the host can dispatch by name —
    /// this mirrors the Python global `WORKERS` dict, now type-checked.
    pub fn register<In, Out, F>(&self, name: impl Into<String>, f: F)
    where
        In: DeserializeOwned + Send + 'static,
        Out: Serialize + Send + 'static,
        F: Fn(In) -> Out + Send + Sync + 'static,
    {
        self.registry.register(name, f);
    }

    /// Start the runtime in the role described by `cfg`.
    ///
    /// * [`Mode::Local`] / [`Mode::Host`] / [`Mode::SlurmHost`] return after
    ///   startup (a background scheduler keeps running).
    /// * [`Mode::Client`] **does not return**: it serves tasks until the host
    ///   disconnects, then exits the process with status 0 — exactly like the
    ///   Python client (`sys.exit(0)`).
    pub fn init(&self, cfg: Config) -> Result<()> {
        if self.started.swap(true, Ordering::SeqCst) {
            return Err(HydraError::AlreadyRunning);
        }
        log::set_quiet(cfg.quiet);

        let local_cpus = cfg
            .num_cpus
            .unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            })
            .max(1);
        let local_gpus = gpu::detect(cfg.num_gpus);

        // ---- Client: serve forever, then exit (never returns) -------------
        if let Mode::Client { address } = &cfg.mode {
            if !cfg.quiet {
                eprint!(
                    "{}",
                    banner::render(VERSION, "client", local_cpus, local_gpus.len())
                );
            }
            let res = cluster::client_serve(
                self.registry.clone(),
                address,
                cfg.port,
                local_cpus,
                local_gpus,
            );
            if let Err(e) = res {
                printlog!("CLIENT ERROR: {e}");
                std::process::exit(1);
            }
            std::process::exit(0);
        }

        // ---- SLURM head: launch clients, then behave as a host ------------
        let mode_label = match &cfg.mode {
            Mode::Local => "local",
            Mode::Host => "host",
            Mode::SlurmHost { .. } => "slurm-host",
            Mode::Client { .. } => unreachable!(),
        };
        if !cfg.quiet {
            eprint!(
                "{}",
                banner::render(VERSION, mode_label, local_cpus, local_gpus.len())
            );
        }

        if let Mode::SlurmHost { nodelist } = &cfg.mode {
            match slurm::spawn_clients(nodelist, cfg.rest_args(), cfg.port, cfg.num_cpus, cfg.num_gpus)
            {
                Ok((head, clients)) => {
                    printlog!("SLURM: this node ({head}) is the host");
                    *self.slurm.lock().unwrap() = Some(clients);
                }
                Err(e) => {
                    printlog!("SLURM ERROR: {e}");
                    self.started.store(false, Ordering::SeqCst);
                    return Err(e);
                }
            }
        }

        // ---- Local / Host / SlurmHost: build core + background threads ----
        // Decide CPU-bound execution strategy. cpus(0) jobs always thread-per-job.
        let use_pool = match cfg.pool {
            PoolMode::Always => true,
            PoolMode::Never => false,
            // The warm pool's per-job hand-off only pays off with real
            // parallelism; on a single core a thread spawn is cheaper.
            PoolMode::Auto => local_cpus > 1,
        };
        let local = Node::local(log::hostname(), local_cpus, local_gpus.clone());
        // Zero-copy typed payloads are only safe when no job can be shipped.
        let local_only = matches!(cfg.mode, Mode::Local);
        let core = Core::new(
            local,
            self.registry.clone(),
            use_pool,
            local_only,
            cfg.worker_stack_bytes,
        );
        *self.inner.lock().unwrap() = Some(core.clone());

        if use_pool {
            scheduler::start_pool(&core, local_cpus);
        }

        printlog!(
            "starting HydraMPP v{VERSION} ({mode_label}, {local_cpus} local CPUs, {} local GPUs {:?}, cpu-pool={})",
            local_gpus.len(),
            local_gpus,
            if use_pool { "on" } else { "off" }
        );

        // Background scheduler.
        {
            let core = core.clone();
            let _ = thread::Builder::new()
                .name("hydra-scheduler".into())
                .spawn(move || scheduler::run_scheduler(core));
        }
        // UDP status monitor (host + local).
        cluster::spawn_status_server(core.clone(), cfg.port);

        // Host / SlurmHost: synchronously gather clients for the accept window.
        if matches!(cfg.mode, Mode::Host | Mode::SlurmHost { .. }) {
            if let Err(e) = cluster::host_accept(&core, cfg.port, cfg.timeout) {
                printlog!("HOST WARN: accept loop ended early: {e}");
            }
        }
        Ok(())
    }

    /// Submit `func(arg)` for execution and return its job id immediately
    /// (non-blocking). Reserves one CPU slot. For more control use [`Hydra::task`].
    pub fn remote<A: Serialize>(&self, func: &str, arg: &A) -> Result<JobId> {
        self.task(func).submit(arg)
    }

    /// Submit `func` over a whole slice of arguments in **one** scheduler lock
    /// cycle and **one** wakeup (vs `remote` in a loop, which costs one of each
    /// per item). Returns the job ids in input order without waiting — collect
    /// them later with [`Hydra::wait`] / [`Hydra::get`], or use [`Hydra::map`] to
    /// submit-wait-collect in a single call. Each job reserves `cpus_each` CPUs
    /// (`0` for I/O- or GPU-bound work).
    pub fn remote_many<A: Serialize>(
        &self,
        func: &str,
        args: &[A],
        cpus_each: usize,
    ) -> Result<Vec<JobId>> {
        let core = self.core()?;
        let mut blobs = Vec::with_capacity(args.len());
        for a in args {
            blobs.push(bincode::serialize(a).map_err(HydraError::serde)?);
        }
        scheduler::submit_batch(&core, func, blobs, cpus_each, 0)
    }

    /// Begin a customized submission (CPU reservation, dispatch-blocking).
    ///
    /// ```no_run
    /// # use hydra_mpp_core::prelude::*;
    /// # let hydra = Hydra::new();
    /// # hydra.register("heavy", |x: u64| x);
    /// # hydra.init(Config::local()).unwrap();
    /// let id = hydra.task("heavy").cpus(4).submit(&42).unwrap();
    /// ```
    pub fn task<'h>(&'h self, func: &str) -> TaskCall<'h> {
        TaskCall {
            hydra: self,
            func: func.to_string(),
            num_cpus: 1,
            num_gpus: 0,
            blocking: false,
        }
    }

    /// Inject an already-computed value as a finished job, returning its id.
    /// Mirrors the Python `put()` (handy for seeding shared constants).
    pub fn put<A: Serialize>(&self, name: &str, value: &A) -> Result<JobId> {
        let core = self.core()?;
        let bytes = bincode::serialize(value).map_err(HydraError::serde)?;
        scheduler::put(&core, name, bytes)
    }

    /// Partition `ids` into `(ready, pending)`.
    ///
    /// Blocks until at least `max` jobs are ready, the pending set empties, or
    /// `timeout` elapses (`None` = wait indefinitely). The global queue is left
    /// untouched — call [`Hydra::get`] to actually retrieve a result.
    pub fn wait(
        &self,
        ids: &[JobId],
        timeout: Option<Duration>,
        max: usize,
    ) -> Result<(Vec<JobId>, Vec<JobId>)> {
        let core = self.core()?;
        Ok(scheduler::wait(&core, ids, timeout, max))
    }

    /// Fetch a job's [`JobResult`]. If the job has finished it is removed from
    /// the queue (matching the Python `get()`, which pops completed entries).
    pub fn get(&self, id: JobId) -> Result<JobResult> {
        let core = self.core()?;
        scheduler::get(&core, id)
    }

    /// Convenience: block until job `id` is finished, then return its value.
    pub fn get_blocking<R: DeserializeOwned>(&self, id: JobId) -> Result<R> {
        let _ = self.wait(&[id], None, 1)?;
        self.get(id)?.value::<R>()
    }

    /// Block until job `id` is finished, then **consume** it and return its
    /// value — downcasting a typed (zero-copy) result or deserializing a bytes
    /// result. This is the retrieval counterpart to [`map_owned`](Self::map_owned):
    /// a typed result is handed back with no deserialization.
    pub fn get_typed<R: DeserializeOwned + 'static>(&self, id: JobId) -> Result<R> {
        let core = self.core()?;
        let _ = self.wait(&[id], None, 1)?;
        scheduler::take_value::<R>(&core, id)
    }

    /// Stream a function over `args` with **bounded memory**. At most `window`
    /// jobs are in flight — and at most `window` results held — at once, so a
    /// fan-out over millions of inputs (or with large per-task results) stays in
    /// a fixed memory envelope instead of materializing every result the way
    /// [`map`](Self::map) does.
    ///
    /// Results are delivered to `f(index, result)` in **input order** as they
    /// complete, and each result is dropped after the callback returns — freeing
    /// its memory before the next is collected. A failing task is passed to `f`
    /// as an `Err` (the batch is not aborted; `f` decides what to do per item).
    ///
    /// This is backpressure: submission is paced to the completion rate, capping
    /// both in-flight work and peak memory. Use it for huge fan-outs or large
    /// results; use [`map`](Self::map) when you want every result returned at once.
    ///
    /// ```no_run
    /// # use hydra_mpp_core::prelude::*;
    /// # let hydra = Hydra::new();
    /// # hydra.register("work", |x: u64| x * 2);
    /// # hydra.init(Config::local()).unwrap();
    /// let inputs: Vec<u64> = (0..1_000_000).collect();
    /// let mut total = 0u64;
    /// hydra.for_each::<u64, u64, _>("work", &inputs, 256, |_i, r| {
    ///     if let Ok(v) = r { total += v; }
    /// }).unwrap();
    /// ```
    pub fn for_each<In, Out, F>(&self, func: &str, args: &[In], window: usize, mut f: F) -> Result<()>
    where
        In: Serialize,
        Out: DeserializeOwned,
        F: FnMut(usize, Result<Out>),
    {
        let n = args.len();
        if n == 0 {
            return Ok(());
        }
        let window = window.clamp(1, n);
        let mut ids: Vec<Option<JobId>> = (0..n).map(|_| None).collect();
        // Prime the window: submit the first `window` jobs.
        let mut submitted = 0usize;
        while submitted < window {
            ids[submitted] = Some(self.remote(func, &args[submitted])?);
            submitted += 1;
        }
        // Deliver in input order; refill one job as each is delivered.
        for deliver in 0..n {
            let id = ids[deliver].take().expect("job was submitted");
            let res = self.get_blocking::<Out>(id);
            f(deliver, res); // result moves into `f` and is dropped when it returns
            if submitted < n {
                ids[submitted] = Some(self.remote(func, &args[submitted])?);
                submitted += 1;
            }
        }
        Ok(())
    }

    /// Run `func` over every element of `args` **in parallel across the whole
    /// cluster**, returning results in input order. Each task reserves one CPU
    /// slot; use [`Hydra::map_cpus`] to change that. A single failing task
    /// aborts the batch with its error.
    ///
    /// This is the batch ergonomic analogue of rayon's `par_iter().map()` — but
    /// the work fans out across nodes (and can pin GPUs per task) instead of
    /// staying on one machine. For fine-grained *local* data parallelism, call
    /// rayon *inside* `func`; the two compose.
    ///
    /// ```no_run
    /// use hydra_mpp_core::prelude::*;
    /// let hydra = Hydra::new();
    /// hydra.register("sq", |x: u64| x * x);
    /// hydra.init(Config::local().quiet(true)).unwrap();
    /// let squares: Vec<u64> = hydra.map("sq", &[1u64, 2, 3, 4]).unwrap();
    /// assert_eq!(squares, vec![1, 4, 9, 16]);
    /// ```
    pub fn map<In, Out>(&self, func: &str, args: &[In]) -> Result<Vec<Out>>
    where
        In: Serialize,
        Out: DeserializeOwned,
    {
        self.map_cpus(func, args, 1)
    }

    /// Like [`Hydra::map`], but each task reserves `cpus_each` CPU slots. Pass
    /// `0` for I/O- or GPU-bound work that should not be throttled by CPU slots.
    pub fn map_cpus<In, Out>(&self, func: &str, args: &[In], cpus_each: usize) -> Result<Vec<Out>>
    where
        In: Serialize,
        Out: DeserializeOwned,
    {
        let core = self.core()?;
        // One lock acquisition + one wakeup for the whole batch (see remote_many).
        let ids = self.remote_many(func, args, cpus_each)?;
        let (_ready, pending) = scheduler::wait(&core, &ids, None, ids.len());
        if !pending.is_empty() {
            return Err(HydraError::Protocol(format!(
                "{} task(s) did not complete (runtime shut down?)",
                pending.len()
            )));
        }
        let mut out = Vec::with_capacity(ids.len());
        for id in ids {
            out.push(self.get(id)?.value::<Out>()?);
        }
        Ok(out)
    }

    /// Like [`map`](Self::map), but **moves** owned arguments in and takes the
    /// zero-copy local fast path: on a local-only runtime each argument is handed
    /// to its task as an in-memory value, skipping the serialize → store →
    /// deserialize round-trip that `map` pays. The win scales with argument size
    /// — for large payloads (sequence chunks, matrices) it removes the dominant
    /// per-task cost. On a distributed runtime (a host with clients) it falls
    /// back to serializing, since the value may cross the wire. Results are
    /// returned in input order; the first failure aborts the batch.
    ///
    /// Use this when you own the inputs and they are large; use `map` when you
    /// only have a borrowed slice or the inputs are small.
    pub fn map_owned<In, Out>(&self, func: &str, args: Vec<In>) -> Result<Vec<Out>>
    where
        In: Serialize + Send + 'static,
        Out: DeserializeOwned + 'static,
    {
        let core = self.core()?;
        let payloads: Vec<job::Payload> = if core.local_only {
            // Zero-copy: move each owned value straight into the queue.
            args.into_iter()
                .map(|a| job::Payload::Typed(Box::new(a)))
                .collect()
        } else {
            // Distributed runtime: the value may be shipped, so serialize it.
            let mut v = Vec::with_capacity(args.len());
            for a in &args {
                v.push(job::Payload::Bytes(
                    bincode::serialize(a).map_err(HydraError::serde)?,
                ));
            }
            v
        };
        let ids = scheduler::submit_payloads(&core, func, payloads, 1, 0)?;
        let (_ready, pending) = scheduler::wait(&core, &ids, None, ids.len());
        if !pending.is_empty() {
            return Err(HydraError::Protocol(format!(
                "{} task(s) did not complete (runtime shut down?)",
                pending.len()
            )));
        }
        let mut out = Vec::with_capacity(ids.len());
        for id in ids {
            // Zero-copy on the result side too: downcast the in-memory value
            // (local) or deserialize the bytes (distributed fallback).
            out.push(scheduler::take_value::<Out>(&core, id)?);
        }
        Ok(out)
    }

    /// A snapshot of every node in the cluster.
    pub fn nodes(&self) -> Result<Vec<NodeInfo>> {
        let core = self.core()?;
        Ok(scheduler::nodes(&core))
    }

    /// Whether the runtime is currently running.
    pub fn is_running(&self) -> bool {
        self.inner
            .lock()
            .unwrap()
            .as_ref()
            .map(|c| c.state.lock().unwrap().running)
            .unwrap_or(false)
    }

    /// Stop the runtime: signal the scheduler, status server, and reader
    /// threads to wind down. Idempotent. After this, [`Hydra::init`] may not be
    /// called again on the same handle.
    pub fn shutdown(&self) {
        if let Some(core) = self.inner.lock().unwrap().as_ref() {
            {
                let mut st = core.state.lock().unwrap();
                st.running = false;
            }
            core.pool_alive.store(false, Ordering::SeqCst);
            core.work_cv.notify_all();
            core.done_cv.notify_all();
            core.pool.1.notify_all();
            printlog!("shutdown requested");
        }
    }

    // -- internal helpers ---------------------------------------------------

    fn core(&self) -> Result<CoreRef> {
        self.inner
            .lock()
            .unwrap()
            .as_ref()
            .cloned()
            .ok_or(HydraError::NotRunning)
    }
}

impl Drop for Hydra {
    fn drop(&mut self) {
        self.shutdown();
    }
}

/// Builder returned by [`Hydra::task`] for a customized submission.
pub struct TaskCall<'h> {
    hydra: &'h Hydra,
    func: String,
    num_cpus: usize,
    num_gpus: usize,
    blocking: bool,
}

impl<'h> TaskCall<'h> {
    /// Reserve `n` CPU slots for this job (the Python `num_cpus=`). The job
    /// still runs on a single worker, but the node's free-CPU count drops by
    /// `n`, throttling co-scheduled work. Pass `0` for GPU- or I/O-bound tasks
    /// that should be gated only by their GPU reservation, not by CPU slots
    /// (note: a job reserving `0` CPUs **and** `0` GPUs is unthrottled).
    pub fn cpus(mut self, n: usize) -> Self {
        self.num_cpus = n;
        self
    }

    /// Reserve `n` GPU devices for this job. The scheduler pins `n` concrete
    /// device ids on whichever node runs it; the task reads them with
    /// [`current_gpus`] / [`cuda_visible_devices`]. A job that asks for more
    /// GPUs than any node has fails fast with a clear error.
    pub fn gpus(mut self, n: usize) -> Self {
        self.num_gpus = n;
        self
    }

    /// If `true`, [`TaskCall::submit`] blocks until the job is *dispatched*
    /// (not until it finishes) — the Python `blocking=True` semantics.
    pub fn blocking(mut self, blocking: bool) -> Self {
        self.blocking = blocking;
        self
    }

    /// Serialize `arg`, enqueue the job, and return its id.
    pub fn submit<A: Serialize>(self, arg: &A) -> Result<JobId> {
        let core = self.hydra.core()?;
        let bytes = bincode::serialize(arg).map_err(HydraError::serde)?;
        let id = scheduler::submit(&core, &self.func, bytes, self.num_cpus, self.num_gpus)?;
        if self.blocking {
            scheduler::await_dispatch(&core, id);
        }
        Ok(id)
    }
}
