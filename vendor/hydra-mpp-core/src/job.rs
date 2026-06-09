//! Jobs and their results.
//!
//! The Python original stores each queue entry as a 6-field list indexed
//! `[finished, func_name, ret, num_cpus, runtime, hostname]`. We keep the same
//! information but model the lifecycle as a proper enum so illegal states are
//! unrepresentable, and we additionally track the GPU devices a job reserved.

use std::any::Any;
use std::sync::Arc;

use serde::de::DeserializeOwned;

use crate::error::{HydraError, Result};
use crate::wire::Outcome;

/// Opaque identifier returned by `remote()` / `put()`.
pub type JobId = u64;

/// A task argument in the queue: either serialized bytes (for jobs that may be
/// shipped to another node) or a typed in-memory value (the zero-copy fast path
/// used by local-only runtimes — see `Hydra::map_owned`). `Typed` never crosses
/// the wire; it only exists when every node is local.
pub(crate) enum Payload {
    Bytes(Vec<u8>),
    Typed(Box<dyn Any + Send>),
}

impl std::fmt::Debug for Payload {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Payload::Bytes(b) => write!(f, "Bytes({} bytes)", b.len()),
            Payload::Typed(_) => f.write_str("Typed(<in-memory>)"),
        }
    }
}

/// Where a job is in its lifecycle.
pub(crate) enum JobState {
    /// Queued, not yet dispatched. Holds the (bytes or typed) argument.
    Pending { args: Payload },
    /// Dispatched to `node`, awaiting its result; pinned to `gpus`.
    Running { node: String, gpus: Vec<usize> },
    /// Completed (successfully or not) on `node`.
    Done {
        node: String,
        gpus: Vec<usize>,
        completed: Completed,
        runtime_secs: f64,
    },
}

/// The stored result of a finished job. `Serialized` is the wire path and every
/// bytes-argument task (`map`, `remote`, remote nodes). `Typed` is the local
/// zero-copy path used by `map_owned`: the worker keeps the result as an
/// in-memory value instead of serializing it, and `get_typed` downcasts it back.
/// `Typed` never crosses the wire.
pub(crate) enum Completed {
    Serialized(Outcome),
    Typed(std::result::Result<Box<dyn Any + Send>, String>),
}

impl Completed {
    /// A serialized error outcome (used by the bytes paths and failure cases).
    pub(crate) fn err(msg: String) -> Self {
        Completed::Serialized(Outcome::Err(msg))
    }
}

// Manual Debug: `Payload::Typed` holds `dyn Any`, which is not `Debug`.
impl std::fmt::Debug for JobState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            JobState::Pending { .. } => f.write_str("Pending"),
            JobState::Running { node, gpus } => {
                write!(f, "Running {{ node: {node:?}, gpus: {gpus:?} }}")
            }
            JobState::Done {
                node, runtime_secs, ..
            } => write!(f, "Done {{ node: {node:?}, runtime_secs: {runtime_secs} }}"),
        }
    }
}

/// Internal queue entry.
#[derive(Debug)]
pub(crate) struct Job {
    pub id: JobId,
    /// Interned task name (shared `Arc` — cloning is a refcount bump, not an alloc).
    pub func: Arc<str>,
    pub num_cpus: usize,
    pub num_gpus: usize,
    pub state: JobState,
}

impl Job {
    pub(crate) fn is_finished(&self) -> bool {
        matches!(self.state, JobState::Done { .. })
    }
}

/// The public, read-only record returned by `get()`.
///
/// Field order mirrors the Python `get()` tuple so existing mental models carry
/// over: `finished, func_name, ret, num_cpus, runtime, hostname`. Here `ret` is
/// returned lazily through [`JobResult::value`] so callers get their own type
/// back instead of raw bytes. GPU reservations are exposed too.
#[derive(Debug, Clone)]
pub struct JobResult {
    /// Index 0 — whether the job has completed.
    pub finished: bool,
    /// Index 1 — the name of the executed task.
    pub func_name: String,
    /// Index 3 — number of CPUs the job reserved.
    pub num_cpus: usize,
    /// Number of GPUs the job reserved.
    pub num_gpus: usize,
    /// The concrete GPU device ids the job was pinned to (empty if none).
    pub gpu_ids: Vec<usize>,
    /// Index 4 — wall-clock run time in seconds (0.0 until finished).
    pub runtime_secs: f64,
    /// Index 5 — the node the job ran on (None until dispatched).
    pub hostname: Option<String>,
    /// Index 2 — the raw serialized return value (None until finished / on error).
    pub(crate) outcome: Option<Outcome>,
}

impl JobResult {
    /// Deserialize the task's return value into `R`.
    ///
    /// Returns `Err` if the job is not finished yet, if the task itself failed,
    /// or if the bytes cannot be decoded as `R`.
    pub fn value<R: DeserializeOwned>(&self) -> Result<R> {
        match &self.outcome {
            None => Err(HydraError::TaskFailed {
                func: self.func_name.clone(),
                message: "job not finished".into(),
            }),
            Some(Outcome::Err(msg)) => Err(HydraError::TaskFailed {
                func: self.func_name.clone(),
                message: msg.clone(),
            }),
            Some(Outcome::Ok(bytes)) => bincode::deserialize(bytes).map_err(HydraError::serde),
        }
    }

    /// `true` if the job finished *and* the task returned successfully.
    pub fn is_ok(&self) -> bool {
        matches!(self.outcome, Some(Outcome::Ok(_)))
    }

    /// The task's error message, if it finished but failed.
    pub fn error(&self) -> Option<&str> {
        match &self.outcome {
            Some(Outcome::Err(m)) => Some(m.as_str()),
            _ => None,
        }
    }
}
