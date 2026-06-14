//! On-the-wire message types.
//!
//! Framing is identical in spirit to the Python original: every message is
//! prefixed with a 4-byte big-endian length. The *payload*, however, is
//! `bincode` rather than Python `pickle`, so both ends must be the Rust build
//! (this is a from-scratch Rust protocol, not pickle-compatible).

use serde::{Deserialize, Serialize};

/// Outcome of running a task: either serialized bytes or an error string.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Outcome {
    /// `bincode`-serialized return value.
    Ok(Vec<u8>),
    /// Human-readable failure message (task returned Err or panicked).
    Err(String),
}

/// Messages exchanged over the TCP control channel between host and clients.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Wire {
    /// First framed message a client sends after connecting. Advertises the
    /// CPU count and the GPU device ids this client offers.
    Handshake {
        cpus: usize,
        gpus: Vec<usize>,
        hostname: String,
    },

    /// Host → client: please run `func(args)` under job `id`, pinned to the
    /// assigned `gpu_ids` (empty for CPU-only tasks).
    Task {
        id: u64,
        func: String,
        args: Vec<u8>,
        num_cpus: usize,
        gpu_ids: Vec<usize>,
    },

    /// Client → host: job `id` finished with `outcome` after `runtime_secs`.
    Result {
        id: u64,
        func: String,
        outcome: Outcome,
        runtime_secs: f64,
    },
}

/// A single node as seen by the status monitor.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeView {
    pub hostname: String,
    pub address: String,
    pub num_cpus: usize,
    pub free_cpus: usize,
    pub num_gpus: usize,
    pub free_gpus: usize,
}

/// A single job as seen by the status monitor.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JobView {
    pub id: u64,
    pub hostname: Option<String>,
    pub finished: bool,
    pub num_cpus: usize,
    pub num_gpus: usize,
    pub runtime_secs: f64,
    pub func: String,
}

/// What the UDP status server returns in response to a `STATUS` datagram.
///
/// `pending`/`running`/`done`/`total` summarize the *entire* queue, while
/// `queue` carries only a bounded sample (running jobs first) so the reply
/// always fits in a single UDP datagram even when there are millions of jobs.
/// `queue_truncated` is set when the sample omits entries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatusSnapshot {
    /// Seconds since the Unix epoch when the snapshot was taken.
    pub epoch_secs: u64,
    pub nodes: Vec<NodeView>,
    pub pending: usize,
    pub running: usize,
    pub done: usize,
    pub total: usize,
    /// Bounded sample of jobs (running first). May omit entries — see
    /// `queue_truncated`.
    pub queue: Vec<JobView>,
    pub queue_truncated: bool,
}

/// Maximum number of job entries included in a [`StatusSnapshot::queue`] sample.
/// Keeps the bincode'd datagram comfortably under the UDP size limit.
pub const STATUS_QUEUE_SAMPLE: usize = 100;

/// The literal datagram a client sends to request a [`StatusSnapshot`].
pub const STATUS_REQUEST: &[u8] = b"STATUS";
