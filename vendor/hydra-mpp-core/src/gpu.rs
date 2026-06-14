//! GPU awareness.
//!
//! HydraMPP treats GPUs as a *reservable resource*, exactly like CPU slots: a
//! node advertises the device indices it owns, a task reserves `num_gpus` of
//! them with [`crate::Hydra::task`]`.gpus(n)`, and the scheduler pins concrete
//! device ids to the job for its lifetime.
//!
//! Crucially, **HydraMPP does not link CUDA**. It never calls a GPU API itself;
//! it only *schedules* and *pins*. Your task does the actual GPU work — whether
//! that is a Rust CUDA crate, PyTorch in a subprocess, or a GPU-accelerated
//! aligner. The runtime tells the task which devices it owns in two race-free
//! ways:
//!
//! * [`current_gpus`] — the device ids assigned to the *currently running*
//!   task (a thread-local, so concurrent GPU jobs never see each other's ids);
//! * [`cuda_visible_devices`] — the same ids as a comma-joined string, ready to
//!   drop into a subprocess: `.env("CUDA_VISIBLE_DEVICES", cuda_visible_devices())`.
//!
//! Detection is dependency-free and degrades gracefully to "no GPUs" when none
//! are present, so the same binary runs unchanged on CPU-only nodes.

use std::cell::RefCell;
use std::process::Command;

thread_local! {
    /// Device ids assigned to the task running on *this* worker thread.
    static ASSIGNED: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) };
}

/// Set the GPU ids for the task about to run on this thread. Called by the
/// scheduler / client just before invoking the task body.
pub(crate) fn set_assigned(ids: &[usize]) {
    ASSIGNED.with(|a| *a.borrow_mut() = ids.to_vec());
}

/// The GPU device ids assigned to the currently running task.
///
/// Returns an empty vector outside a GPU task (or for a task that reserved no
/// GPUs). Call this from inside a registered task to learn which devices it
/// owns:
///
/// ```ignore
/// hydra.register("infer", |batch: Vec<f32>| {
///     let devices = hydra_mpp_core::current_gpus(); // e.g. [2]
///     // ... run your kernel on `devices` ...
///     batch.len()
/// });
/// ```
pub fn current_gpus() -> Vec<usize> {
    ASSIGNED.with(|a| a.borrow().clone())
}

/// The currently assigned GPU ids formatted as a `CUDA_VISIBLE_DEVICES` string
/// (e.g. `"2,3"`). Empty string if none are assigned. Hand this to a subprocess:
///
/// ```no_run
/// use std::process::Command;
/// Command::new("my_gpu_tool")
///     .env("CUDA_VISIBLE_DEVICES", hydra_mpp_core::cuda_visible_devices())
///     .status()
///     .unwrap();
/// ```
pub fn cuda_visible_devices() -> String {
    current_gpus()
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// Parse a `CUDA_VISIBLE_DEVICES`-style value into device indices.
///
/// Handles the common integer form (`"0,1,3"`). If entries are non-numeric
/// (e.g. GPU UUIDs), falls back to logical slots `0..count` so scheduling still
/// works.
fn parse_visible(value: &str) -> Option<Vec<usize>> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Some(Vec::new());
    }
    let parts: Vec<&str> = trimmed.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
    if parts.is_empty() {
        return Some(Vec::new());
    }
    let ids: Result<Vec<usize>, _> = parts.iter().map(|s| s.parse::<usize>()).collect();
    match ids {
        Ok(v) => Some(v),
        // UUIDs or other tokens: use logical slots 0..count.
        Err(_) => Some((0..parts.len()).collect()),
    }
}

/// Query NVIDIA GPUs via `nvidia-smi` (no driver library needed at build time).
fn nvidia_smi_indices() -> Option<Vec<usize>> {
    let out = Command::new("nvidia-smi")
        .args(["--query-gpu=index", "--format=csv,noheader,nounits"])
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let ids: Vec<usize> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| l.trim().parse::<usize>().ok())
        .collect();
    if ids.is_empty() {
        None
    } else {
        Some(ids)
    }
}

/// Decide which GPU device ids this node should offer.
///
/// Resolution order:
/// 1. an explicit `--hydra-gpus N` override (honoring `CUDA_VISIBLE_DEVICES`
///    ordering if that env var is set, otherwise `0..N`);
/// 2. the `CUDA_VISIBLE_DEVICES` environment variable (set by SLURM `--gres`
///    allocations, or by the user);
/// 3. `nvidia-smi` enumeration;
/// 4. otherwise no GPUs.
pub(crate) fn detect(override_count: Option<usize>) -> Vec<usize> {
    let env_ids = std::env::var("CUDA_VISIBLE_DEVICES")
        .ok()
        .and_then(|v| parse_visible(&v));

    if let Some(n) = override_count {
        return match env_ids {
            Some(ids) if !ids.is_empty() => ids.into_iter().take(n).collect(),
            _ => (0..n).collect(),
        };
    }
    if let Some(ids) = env_ids {
        if !ids.is_empty() {
            return ids;
        }
    }
    nvidia_smi_indices().unwrap_or_default()
}
