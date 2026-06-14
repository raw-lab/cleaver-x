//! Error type for the HydraMPP engine.
//!
//! Hand-written (no `thiserror`) to keep the dependency tree minimal — the crate
//! pulls in nothing just to format these messages.

use std::fmt;

/// Result alias used throughout the crate's public API.
pub type Result<T> = std::result::Result<T, HydraError>;

/// Everything that can go wrong inside HydraMPP.
#[derive(Debug)]
pub enum HydraError {
    /// `init()` was called twice without a `shutdown()` in between.
    AlreadyRunning,
    /// An operation needed a running runtime but none was started.
    NotRunning,
    /// A task name was submitted that was never registered.
    UnknownTask(String),
    /// A job id was requested that the queue does not know about.
    UnknownJob(u64),
    /// Serialization/deserialization of task arguments or results failed.
    Serde(String),
    /// The remote task ran but returned an error / panicked.
    TaskFailed { func: String, message: String },
    /// A networking / transport level failure.
    Net(std::io::Error),
    /// A handshake or wire-protocol violation.
    Protocol(String),
    /// SLURM auto-configuration failed.
    Slurm(String),
}

impl fmt::Display for HydraError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HydraError::AlreadyRunning => f.write_str(
                "HydraMPP is already running; call shutdown() before re-initializing",
            ),
            HydraError::NotRunning => f.write_str("HydraMPP is not running; call init() first"),
            HydraError::UnknownTask(t) => write!(
                f,
                "task '{t}' is not registered; register it before calling remote()"
            ),
            HydraError::UnknownJob(id) => write!(f, "unknown job id: {id}"),
            HydraError::Serde(m) => write!(f, "serialization error: {m}"),
            HydraError::TaskFailed { func, message } => {
                write!(f, "remote task '{func}' failed: {message}")
            }
            HydraError::Net(e) => write!(f, "network error: {e}"),
            HydraError::Protocol(m) => write!(f, "protocol error: {m}"),
            HydraError::Slurm(m) => write!(f, "SLURM configuration error: {m}"),
        }
    }
}

impl std::error::Error for HydraError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            HydraError::Net(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for HydraError {
    fn from(e: std::io::Error) -> Self {
        HydraError::Net(e)
    }
}

impl HydraError {
    pub(crate) fn serde<E: fmt::Display>(e: E) -> Self {
        HydraError::Serde(e.to_string())
    }
}
