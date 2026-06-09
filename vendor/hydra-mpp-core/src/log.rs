//! Tiny dependency-free logger that mirrors the Python original's
//! `HYDRA HH:MM:SS [hostname]:` prefix. Lines go to stderr by default; if a log
//! directory has been registered they are also appended to `stdout.log`.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;
use std::time::{SystemTime, UNIX_EPOCH};

static QUIET: AtomicBool = AtomicBool::new(false);
static LOG_DIR: Mutex<Option<PathBuf>> = Mutex::new(None);

/// Silence (or un-silence) HydraMPP's own log chatter. Task output is untouched.
pub fn set_quiet(quiet: bool) {
    QUIET.store(quiet, Ordering::Relaxed);
}

/// Register a directory; future log lines are appended to `<dir>/stdout.log`.
#[allow(dead_code)] // public utility; not called internally
pub fn set_log_dir(dir: PathBuf) {
    let _ = std::fs::create_dir_all(&dir);
    *LOG_DIR.lock().unwrap() = Some(dir);
}

/// Current wall-clock time formatted as `HH:MM:SS` (UTC, no tz database needed).
pub fn hms() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let day = secs % 86_400;
    format!("{:02}:{:02}:{:02}", day / 3600, (day % 3600) / 60, day % 60)
}

/// The local hostname (best effort, falls back to "localhost").
pub fn hostname() -> String {
    gethostname::gethostname()
        .into_string()
        .unwrap_or_else(|_| "localhost".to_string())
}

#[doc(hidden)]
pub fn _log(msg: &str) {
    if QUIET.load(Ordering::Relaxed) {
        return;
    }
    let line = format!("HYDRA {} [{}]:\t{}", hms(), hostname(), msg);
    eprintln!("{line}");
    if let Some(dir) = LOG_DIR.lock().unwrap().as_ref() {
        if let Ok(mut f) = OpenOptions::new()
            .create(true)
            .append(true)
            .open(dir.join("stdout.log"))
        {
            let _ = writeln!(f, "{line}");
        }
    }
}

/// `printlog!`-style macro: same call sites as the Python `printlog(...)`.
#[macro_export]
macro_rules! printlog {
    ($($arg:tt)*) => {{
        $crate::log::_log(&format!($($arg)*));
    }};
}
