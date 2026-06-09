//! The task registry.
//!
//! HydraMPP dispatches work *by name*: every node registers the same functions
//! under the same names (the Python original keeps a global `WORKERS` dict). We
//! reproduce that model in a type-safe way. A registered handler
//! `Fn(In) -> Out` is wrapped into a type-erased
//! `Fn(&[u8]) -> Outcome` that:
//!   1. deserializes the argument bytes into `In`,
//!   2. runs the handler inside `catch_unwind` (so a panicking task is isolated,
//!      exactly like the Python `try/except` around each call),
//!   3. serializes `Out` back into bytes.

use std::any::Any;
use std::collections::HashMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, RwLock};

use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::wire::Outcome;

/// A type-erased handler over serialized argument bytes (the wire path).
pub(crate) type ErasedTask = Arc<dyn Fn(&[u8]) -> Outcome + Send + Sync>;

/// A type-erased handler over a typed, boxed argument (the local zero-copy
/// path): it skips argument *de*serialization and returns the result as a typed
/// in-memory value — no serialization on either side. The result is downcast
/// back by `get_typed`. This is what makes `map_owned` end-to-end zero-copy.
pub(crate) type ErasedTaskTyped =
    Arc<dyn Fn(Box<dyn Any + Send>) -> std::result::Result<Box<dyn Any + Send>, String> + Send + Sync>;

/// One registered task: both dispatch paths share the same handler `Arc`.
#[derive(Clone)]
struct Handlers {
    bytes: ErasedTask,
    typed: ErasedTaskTyped,
}

/// Thread-safe map of task name → handlers. Cloning is cheap (`Arc`).
#[derive(Clone, Default)]
pub(crate) struct Registry {
    inner: Arc<RwLock<HashMap<Arc<str>, Handlers>>>,
}

impl Registry {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Register a typed handler under `name`.
    pub(crate) fn register<In, Out, F>(&self, name: impl Into<String>, f: F)
    where
        In: DeserializeOwned + Send + 'static,
        Out: Serialize + Send + 'static,
        F: Fn(In) -> Out + Send + Sync + 'static,
    {
        let f = Arc::new(f);
        let f_bytes = f.clone();
        let bytes: ErasedTask = Arc::new(move |b: &[u8]| {
            let arg: In = match bincode::deserialize(b) {
                Ok(a) => a,
                Err(e) => return Outcome::Err(format!("argument decode failed: {e}")),
            };
            // Isolate panics so one bad task cannot kill the worker thread.
            let f = f_bytes.clone();
            let result = catch_unwind(AssertUnwindSafe(|| f(arg)));
            match result {
                Ok(out) => match bincode::serialize(&out) {
                    Ok(b) => Outcome::Ok(b),
                    Err(e) => Outcome::Err(format!("result encode failed: {e}")),
                },
                Err(payload) => Outcome::Err(format!("task panicked: {}", panic_msg(payload))),
            }
        });
        // Zero-copy local path: run on the in-memory value and return the
        // result as a typed value (no argument decode, no result encode).
        let f_typed = f.clone();
        let typed: ErasedTaskTyped = Arc::new(move |boxed: Box<dyn Any + Send>| {
            let arg: In = match boxed.downcast::<In>() {
                Ok(a) => *a,
                Err(_) => return Err("internal: typed argument type mismatch".into()),
            };
            let f = f_typed.clone();
            match catch_unwind(AssertUnwindSafe(|| f(arg))) {
                Ok(out) => Ok(Box::new(out) as Box<dyn Any + Send>),
                Err(payload) => Err(format!("task panicked: {}", panic_msg(payload))),
            }
        });
        self.inner
            .write()
            .unwrap()
            .insert(Arc::from(name.into()), Handlers { bytes, typed });
    }

    /// Return the registry's own interned `Arc<str>` for `name`, if registered.
    /// Submitters clone this instead of allocating a fresh `String` per job, so a
    /// batch of N jobs costs one lookup and N refcount bumps, not N allocations.
    pub(crate) fn intern(&self, name: &str) -> Option<Arc<str>> {
        self.inner
            .read()
            .unwrap()
            .get_key_value(name)
            .map(|(k, _)| k.clone())
    }

    /// Look up the bytes handler by name (cloned `Arc`, cheap).
    pub(crate) fn get(&self, name: &str) -> Option<ErasedTask> {
        self.inner.read().unwrap().get(name).map(|h| h.bytes.clone())
    }

    /// Look up the typed (zero-copy) handler by name (cloned `Arc`, cheap).
    pub(crate) fn get_typed(&self, name: &str) -> Option<ErasedTaskTyped> {
        self.inner.read().unwrap().get(name).map(|h| h.typed.clone())
    }
}

/// Best-effort extraction of a panic message from the unwind payload.
fn panic_msg(payload: Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown panic".to_string()
    }
}
