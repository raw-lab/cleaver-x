//! Cross-node / multi-core execution backend over the vendored HydraMPP engine.
//! Compiled only with `--features hydra`.
//!
//! This mirrors [`crate::parallel::run`] at the call site, so the CLI swaps
//! backends with no other changes. Each input file becomes one HydraMPP task;
//! the `cleaver-core` engine runs unchanged on whichever worker node the
//! scheduler places it on. Because a scheduler cannot ship a closure to a
//! remote process, the task is a *named* function (registered by name, see
//! [`crate::jobs`]) operating on `serde`-serializable jobs.

#![cfg(feature = "hydra")]

use anyhow::{anyhow, Context, Result};
use serde::{de::DeserializeOwned, Serialize};

use hydra_mpp::prelude::*;

/// Build a HydraMPP [`Config`] from the CLI's cluster flags:
///   * `--head`            run as the head node (workers connect in)
///   * `--client ADDR`     join a head as a worker/client
///   * neither             local mode (single box); `--sim-gpus N` fakes N GPUs
///     so device scheduling can be exercised without real hardware.
pub fn build_config(
    head: bool,
    client: Option<&str>,
    cpus: Option<usize>,
    sim_gpus: Option<usize>,
    port: Option<u16>,
    quiet: bool,
) -> Config {
    let mut cfg = match client {
        Some(addr) => Config::client(addr.to_string()),
        None if head => Config::host(),
        None => Config::local(),
    };
    if let Some(c) = cpus {
        cfg = cfg.num_cpus(c);
    }
    if let Some(g) = sim_gpus {
        cfg = cfg.num_gpus(g);
    }
    if let Some(p) = port {
        cfg = cfg.port(p);
    }
    cfg.quiet(quiet)
}

/// Run one task per job on a HydraMPP cluster, returning results in input order.
///
/// `name` is the registered function name; `task` is the function compiled into
/// this binary (so every worker has it). `cpus`/`gpus` are the per-task resource
/// requests the scheduler honours — request 1 GPU for `stats` to have a device
/// pinned for the task; the request is clamped to what the cluster advertises so
/// a CPU-only box (or one with fewer GPUs) does not fail fast.
pub fn run<J, O>(
    name: &str,
    task: fn(J) -> O,
    jobs: Vec<J>,
    cpus: usize,
    gpus: usize,
    cfg: Config,
) -> Result<Vec<O>>
where
    J: Serialize + DeserializeOwned + Send + 'static,
    O: Serialize + DeserializeOwned + Send + 'static,
{
    let h = Hydra::new();
    h.register(name, task);
    h.init(cfg).context("starting HydraMPP")?;

    // Clamp the GPU request to the cluster's capacity so scheduling is
    // demonstrated rather than rejected.
    let avail = h
        .nodes()
        .ok()
        .and_then(|ns| ns.iter().map(|n| n.total_gpus).max())
        .unwrap_or(0);
    let g = gpus.min(avail);

    let mut ids = Vec::with_capacity(jobs.len());
    for j in &jobs {
        let id = h
            .task(name)
            .cpus(cpus.max(1))
            .gpus(g)
            .submit(j)
            .context("submitting task")?;
        ids.push(id);
    }

    // Block until every task is finished, then collect in submission order.
    let _ = h.wait(&ids, None, ids.len());
    let mut out = Vec::with_capacity(ids.len());
    for id in ids {
        let res = h.get(id).context("fetching result")?;
        if let Some(e) = res.error() {
            let msg = format!("task '{}' failed: {e}", res.func_name);
            h.shutdown();
            return Err(anyhow!(msg));
        }
        out.push(res.value::<O>().context("decoding task result")?);
    }
    h.shutdown();
    Ok(out)
}
