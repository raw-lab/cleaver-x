//! SLURM auto-wiring.
//!
//! Faithful to the Python original's `slurm()`: when the head node is started
//! with `--hydra-slurm <nodelist>`, HydraMPP
//!   1. expands the nodelist with `scontrol show hostnames <nodelist>`,
//!   2. resolves the head IP with `srun … -w <node0> hostname --ip-address`,
//!   3. `srun`s a copy of *this* program on every other node, replaying the
//!      original arguments and appending `--hydra-client <head_ip>`,
//!   4. becomes the host itself.
//!
//! Child stdout/stderr are redirected to `tmp-hydra/<date>/<node>.std{out,err}`
//! exactly as the Python version does, so a misbehaving client can be inspected.

use std::fs::{self, File};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::error::{HydraError, Result};
use crate::printlog;

/// Handles to the `srun` client processes spawned by the head node. Dropping
/// this struct does *not* kill the clients (they live for the job); it simply
/// stops us from tracking them. Kept so callers may `wait`/inspect if desired.
#[derive(Default)]
pub struct SlurmClients {
    pub children: Vec<Child>,
}

/// Directory `tmp-hydra/<YYYY-MM-DD_HH.MM>` used for per-node client logs,
/// mirroring the Python `TMP_HYDRA` path.
fn tmp_hydra_dir() -> PathBuf {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Cheap civil-time breakdown (UTC) — good enough for a log folder name.
    let days = secs / 86_400;
    let day = secs % 86_400;
    let (h, mi) = (day / 3600, (day % 3600) / 60);
    // Convert epoch-days to Y-M-D (proleptic Gregorian).
    let (y, m, d) = civil_from_days(days as i64);
    PathBuf::from("tmp-hydra").join(format!("{y:04}-{m:02}-{d:02}_{h:02}.{mi:02}"))
}

/// Days-since-epoch → (year, month, day). Howard Hinnant's algorithm.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// Expand a SLURM nodelist into individual hostnames via `scontrol`.
fn expand_nodelist(nodelist: &str) -> Result<Vec<String>> {
    let out = Command::new("scontrol")
        .args(["show", "hostnames", nodelist])
        .output()
        .map_err(|e| HydraError::Slurm(format!("failed to run scontrol: {e}")))?;
    if !out.status.success() {
        return Err(HydraError::Slurm(format!(
            "scontrol exited with {}",
            out.status
        )));
    }
    let nodes: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();
    if nodes.is_empty() {
        return Err(HydraError::Slurm("scontrol returned no hostnames".into()));
    }
    Ok(nodes)
}

/// Resolve the head node's IP with `srun … hostname --ip-address`.
fn resolve_head_ip(head_node: &str) -> Result<String> {
    let out = Command::new("srun")
        .args(["--nodes=1", "--ntasks=1", "-w", head_node, "hostname", "--ip-address"])
        .output()
        .map_err(|e| HydraError::Slurm(format!("failed to srun hostname: {e}")))?;
    let ip = String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string();
    if ip.is_empty() {
        return Err(HydraError::Slurm("could not determine head IP".into()));
    }
    Ok(ip)
}

/// Bring up SLURM clients from the head node.
///
/// `passthrough` is the original program's non-HydraMPP arguments (so the
/// clients run the same workload), `port`/`cpus` are forwarded to each client.
/// Returns once every client has been `srun`-launched; the head should then
/// enter host mode. The head node's hostname (nodes\[0]) is returned so the
/// caller can log it.
pub fn spawn_clients(
    nodelist: &str,
    passthrough: &[String],
    port: u16,
    cpus: Option<usize>,
    gpus: Option<usize>,
) -> Result<(String, SlurmClients)> {
    printlog!("SLURM: expanding nodelist '{nodelist}'");
    let nodes = expand_nodelist(nodelist)?;
    printlog!("SLURM: node list = {nodes:?}");

    let head_node = nodes[0].clone();
    let head_ip = resolve_head_ip(&head_node)?;
    printlog!("SLURM: head node {head_node} -> {head_ip}");

    let log_dir = tmp_hydra_dir();
    if nodes.len() > 1 {
        fs::create_dir_all(&log_dir)
            .map_err(|e| HydraError::Slurm(format!("cannot create {}: {e}", log_dir.display())))?;
    }

    // The executable to relaunch on every worker node (this very binary).
    let exe = std::env::current_exe()
        .map_err(|e| HydraError::Slurm(format!("cannot find current exe: {e}")))?;

    let mut clients = SlurmClients::default();
    for node in &nodes[1..] {
        let mut cmd = Command::new("srun");
        cmd.args(["--nodes=1", "--ntasks=1", "-w", node]);
        cmd.arg(&exe);
        cmd.args(passthrough);
        cmd.args(["--hydra-client", &head_ip]);
        cmd.args(["--hydra-port", &port.to_string()]);
        if let Some(n) = cpus {
            cmd.args(["--hydra-cpus", &n.to_string()]);
        }
        if let Some(n) = gpus {
            cmd.args(["--hydra-gpus", &n.to_string()]);
        }

        let stdout = File::create(log_dir.join(format!("{node}.stdout")))
            .map_err(HydraError::Net)?;
        let stderr = File::create(log_dir.join(format!("{node}.stderr")))
            .map_err(HydraError::Net)?;
        cmd.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

        match cmd.spawn() {
            Ok(child) => {
                printlog!("SLURM: launched client on {node}");
                clients.children.push(child);
            }
            Err(e) => printlog!("SLURM WARN: failed to launch client on {node}: {e}"),
        }
    }

    printlog!("SLURM: setting node {head_node} as host");
    Ok((head_node, clients))
}
