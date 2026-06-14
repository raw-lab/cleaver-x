//! Compute nodes.
//!
//! `nodes[0]` is always the local node. Remote nodes carry a TCP link the
//! scheduler uses to ship tasks; the matching reader thread feeds results back.
//!
//! Each node tracks two reservable resources: **CPU slots** (a free count) and
//! **GPU devices** (concrete device ids, each free or busy). Tracking ids — not
//! just a count — lets the scheduler pin specific devices to a job so the task
//! can set `CUDA_VISIBLE_DEVICES` correctly.

use std::net::TcpStream;
use std::sync::{Arc, Mutex};

/// How the scheduler reaches a node.
pub(crate) enum NodeLink {
    /// The driver's own machine — work runs in local threads.
    Local,
    /// A connected client — work is shipped over this TCP stream.
    ///
    /// The stream is wrapped in a `Mutex` because the scheduler thread writes
    /// `Task` frames while a dedicated reader thread (holding a `try_clone`)
    /// reads `Result` frames.
    Remote { write: Arc<Mutex<TcpStream>> },
}

/// A compute node in the cluster.
pub(crate) struct Node {
    pub hostname: String,
    pub address: String,
    pub total_cpus: usize,
    pub free_cpus: usize,
    /// Device ids this node owns (e.g. `[0, 1, 2, 3]`).
    pub gpu_ids: Vec<usize>,
    /// Parallel to `gpu_ids`: `true` if that device is currently free.
    pub gpu_free: Vec<bool>,
    pub link: NodeLink,
}

impl Node {
    pub(crate) fn local(hostname: String, cpus: usize, gpu_ids: Vec<usize>) -> Self {
        let gpu_free = vec![true; gpu_ids.len()];
        Node {
            hostname,
            address: "local".into(),
            total_cpus: cpus,
            free_cpus: cpus,
            gpu_ids,
            gpu_free,
            link: NodeLink::Local,
        }
    }

    pub(crate) fn remote(
        hostname: String,
        address: String,
        cpus: usize,
        gpu_ids: Vec<usize>,
        write: Arc<Mutex<TcpStream>>,
    ) -> Self {
        let gpu_free = vec![true; gpu_ids.len()];
        Node {
            hostname,
            address,
            total_cpus: cpus,
            free_cpus: cpus,
            gpu_ids,
            gpu_free,
            link: NodeLink::Remote { write },
        }
    }

    /// Total number of GPUs this node owns.
    pub(crate) fn total_gpus(&self) -> usize {
        self.gpu_ids.len()
    }

    /// Number of GPUs currently free.
    pub(crate) fn free_gpus(&self) -> usize {
        self.gpu_free.iter().filter(|&&f| f).count()
    }

    /// Reserve `k` free GPUs, returning their device ids (or `None` if fewer
    /// than `k` are free). Marks the chosen devices busy.
    pub(crate) fn reserve_gpus(&mut self, k: usize) -> Option<Vec<usize>> {
        if k == 0 {
            return Some(Vec::new());
        }
        if self.free_gpus() < k {
            return None;
        }
        let mut chosen = Vec::with_capacity(k);
        for (free, &id) in self.gpu_free.iter_mut().zip(self.gpu_ids.iter()) {
            if *free {
                *free = false;
                chosen.push(id);
                if chosen.len() == k {
                    break;
                }
            }
        }
        Some(chosen)
    }

    /// Return the given device ids to the free pool.
    pub(crate) fn release_gpus(&mut self, ids: &[usize]) {
        for id in ids {
            if let Some(pos) = self.gpu_ids.iter().position(|g| g == id) {
                self.gpu_free[pos] = true;
            }
        }
    }
}

/// Public, cloneable snapshot of a node (returned by `Hydra::nodes()`).
#[derive(Debug, Clone)]
pub struct NodeInfo {
    pub hostname: String,
    pub address: String,
    pub total_cpus: usize,
    pub free_cpus: usize,
    pub total_gpus: usize,
    pub free_gpus: usize,
    /// The device ids this node owns.
    pub gpu_ids: Vec<usize>,
    pub is_local: bool,
}
