//! In-node parallel driver. Splitting N files is embarrassingly parallel, so
//! each file is an independent task on a work-stealing pool. `threads = 0`
//! uses all logical cores. Results are returned in input order.
//!
//! This is the local backend; the cross-node backend (HydraMPP) lives in
//! `hydra.rs` behind `--features hydra` and mirrors this signature.

use anyhow::{Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

pub fn run<I, R, F>(items: &[I], threads: usize, f: F) -> Result<Vec<R>>
where
    I: Sync,
    R: Send,
    F: Fn(&I) -> R + Sync + Send,
{
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("building worker pool")?;
    Ok(pool.install(|| items.par_iter().map(&f).collect()))
}
