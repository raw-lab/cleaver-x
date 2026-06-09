//! Work units shared by both execution backends.
//!
//! A distributed scheduler cannot ship an arbitrary Rust closure to a remote
//! worker — the work has to be a *named function over serializable data* that is
//! already compiled into the worker binary. So every unit of work here is a
//! plain `fn(Job) -> Outcome` plus `serde` structs, and the **same** function is
//! driven by the in-node rayon backend ([`crate::parallel`]) and the cross-node
//! HydraMPP backend ([`crate::hydra`]). That is what makes `--features hydra` a
//! drop-in: the `cleaver-core` engine runs unchanged whether a task lands on
//! this core or on a remote node.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use cleaver_core::{align, chunk_file, formats, genome, Config, Format, Match, Mode};

use crate::gpu;

/// `Format` lives in core and is intentionally `serde`-free, so jobs carry the
/// forced format as a short string and resolve it on the worker.
pub fn fmt_from_str(s: &str) -> Option<Format> {
    match s.to_ascii_lowercase().as_str() {
        "fasta" | "fa" | "fna" | "ffn" | "faa" => Some(Format::Fasta),
        "fastq" | "fq" => Some(Format::Fastq),
        "sam" => Some(Format::Sam),
        "bam" => Some(Format::Bam),
        _ => None,
    }
}

pub fn fmt_name(f: Format) -> &'static str {
    match f {
        Format::Fasta => "fasta",
        Format::Fastq => "fastq",
        Format::Sam => "sam",
        Format::Bam => "bam",
    }
}

// ----------------------------------------------------------------------------
// split
// ----------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplitJob {
    pub input: PathBuf,
    pub outdir: PathBuf,
    /// Max bytes per chunk for sequence formats (FASTA/FASTQ).
    pub chunk_size: u64,
    /// Records per chunk for alignment formats (SAM/BAM).
    pub records: u64,
    /// Forced format name, or `None` to auto-detect.
    pub format: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplitOutcome {
    pub input: PathBuf,
    pub ok: bool,
    pub format: String,
    pub chunks: usize,
    /// Empty when `ok`.
    pub error: String,
}

/// Split one input into record-aware chunks. Pure function of its `SplitJob`,
/// so it is safe to register with HydraMPP and run on any worker.
pub fn do_split(job: SplitJob) -> SplitOutcome {
    let run = || -> anyhow::Result<(Format, usize)> {
        let fmt = match job.format.as_deref().and_then(fmt_from_str) {
            Some(f) => f,
            None => formats::detect(&job.input)?,
        };
        let n = match fmt {
            Format::Fasta => {
                let cfg = Config {
                    chunk_size: job.chunk_size,
                    mode: Mode::Delimiter { delim: b">".to_vec(), matching: Match::StartsWith },
                };
                chunk_file(&job.input, &job.outdir, &cfg)?.chunks.len()
            }
            Format::Fastq => {
                let cfg = Config { chunk_size: job.chunk_size, mode: Mode::Lines { n: 4 } };
                chunk_file(&job.input, &job.outdir, &cfg)?.chunks.len()
            }
            Format::Sam => align::split_sam(&job.input, &job.outdir, job.records)?.len(),
            Format::Bam => align::split_bam(&job.input, &job.outdir, job.records)?.len(),
        };
        Ok((fmt, n))
    };
    match run() {
        Ok((fmt, chunks)) => SplitOutcome {
            input: job.input,
            ok: true,
            format: fmt_name(fmt).to_string(),
            chunks,
            error: String::new(),
        },
        Err(e) => SplitOutcome {
            input: job.input,
            ok: false,
            format: String::new(),
            chunks: 0,
            error: format!("{e:#}"),
        },
    }
}

// ----------------------------------------------------------------------------
// stats (the GPU-parallel workload)
// ----------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatsJob {
    pub input: PathBuf,
    pub format: Option<String>,
}

/// Full genome statistics for one file (flat scalars so it ships over HydraMPP).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatsOutcome {
    pub input: PathBuf,
    pub ok: bool,
    pub n_seqs: u64,
    pub total_bp: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub mean_len: f64,
    pub median_len: u64,
    pub n50: u64,
    pub l50: u64,
    pub n90: u64,
    pub l90: u64,
    pub gc_percent: f64,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub other: u64,
    /// Which backend tallied base composition: `cpu` or `gpu:<id>`.
    pub device: String,
    pub error: String,
}

/// The device this task was pinned to. Under HydraMPP the scheduler sets a
/// thread-local set of GPU ids before invoking the task; without HydraMPP we
/// fall back to device 0 when the GPU kernel is compiled in, else CPU.
fn pinned_device() -> Option<usize> {
    #[cfg(feature = "hydra")]
    if let Some(d) = hydra_mpp_core::current_gpus().first().copied() {
        return Some(d);
    }
    #[cfg(feature = "gpu")]
    {
        return Some(0);
    }
    #[allow(unreachable_code)]
    None
}

/// Compute full genome statistics for one file. Base composition is tallied on
/// the GPU when the kernel is compiled in and a device is pinned (the same
/// reduction as `compute::count_bytes`), and reused by the length pass;
/// otherwise lengths and composition are tallied together in one CPU pass. The
/// N/L assembly metrics come from `rustyomestats::stats::compute_nl`.
pub fn do_stats(job: StatsJob) -> StatsOutcome {
    let device = pinned_device();
    let run = || -> anyhow::Result<(genome::GenomeStats, String)> {
        let fmt = match job.format.as_deref().and_then(fmt_from_str) {
            Some(f) => f,
            None => formats::detect(&job.input)?,
        };
        if matches!(fmt, Format::Sam | Format::Bam) {
            anyhow::bail!("stats expects FASTA/FASTQ (got {})", fmt_name(fmt));
        }

        // GPU base-composition path when compiled in and a device is pinned;
        // its counts are then reused so the length pass skips re-counting.
        let (precomputed, backend) = if gpu::compute_enabled() && device.is_some() {
            match gpu::count_file_accel(&job.input, fmt, device) {
                Ok((c, b)) => (Some(c), b),
                Err(e) => {
                    eprintln!("  warning: GPU count failed ({e:#}); using CPU");
                    (None, "cpu")
                }
            }
        } else {
            (None, "cpu")
        };

        let gs = genome::genome_stats_file(&job.input, fmt, precomputed)?;
        let dev = match (backend, device) {
            ("gpu", Some(d)) => format!("gpu:{d}"),
            // Scheduler pinned a device but the CUDA kernel isn't compiled in
            // (built without --features gpu): work ran on CPU, device shown.
            ("cpu", Some(d)) => format!("cpu(gpu:{d})"),
            _ => "cpu".to_string(),
        };
        Ok((gs, dev))
    };
    match run() {
        Ok((s, device)) => StatsOutcome {
            input: job.input,
            ok: true,
            n_seqs: s.n_seqs,
            total_bp: s.total_bp,
            min_len: s.min_len,
            max_len: s.max_len,
            mean_len: s.mean_len,
            median_len: s.median_len,
            n50: s.n50(),
            l50: s.l50(),
            n90: s.n90(),
            l90: s.l90(),
            gc_percent: s.gc_percent,
            a: s.counts.a,
            c: s.counts.c,
            g: s.counts.g,
            t: s.counts.t,
            n: s.counts.n,
            other: s.counts.other,
            device,
            error: String::new(),
        },
        Err(e) => StatsOutcome {
            input: job.input,
            ok: false,
            n_seqs: 0,
            total_bp: 0,
            min_len: 0,
            max_len: 0,
            mean_len: 0.0,
            median_len: 0,
            n50: 0,
            l50: 0,
            n90: 0,
            l90: 0,
            gc_percent: 0.0,
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
            other: 0,
            device: "cpu".to_string(),
            error: format!("{e:#}"),
        },
    }
}
