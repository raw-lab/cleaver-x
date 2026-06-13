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

use cleaver_core::{align, annotation, chunk_file, count, formats, genome, Config, Format, Match, Mode};

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
    pub gc_percent: f64,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub other: u64,
    /// Which backend produced the base composition: `cpu` or `gpu:<id>`.
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

/// Compute genome statistics for one file. Base composition runs on the GPU
/// when available (else CPU); the N50/length pass is CPU.
pub fn do_stats(job: StatsJob) -> StatsOutcome {
    let device = pinned_device();
    let run = || -> anyhow::Result<(genome::GenomeStats, String)> {
        let fmt = match job.format.as_deref().and_then(fmt_from_str) {
            Some(f) => f,
            None => formats::detect(&job.input)?,
        };
        if !fmt.is_sequence() {
            anyhow::bail!("stats expects FASTA/FASTQ input, got {}", fmt_name(fmt));
        }
        let (counts, backend) = gpu::count_file_accel(&job.input, fmt, device)?;
        let stats = genome::genome_stats_file(&job.input, fmt, Some(counts))?;
        let dev = match (backend, device) {
            ("gpu", Some(d)) => format!("gpu:{d}"),
            ("cpu", Some(d)) => format!("cpu(gpu:{d})"),
            _ => "cpu".to_string(),
        };
        Ok((stats, dev))
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

// ----------------------------------------------------------------------------
// count (featureCounts / htseq-style read counting)
// ----------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CountJob {
    pub input: PathBuf,
    pub annotation: PathBuf,
    pub feature_type: String,
    pub group_by: String,
    pub stranded: u8,
    pub min_mapq: u8,
    pub count_multimappers: bool,
    pub allow_multi_overlap: bool,
    /// 0 = union, 1 = intersection-strict, 2 = intersection-nonempty.
    pub mode: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CountOutcome {
    pub input: PathBuf,
    pub ok: bool,
    pub sample: String,
    /// Number of indexed feature lines (of the requested type).
    pub n_features: u64,
    /// Meta-feature ids in annotation first-seen order (shared across inputs).
    pub gene_ids: Vec<String>,
    pub counts: Vec<u64>,
    pub assigned: u64,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub unmapped: u64,
    pub low_mapq: u64,
    pub multimapping: u64,
    pub total: u64,
    pub error: String,
}

fn mode_from_u8(m: u8) -> annotation::Mode {
    match m {
        1 => annotation::Mode::IntersectionStrict,
        2 => annotation::Mode::IntersectionNonempty,
        _ => annotation::Mode::Union,
    }
}

/// Count one alignment file against its annotation. Pure over `CountJob`: the
/// annotation is (re)parsed here so the unit is self-contained on any worker.
pub fn do_count(job: CountJob) -> CountOutcome {
    let sample = job
        .input
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
        .to_string();
    let run = || -> anyhow::Result<(Vec<String>, u64, count::FileCounts)> {
        let fmt = formats::detect(&job.input)?;
        if !fmt.is_alignment() {
            anyhow::bail!("count expects SAM/BAM input, got {}", fmt_name(fmt));
        }
        let ann = annotation::Annotation::from_path(&job.annotation, &job.feature_type, &job.group_by)?;
        let params = count::CountParams {
            stranded: job.stranded,
            min_mapq: job.min_mapq,
            count_multimappers: job.count_multimappers,
            allow_multi_overlap: job.allow_multi_overlap,
            mode: mode_from_u8(job.mode),
            primary_only: true,
        };
        let fc = count::count_file(&job.input, fmt, &ann, &params)?;
        Ok((ann.gene_ids().to_vec(), ann.n_features, fc))
    };
    match run() {
        Ok((gene_ids, n_features, fc)) => CountOutcome {
            input: job.input,
            ok: true,
            sample,
            n_features,
            gene_ids,
            counts: fc.counts,
            assigned: fc.assigned,
            no_feature: fc.no_feature,
            ambiguous: fc.ambiguous,
            unmapped: fc.unmapped,
            low_mapq: fc.low_mapq,
            multimapping: fc.multimapping,
            total: fc.total,
            error: String::new(),
        },
        Err(e) => CountOutcome {
            input: job.input,
            ok: false,
            sample,
            n_features: 0,
            gene_ids: Vec::new(),
            counts: Vec::new(),
            assigned: 0,
            no_feature: 0,
            ambiguous: 0,
            unmapped: 0,
            low_mapq: 0,
            multimapping: 0,
            total: 0,
            error: format!("{e:#}"),
        },
    }
}
