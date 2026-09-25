//! Cleaver — a fast, scalable, record-aware splitter and converter for the
//! common sequence and alignment formats (FASTA/FASTQ/SAM/BAM).
//!
//! Work fans out one-file-per-task. The default build runs tasks on an in-node
//! work-stealing pool (rayon); `--features hydra` runs the very same tasks on a
//! HydraMPP cluster (multi-core or cross-node), with `stats` offloaded to a GPU
//! when one is scheduled.

mod doctor;
mod gpu;
mod jobs;

#[cfg(feature = "hydra")]
mod hydra;
#[cfg(not(feature = "hydra"))]
mod parallel;

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{ArgAction, Parser, Subcommand, ValueEnum};

use cleaver::{align, compute::Counts, parse_size, Format};

const BANNER: &str = r"
   ___ _
  / __\ | ___  __ ___   _____ _ __
 / /  | |/ _ \/ _` \ \ / / _ \ '__|
/ /___| |  __/ (_| |\ V /  __/ |
\____/|_|\___|\__,_| \_/ \___|_|   record-aware FASTA/FASTQ/SAM/BAM
";

/// Cluster controls for the HydraMPP backend (only present with `--features hydra`).
#[cfg(feature = "hydra")]
#[derive(clap::Args)]
struct Cluster {
    /// Run as a HydraMPP head node (workers connect to this process).
    #[arg(long, global = true)]
    head: bool,
    /// Join a HydraMPP head at ADDR (host or host:port) as a worker/client.
    #[arg(long, value_name = "ADDR", global = true)]
    client: Option<String>,
    /// CPUs to advertise to HydraMPP (default: all logical cores).
    #[arg(long, value_name = "N", global = true)]
    cpus: Option<usize>,
    /// GPUs to advertise/schedule; simulates N GPUs on a CPU-only box for testing.
    #[arg(long, value_name = "N", global = true)]
    sim_gpus: Option<usize>,
    /// HydraMPP TCP port (head/client).
    #[arg(long, value_name = "PORT", global = true)]
    port: Option<u16>,
}

#[derive(Parser)]
#[command(
    name = "cleaver",
    version,
    about = "Split, convert, profile & count FASTA/FASTQ/SAM/BAM without cutting records.",
    before_help = BANNER,
    disable_version_flag = true
)]
struct Cli {
    /// Print version and exit.
    #[arg(short = 'v', long = "version", action = ArgAction::Version)]
    version: Option<bool>,
    #[cfg(feature = "hydra")]
    #[command(flatten)]
    cluster: Cluster,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Split one or more files into record-aware chunks (format auto-detected).
    #[command(before_help = BANNER)]
    Split {
        /// Inputs: FASTA (.fasta/.fa/.fna/.ffn/.faa/.frn), FASTQ (.fastq/.fq),
        /// SAM (.sam), BAM (.bam). `.gz` is handled transparently.
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Output directory (chunks are named <stem>.NNNNN.<ext>).
        #[arg(short = 'o', long)]
        outdir: PathBuf,
        /// Max chunk size for sequence formats (FASTA/FASTQ): e.g. 1G, 256M, 50Mi.
        #[arg(short = 'c', long, default_value = "1G")]
        chunk_size: String,
        /// Records per chunk for alignment formats (SAM/BAM).
        #[arg(short = 'r', long, default_value_t = 1_000_000)]
        records: u64,
        /// Worker threads for the in-node backend. 0 = all logical cores.
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
        /// Force a format instead of auto-detecting.
        #[arg(long, value_enum)]
        format: Option<FmtArg>,
    },
    /// Convert between formats: SAM<->BAM (with mapped/unmapped), FASTQ->FASTA.
    #[command(before_help = BANNER)]
    Convert {
        /// Input file.
        input: PathBuf,
        /// Output file (extension selects the target format).
        output: PathBuf,
        /// Alignment only: write mapped and unmapped records to separate files
        /// (<out_stem>.mapped.<ext> and <out_stem>.unmapped.<ext>).
        #[arg(long, conflicts_with_all = ["mapped_only", "unmapped_only"])]
        partition: bool,
        /// Alignment only: keep mapped records only.
        #[arg(long, conflicts_with = "unmapped_only")]
        mapped_only: bool,
        /// Alignment only: keep unmapped records only.
        #[arg(long)]
        unmapped_only: bool,
    },
    /// Genome/assembly statistics (N50/L50/N90, lengths, GC%) for FASTA/FASTQ.
    #[command(before_help = BANNER)]
    Stats {
        /// Input FASTA/FASTQ files (`.gz` transparent).
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Force a format instead of auto-detecting.
        #[arg(long, value_enum)]
        format: Option<FmtArg>,
        /// Worker threads for the in-node backend. 0 = all logical cores.
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },
    /// Count reads per feature from BAM/SAM against a GTF/GFF (featureCounts/htseq-style).
    #[command(before_help = BANNER)]
    Count {
        /// One or more alignment files (BAM/SAM; `.gz` SAM transparent).
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// Annotation file (GTF or GFF3).
        #[arg(short = 'a', long)]
        annotation: PathBuf,
        /// Output count matrix (TSV). A `<output>.summary` is written alongside.
        #[arg(short = 'o', long)]
        output: PathBuf,
        /// Feature type(s) in column 3 to count. Multiple `;`-separated types
        /// (e.g. "exon;intron;xine") enable VERSE multi-feature counting.
        #[arg(short = 't', long, default_value = "exon")]
        feature_type: String,
        /// VERSE scheme when multiple feature types are given.
        #[arg(long, value_enum, default_value_t = SchemeArg::Independent)]
        scheme: SchemeArg,
        /// Attribute used to group features into meta-features (e.g. gene_id).
        #[arg(short = 'g', long, default_value = "gene_id")]
        group_by: String,
        /// Strandedness: 0 unstranded, 1 stranded, 2 reverse-stranded.
        #[arg(short = 's', long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=2))]
        stranded: u8,
        /// Minimum MAPQ to count an alignment (0 keeps all).
        #[arg(short = 'Q', long, default_value_t = 0)]
        min_mapq: u8,
        /// Count multi-mapping reads (NH > 1) instead of discarding them.
        #[arg(short = 'M', long)]
        count_multimappers: bool,
        /// Count primary alignments only (excludes FLAG 0x100), like
        /// featureCounts' --primary. Off by default, as in featureCounts.
        #[arg(long)]
        primary: bool,
        /// Count reads overlapping >1 meta-feature against all of them.
        #[arg(short = 'O', long)]
        allow_multi_overlap: bool,
        /// Overlap-resolution mode (htseq-style; alias for the common -z values).
        #[arg(long, value_enum, default_value_t = ModeArg::Union)]
        mode: ModeArg,
        /// VERSE assignment mode -z: 0 featureCounts, 1 union, 2 strict,
        /// 3 nonempty, 4 union-strict, 5 cover-length. Overrides --mode.
        #[arg(short = 'z', long, value_parser = clap::value_parser!(u8).range(0..=5))]
        assign_mode: Option<u8>,
        /// PE: only count fragments with both ends mapped.
        #[arg(short = 'B', long)]
        require_both_ends: bool,
        /// PE: exclude chimeric reads (supplementary / mate on another chr).
        #[arg(short = 'C', long)]
        exclude_chimeric: bool,
        /// PE: check fragment length is within [-d, -D].
        #[arg(short = 'P', long)]
        check_pe_dist: bool,
        /// PE: minimum fragment length (with -P).
        #[arg(short = 'd', long, default_value_t = 50)]
        min_frag_len: u32,
        /// PE: maximum fragment length (with -P).
        #[arg(short = 'D', long, default_value_t = 600)]
        max_frag_len: u32,
        /// Worker threads (one file per task) for the in-node backend. 0 = all cores.
        #[arg(short = 'T', long, default_value_t = 0)]
        threads: usize,
    },
    /// Preprocess FASTQ (a pure-Rust fastp): adapter/quality/polyX trimming + filtering.
    #[command(before_help = BANNER)]
    Fastp {
        /// Read1 input FASTQ (`.gz` transparent).
        #[arg(short = 'i', long = "in1")]
        in1: PathBuf,
        /// Read1 output FASTQ (`.gz` to compress).
        #[arg(short = 'o', long = "out1")]
        out1: PathBuf,
        /// Read2 input (paired-end).
        #[arg(short = 'I', long = "in2")]
        in2: Option<PathBuf>,
        /// Read2 output (paired-end).
        #[arg(short = 'O', long = "out2")]
        out2: Option<PathBuf>,
        /// JSON report path.
        #[arg(short = 'j', long, default_value = "fastp.json")]
        json: PathBuf,

        // ---- adapter trimming ----
        /// Disable adapter trimming.
        #[arg(short = 'A', long)]
        disable_adapter_trimming: bool,
        /// Adapter sequence for read1.
        #[arg(short = 'a', long)]
        adapter_sequence: Option<String>,
        /// Adapter sequence for read2 (PE).
        #[arg(long)]
        adapter_sequence_r2: Option<String>,
        /// FASTA of adapters to trim from both reads.
        #[arg(long)]
        adapter_fasta: Option<PathBuf>,
        /// PE: detect adapters by overlap analysis.
        #[arg(short = '2', long)]
        detect_adapter_for_pe: bool,

        // ---- fixed/global trimming ----
        /// Trim N bases from read1 5' end.
        #[arg(short = 'f', long, default_value_t = 0)]
        trim_front1: usize,
        /// Trim N bases from read1 3' end.
        #[arg(short = 't', long, default_value_t = 0)]
        trim_tail1: usize,
        /// Cap read1 length (keep first N bases; 0 = no cap).
        #[arg(short = 'b', long, default_value_t = 0)]
        max_len1: usize,
        /// Trim N bases from read2 5' end.
        #[arg(short = 'F', long, default_value_t = 0)]
        trim_front2: usize,
        /// Trim N bases from read2 3' end.
        #[arg(short = 'T', long, default_value_t = 0)]
        trim_tail2: usize,
        /// Cap read2 length.
        #[arg(short = 'B', long, default_value_t = 0)]
        max_len2: usize,

        // ---- polyG / polyX ----
        /// Trim polyG tails (common on NovaSeq/NextSeq).
        #[arg(short = 'g', long)]
        trim_poly_g: bool,
        /// Minimum length to detect polyG.
        #[arg(long, default_value_t = 10)]
        poly_g_min_len: usize,
        /// Trim polyX tails.
        #[arg(short = 'x', long)]
        trim_poly_x: bool,
        /// Minimum length to detect polyX.
        #[arg(long, default_value_t = 10)]
        poly_x_min_len: usize,

        // ---- sliding-window quality cutting ----
        /// Sliding-window quality cut from the 5' end.
        #[arg(short = '5', long)]
        cut_front: bool,
        /// Sliding-window quality cut from the 3' end.
        #[arg(short = '3', long)]
        cut_tail: bool,
        /// Sliding-window quality cut, dropping the window and all 3' of it.
        #[arg(short = 'r', long)]
        cut_right: bool,
        /// Window size for quality cutting.
        #[arg(short = 'W', long, default_value_t = 4)]
        cut_window_size: usize,
        /// Mean-quality threshold for quality cutting.
        #[arg(short = 'M', long, default_value_t = 20)]
        cut_mean_quality: u8,

        // ---- quality filtering ----
        /// Disable quality filtering.
        #[arg(short = 'Q', long)]
        disable_quality_filtering: bool,
        /// A base is "qualified" at or above this phred.
        #[arg(short = 'q', long, default_value_t = 15)]
        qualified_quality_phred: u8,
        /// Discard a read if more than this percent of bases are unqualified.
        #[arg(short = 'u', long, default_value_t = 40.0)]
        unqualified_percent_limit: f64,
        /// Discard a read with more than this many N bases.
        #[arg(short = 'n', long, default_value_t = 5)]
        n_base_limit: usize,
        /// Discard a read whose mean quality is below this (0 = off).
        #[arg(short = 'e', long, default_value_t = 0)]
        average_qual: u8,

        // ---- length / complexity filtering ----
        /// Disable length filtering.
        #[arg(short = 'L', long)]
        disable_length_filtering: bool,
        /// Discard reads shorter than this.
        #[arg(short = 'l', long, default_value_t = 15)]
        length_required: usize,
        /// Discard reads longer than this (0 = off).
        #[arg(long, default_value_t = 0)]
        length_limit: usize,
        /// Enable the low-complexity filter.
        #[arg(short = 'y', long)]
        low_complexity_filter: bool,
        /// Complexity threshold percent (0..100).
        #[arg(short = 'Y', long, default_value_t = 30.0)]
        complexity_threshold: f64,

        // ---- PE overlap ----
        /// PE: correct mismatched bases in the overlapped region.
        #[arg(short = 'c', long)]
        correction: bool,
        /// Minimum overlap length for PE overlap analysis.
        #[arg(long, default_value_t = 30)]
        overlap_len_require: usize,
        /// Max mismatched bases in a PE overlap.
        #[arg(long, default_value_t = 5)]
        overlap_diff_limit: usize,
        /// Max percent mismatched bases in a PE overlap.
        #[arg(long, default_value_t = 20.0)]
        overlap_diff_percent_limit: f64,

        // ---- misc ----
        /// Input qualities are phred+64.
        #[arg(short = '6', long)]
        phred64: bool,
        /// Stop after this many reads/pairs (0 = all).
        #[arg(long, default_value_t = 0)]
        reads_to_process: u64,
    },

    /// Demultiplex ONT reads by barcode (pure-Rust barbell): fit-align barcodes
    /// at both ends, assign, trim, and split into per-barcode FASTQ.
    #[command(before_help = BANNER)]
    Demux {
        /// Input reads (FASTQ; `.gz` transparent).
        #[arg(short = 'i', long)]
        input: PathBuf,
        /// Output folder for per-barcode FASTQ files.
        #[arg(short = 'o', long)]
        output: PathBuf,
        /// Barcode FASTA (barbell example barcode sets work here).
        #[arg(short = 'q', long)]
        queries: PathBuf,
        /// Bases at each read end to scan for a barcode.
        #[arg(long, default_value_t = 150)]
        window: usize,
        /// Minimum identity (0..1) to accept a barcode match.
        #[arg(long = "min-score", default_value_t = 0.75)]
        min_score: f64,
        /// Minimum identity margin of the best over the second-best barcode.
        #[arg(long = "min-score-diff", default_value_t = 0.05)]
        min_score_diff: f64,
        /// Hard cap on edit distance (overrides --min-score when set).
        #[arg(long = "max-errors")]
        max_errors: Option<usize>,
        /// Annotate/split only; do not trim the barcode from reads.
        #[arg(long = "no-trim")]
        no_trim: bool,
        /// Only scan the 5' end (default scans both ends).
        #[arg(long = "one-end")]
        one_end: bool,
    },

    /// samtools-compatible SAM/BAM utilities: view, sort, index, fastq,
    /// fasta, flagstat, idxstats, merge, faidx, depth.
    #[command(before_help = BANNER)]
    Samtools {
        #[command(subcommand)]
        cmd: SamCmd,
    },

    /// Self-test: verify the install and that every code path works.
    #[command(before_help = BANNER)]
    Doctor,
    /// Show version, backends, detected GPUs, and supported formats.
    #[command(before_help = BANNER)]
    Info,
    /// Show the banner, version, and build configuration.
    #[command(before_help = BANNER)]
    Version,
}

/// Parse a SAM FLAG value in decimal or `0x` hexadecimal (as samtools accepts).
fn parse_flag(s: &str) -> std::result::Result<u16, String> {
    let s = s.trim();
    let r = if let Some(h) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u16::from_str_radix(h, 16)
    } else {
        s.parse::<u16>()
    };
    r.map_err(|_| format!("invalid FLAG '{s}' (use a decimal or 0x-prefixed value)"))
}

#[derive(Subcommand)]
enum SamCmd {
    /// Filter/print alignments by flag/MAPQ/region/subsample (`samtools view`).
    #[command(before_help = BANNER)]
    View {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Write BAM instead of SAM.
        #[arg(short = 'b', long)]
        bam: bool,
        /// Print only the count of matching records.
        #[arg(short = 'c', long)]
        count: bool,
        /// Write the header only.
        #[arg(short = 'H', long = "header-only")]
        header_only: bool,
        /// Only output records with ALL of these flag bits set (decimal or 0x).
        #[arg(short = 'f', value_parser = parse_flag, default_value_t = 0)]
        require: u16,
        /// Do not output records with ANY of these flag bits set (decimal or 0x).
        #[arg(short = 'F', value_parser = parse_flag, default_value_t = 0)]
        exclude: u16,
        /// Skip records with mapping quality below this value.
        #[arg(short = 'q', default_value_t = 0)]
        min_mapq: u8,
        /// Region filter: `chr` or `chr:beg-end` (linear scan; no index needed).
        #[arg(short = 'r', long)]
        region: Option<String>,
        /// Subsample: FLOAT where the integer part is the seed and the fraction
        /// is the fraction of records to keep (e.g. `42.1` keeps ~10%).
        #[arg(short = 's', long)]
        subsample: Option<f64>,
    },

    /// Sort by coordinate (default) or read name with `-n` (`samtools sort`).
    #[command(before_help = BANNER)]
    Sort {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Sort by read name instead of coordinate.
        #[arg(short = 'n', long = "by-name")]
        by_name: bool,
        /// Output format: `bam` (default) or `sam`.
        #[arg(short = 'O', long, default_value = "bam")]
        format: String,
        /// Thread count (accepted for compatibility; this build is single-core).
        #[arg(short = '@', long = "threads", default_value_t = 1)]
        threads: usize,
    },

    /// Build a BAI index for a coordinate-sorted BAM (`samtools index`).
    #[command(before_help = BANNER)]
    Index {
        /// Coordinate-sorted BAM.
        input: PathBuf,
        /// Output index path (default: `<input>.bai`).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },

    /// Extract reads to FASTQ, flag-filtered (`samtools fastq`; `-f 4` = unmapped).
    #[command(before_help = BANNER)]
    Fastq {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Single/interleaved output (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Require ALL of these flag bits (decimal or 0x). `-f 4` keeps unmapped.
        #[arg(short = 'f', value_parser = parse_flag, default_value_t = 0)]
        require: u16,
        /// Exclude records with ANY of these flag bits (decimal or 0x).
        #[arg(short = 'F', value_parser = parse_flag, default_value_t = 0)]
        exclude: u16,
        /// Skip records below this mapping quality.
        #[arg(short = 'q', default_value_t = 0)]
        min_mapq: u8,
        /// READ1 output file (de-interleave paired reads).
        #[arg(long = "r1")]
        out1: Option<PathBuf>,
        /// READ2 output file (de-interleave paired reads).
        #[arg(long = "r2")]
        out2: Option<PathBuf>,
        /// Singleton output (paired reads whose mate is unmapped).
        #[arg(short = 's', long)]
        singleton: Option<PathBuf>,
        /// Output for reads that are not part of a pair.
        #[arg(long = "r0")]
        out0: Option<PathBuf>,
        /// Do not append `/1` or `/2` to read names.
        #[arg(short = 'n', long = "no-suffix")]
        no_suffix: bool,
    },

    /// Extract reads to FASTA, flag-filtered (`samtools fasta`).
    #[command(before_help = BANNER)]
    Fasta {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Single/interleaved output (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Require ALL of these flag bits (decimal or 0x).
        #[arg(short = 'f', value_parser = parse_flag, default_value_t = 0)]
        require: u16,
        /// Exclude records with ANY of these flag bits (decimal or 0x).
        #[arg(short = 'F', value_parser = parse_flag, default_value_t = 0)]
        exclude: u16,
        /// Skip records below this mapping quality.
        #[arg(short = 'q', default_value_t = 0)]
        min_mapq: u8,
        /// READ1 output file.
        #[arg(long = "r1")]
        out1: Option<PathBuf>,
        /// READ2 output file.
        #[arg(long = "r2")]
        out2: Option<PathBuf>,
        /// Singleton output.
        #[arg(short = 's', long)]
        singleton: Option<PathBuf>,
        /// Output for unpaired reads.
        #[arg(long = "r0")]
        out0: Option<PathBuf>,
        /// Do not append `/1` or `/2` to read names.
        #[arg(short = 'n', long = "no-suffix")]
        no_suffix: bool,
    },

    /// Count records by SAM flag (`samtools flagstat`).
    #[command(before_help = BANNER)]
    Flagstat {
        /// Input SAM/BAM.
        input: PathBuf,
    },

    /// Per-reference mapped/unmapped read counts (`samtools idxstats`).
    #[command(before_help = BANNER)]
    Idxstats {
        /// Input SAM/BAM.
        input: PathBuf,
    },

    /// Merge and coordinate-sort SAM/BAM files (`samtools merge`).
    #[command(before_help = BANNER)]
    Merge {
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Output format: `bam` (default) or `sam`.
        #[arg(short = 'O', long, default_value = "bam")]
        format: String,
        /// Input SAM/BAM files (two or more).
        #[arg(required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,
    },

    /// Index a FASTA (`.fai`) or extract regions (`samtools faidx`).
    #[command(before_help = BANNER)]
    Faidx {
        /// Input FASTA.
        fasta: PathBuf,
        /// Output for extracted regions (default: stdout). Ignored when building
        /// the index (which is always written to `<fasta>.fai`).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Regions to extract (`chr` or `chr:beg-end`). If omitted, writes `.fai`.
        regions: Vec<String>,
    },

    /// Per-position read depth (`samtools depth`).
    #[command(before_help = BANNER)]
    Depth {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Skip reads below this mapping quality.
        #[arg(short = 'Q', long = "min-mapq", default_value_t = 0)]
        min_mapq: u8,
        /// Region filter: `chr` or `chr:beg-end`.
        #[arg(short = 'r', long)]
        region: Option<String>,
    },

    /// Text pileup of every covered position (`samtools mpileup`).
    #[command(before_help = BANNER)]
    Mpileup {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Reference FASTA (enables `.`/`,` matches and mismatch letters).
        #[arg(short = 'f', long = "fasta-ref")]
        reference: Option<PathBuf>,
        /// Skip reads below this mapping quality.
        #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
        min_mapq: u8,
        /// Skip bases below this base quality.
        #[arg(short = 'Q', long = "min-BQ", default_value_t = 13)]
        min_baseq: u8,
        /// Region filter: `chr` or `chr:beg-end`.
        #[arg(short = 'r', long)]
        region: Option<String>,
    },

    /// Alias of `mpileup` (`samtools pileup`).
    #[command(before_help = BANNER)]
    Pileup {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Reference FASTA.
        #[arg(short = 'f', long = "fasta-ref")]
        reference: Option<PathBuf>,
        /// Skip reads below this mapping quality.
        #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
        min_mapq: u8,
        /// Skip bases below this base quality.
        #[arg(short = 'Q', long = "min-BQ", default_value_t = 13)]
        min_baseq: u8,
        /// Region filter: `chr` or `chr:beg-end`.
        #[arg(short = 'r', long)]
        region: Option<String>,
    },

    /// Per-reference coverage summary table (`samtools coverage`).
    #[command(before_help = BANNER)]
    Coverage {
        /// Input SAM/BAM.
        input: PathBuf,
        /// Output path (default: stdout).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Region filter: `chr` or `chr:beg-end`.
        #[arg(short = 'r', long)]
        region: Option<String>,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum FmtArg {
    Fasta,
    Fastq,
    Sam,
    Bam,
}

/// htseq overlap-resolution mode (`cleaver count --mode`).
#[derive(Clone, Copy, ValueEnum)]
enum ModeArg {
    /// Any feature overlapping any base (htseq union / featureCounts default).
    Union,
    /// Meta-features whose features cover every base of the read.
    Strict,
    /// As strict, but bases covered by no feature are ignored.
    Nonempty,
}

/// VERSE quantification scheme for multiple feature types.
#[derive(Clone, Copy, ValueEnum, PartialEq)]
enum SchemeArg {
    /// Count each feature type separately (one matrix per type).
    Independent,
    /// Assign each read to the first feature type (in `-t` order) that matches.
    Hierarchical,
}

impl ModeArg {
    /// Map to the VERSE `-z` scheme (1 union, 2 strict, 3 nonempty).
    fn as_u8(self) -> u8 {
        match self {
            ModeArg::Union => 1,
            ModeArg::Strict => 2,
            ModeArg::Nonempty => 3,
        }
    }
}

impl From<FmtArg> for Format {
    fn from(f: FmtArg) -> Self {
        match f {
            FmtArg::Fasta => Format::Fasta,
            FmtArg::Fastq => Format::Fastq,
            FmtArg::Sam => Format::Sam,
            FmtArg::Bam => Format::Bam,
        }
    }
}

#[cfg_attr(feature = "hydra", allow(unused_variables))]
fn main() -> Result<()> {
    let cli = Cli::parse();
    #[cfg(feature = "hydra")]
    let cluster = cli.cluster;

    match cli.cmd {
        Cmd::Split { inputs, outdir, chunk_size, records, threads, format } => {
            let chunk_size = parse_size(&chunk_size)
                .with_context(|| format!("invalid --chunk-size '{chunk_size}'"))?;
            std::fs::create_dir_all(&outdir)
                .with_context(|| format!("creating output dir '{}'", outdir.display()))?;
            let forced = format.map(|f| jobs::fmt_name(Format::from(f)).to_string());
            let work: Vec<jobs::SplitJob> = inputs
                .iter()
                .map(|i| jobs::SplitJob {
                    input: i.clone(),
                    outdir: outdir.clone(),
                    chunk_size,
                    records,
                    format: forced.clone(),
                })
                .collect();

            #[cfg(feature = "hydra")]
            let results = {
                let cfg = hydra::build_config(
                    cluster.head,
                    cluster.client.as_deref(),
                    cluster.cpus,
                    cluster.sim_gpus,
                    cluster.port,
                    true,
                );
                hydra::run("cleaver_split", jobs::do_split, work, 1, 0, cfg)?
            };
            #[cfg(not(feature = "hydra"))]
            let results = parallel::run(&work, threads, |j| jobs::do_split(j.clone()))?;

            report_split(&results, &outdir)
        }

        Cmd::Convert { input, output, partition, mapped_only, unmapped_only } => {
            if partition {
                let m = with_tag(&output, "mapped");
                let u = with_tag(&output, "unmapped");
                let (nm, nu) = align::convert_alignment_split(&input, &m, &u)?;
                println!(
                    "partitioned {} -> {} ({nm} mapped) + {} ({nu} unmapped)",
                    input.display(),
                    m.display(),
                    u.display()
                );
            } else {
                let filter = if mapped_only {
                    align::MapFilter::Mapped
                } else if unmapped_only {
                    align::MapFilter::Unmapped
                } else {
                    align::MapFilter::All
                };
                let (n, kind) = align::convert(&input, &output, filter)?;
                println!("converted {} -> {} ({kind}, {n} records)", input.display(), output.display());
            }
            Ok(())
        }

        Cmd::Stats { inputs, format, threads } => {
            let forced = format.map(|f| jobs::fmt_name(Format::from(f)).to_string());
            let work: Vec<jobs::StatsJob> =
                inputs.iter().map(|i| jobs::StatsJob { input: i.clone(), format: forced.clone() }).collect();

            #[cfg(feature = "hydra")]
            let results = {
                let cfg = hydra::build_config(
                    cluster.head,
                    cluster.client.as_deref(),
                    cluster.cpus,
                    cluster.sim_gpus,
                    cluster.port,
                    true,
                );
                // Request 1 GPU/task so the scheduler pins a device for the kernel.
                hydra::run("cleaver_stats", jobs::do_stats, work, 1, 1, cfg)?
            };
            #[cfg(not(feature = "hydra"))]
            let results = parallel::run(&work, threads, |j| jobs::do_stats(j.clone()))?;

            report_stats(&results)
        }

        Cmd::Count {
            inputs,
            annotation,
            output,
            feature_type,
            scheme,
            group_by,
            stranded,
            min_mapq,
            count_multimappers,
            primary,
            allow_multi_overlap,
            mode,
            assign_mode,
            require_both_ends,
            exclude_chimeric,
            check_pe_dist,
            min_frag_len,
            max_frag_len,
            threads,
        } => {
            let z = assign_mode.unwrap_or(mode.as_u8());
            let types: Vec<String> = feature_type
                .split(';')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            if types.is_empty() {
                anyhow::bail!("no feature type given (-t)");
            }

            #[cfg(feature = "hydra")]
            let cfg = || {
                hydra::build_config(cluster.head, cluster.client.as_deref(), cluster.cpus, cluster.sim_gpus, cluster.port, true)
            };

            if types.len() > 1 && scheme == SchemeArg::Hierarchical {
                // one pass per file, trying feature types in priority order
                let work: Vec<jobs::HierJob> = inputs
                    .iter()
                    .map(|i| jobs::HierJob {
                        input: i.clone(),
                        annotation: annotation.clone(),
                        feature_types: types.clone(),
                        group_by: group_by.clone(),
                        stranded, min_mapq, count_multimappers, primary, allow_multi_overlap, mode: z,
                        require_both_ends, exclude_chimeric, check_pe_dist, min_frag_len, max_frag_len,
                    })
                    .collect();
                #[cfg(feature = "hydra")]
                let results = hydra::run("cleaver_count_hier", jobs::do_count_hier, work, 1, 0, cfg())?;
                #[cfg(not(feature = "hydra"))]
                let results = parallel::run(&work, threads, |j| jobs::do_count_hier(j.clone()))?;
                report_count_hier(&results, &output, &types)
            } else {
                // single type, or independent multi-type: one CountJob per (file, type)
                let mut work: Vec<jobs::CountJob> = Vec::new();
                for t in &types {
                    for i in &inputs {
                        work.push(jobs::CountJob {
                            input: i.clone(),
                            annotation: annotation.clone(),
                            feature_type: t.clone(),
                            group_by: group_by.clone(),
                            stranded, min_mapq, count_multimappers, primary, allow_multi_overlap, mode: z,
                            require_both_ends, exclude_chimeric, check_pe_dist, min_frag_len, max_frag_len,
                        });
                    }
                }
                #[cfg(feature = "hydra")]
                let results = hydra::run("cleaver_count", jobs::do_count, work, 1, 0, cfg())?;
                #[cfg(not(feature = "hydra"))]
                let results = parallel::run(&work, threads, |j| jobs::do_count(j.clone()))?;

                if types.len() == 1 {
                    report_count(&results, &output)
                } else {
                    // group by feature type, one matrix per type
                    let mut failures = 0;
                    for t in &types {
                        let sub: Vec<jobs::CountOutcome> =
                            results.iter().filter(|r| &r.feature_type == t).cloned().collect();
                        let out = typed_output(&output, t);
                        println!("[{t}]");
                        if let Err(e) = report_count(&sub, &out) {
                            eprintln!("  {t}: {e}");
                            failures += 1;
                        }
                    }
                    if failures > 0 {
                        anyhow::bail!("{failures} feature type(s) failed");
                    }
                    Ok(())
                }
            }
        }

        Cmd::Fastp {
            in1, out1, in2, out2, json,
            disable_adapter_trimming, adapter_sequence, adapter_sequence_r2, adapter_fasta, detect_adapter_for_pe,
            trim_front1, trim_tail1, max_len1, trim_front2, trim_tail2, max_len2,
            trim_poly_g, poly_g_min_len, trim_poly_x, poly_x_min_len,
            cut_front, cut_tail, cut_right, cut_window_size, cut_mean_quality,
            disable_quality_filtering, qualified_quality_phred, unqualified_percent_limit, n_base_limit, average_qual,
            disable_length_filtering, length_required, length_limit, low_complexity_filter, complexity_threshold,
            correction, overlap_len_require, overlap_diff_limit, overlap_diff_percent_limit,
            phred64, reads_to_process,
        } => {
            let adapter_fasta_seqs = match &adapter_fasta {
                Some(p) => read_fasta_seqs(p)?,
                None => Vec::new(),
            };
            let params = cleaver::fastp::FastpParams {
                disable_adapter: disable_adapter_trimming,
                adapter_r1: adapter_sequence.map(|s| s.into_bytes()),
                adapter_r2: adapter_sequence_r2.map(|s| s.into_bytes()),
                adapter_fasta: adapter_fasta_seqs,
                detect_adapter_for_pe,
                trim_front1, trim_tail1, max_len1, trim_front2, trim_tail2, max_len2,
                trim_poly_g, poly_g_min_len, trim_poly_x, poly_x_min_len,
                cut_front, cut_tail, cut_right, cut_window_size, cut_mean_quality,
                cut_front_window_size: None, cut_front_mean_quality: None,
                cut_tail_window_size: None, cut_tail_mean_quality: None,
                cut_right_window_size: None, cut_right_mean_quality: None,
                disable_quality_filtering, qualified_quality_phred, unqualified_percent_limit, n_base_limit, average_qual,
                disable_length_filtering, length_required, length_limit, low_complexity_filter, complexity_threshold,
                correction, overlap_len_require, overlap_diff_limit, overlap_diff_percent_limit,
                phred64, reads_to_process,
            };
            let report = match (&in2, &out2) {
                (Some(i2), Some(o2)) => cleaver::fastp::run_pe(&in1, &out1, i2, o2, &params)?,
                (None, None) => cleaver::fastp::run_se(&in1, &out1, &params)?,
                _ => anyhow::bail!("paired-end requires both --in2 and --out2"),
            };
            finish_fastp(report, &json)
        }

        Cmd::Demux { input, output, queries, window, min_score, min_score_diff, max_errors, no_trim, one_end } => {
            let barcodes = cleaver::demux::read_barcodes(&queries)?;
            let params = cleaver::demux::DemuxParams {
                window,
                min_score,
                min_score_diff,
                max_errors,
                trim: !no_trim,
                both_ends: !one_end,
            };
            let stats = cleaver::demux::demux_file(&input, &output, &barcodes, &params)?;
            let rate = if stats.total > 0 { stats.classified as f64 / stats.total as f64 * 100.0 } else { 0.0 };
            println!("  {} barcode(s) scanned over {} reads", barcodes.len(), stats.total);
            println!("  classified  : {:>10} ({:>5.1}%)", stats.classified, rate);
            println!("  unclassified: {:>10}", stats.unclassified);
            for (name, n) in &stats.per_barcode {
                println!("    {:<16} {:>10}", name, n);
            }
            println!("output -> {}/", output.display());
            Ok(())
        }

        Cmd::Samtools { cmd } => run_samtools(cmd),

        Cmd::Doctor => doctor::run(),

        Cmd::Info => {
            cmd_info();
            Ok(())
        }

        Cmd::Version => {
            print_version();
            Ok(())
        }
    }
}

fn fmt_is_bam(fmt: &str) -> Result<bool> {
    match fmt.to_ascii_lowercase().as_str() {
        "bam" | "b" => Ok(true),
        "sam" | "s" => Ok(false),
        other => anyhow::bail!("unknown output format '{other}' (use sam or bam)"),
    }
}

fn run_samtools(cmd: SamCmd) -> Result<()> {
    use cleaver::samtools as st;
    match cmd {
        SamCmd::View {
            input, output, bam, count, header_only, require, exclude, min_mapq, region, subsample,
        } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let (frac, seed) = match subsample {
                Some(s) => (Some(s.fract().max(0.0)), s.floor().max(0.0) as u64),
                None => (None, 0),
            };
            let p = st::ViewParams {
                filter: st::Filter { require, exclude, min_mapq },
                count,
                header_only,
                no_header: false,
                region,
                subsample: frac,
                seed,
                output,
                bam,
            };
            let n = st::view(&input, in_fmt, &p)?;
            if !count && !header_only {
                eprintln!("{n} record(s) written");
            }
            Ok(())
        }
        SamCmd::Sort { input, output, by_name, format, threads: _ } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let bam = fmt_is_bam(&format)?;
            let n = st::sort(&input, in_fmt, output.as_deref(), by_name, bam)?;
            eprintln!("sorted {n} record(s)");
            Ok(())
        }
        SamCmd::Index { input, output } => {
            let out = st::index(&input, output.as_deref())?;
            eprintln!("wrote {}", out.display());
            Ok(())
        }
        SamCmd::Fastq { input, output, require, exclude, min_mapq, out1, out2, singleton, out0, no_suffix } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let p = st::FastxParams {
                filter: st::Filter { require, exclude, min_mapq },
                fasta: false,
                out1,
                out2,
                singleton,
                out0,
                output,
                no_suffix,
            };
            let s = st::fastx(&input, in_fmt, &p)?;
            eprintln!(
                "{} read(s): R1 {}, R2 {}, singleton {}, single {}",
                s.total, s.read1, s.read2, s.singletons, s.single
            );
            Ok(())
        }
        SamCmd::Fasta { input, output, require, exclude, min_mapq, out1, out2, singleton, out0, no_suffix } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let p = st::FastxParams {
                filter: st::Filter { require, exclude, min_mapq },
                fasta: true,
                out1,
                out2,
                singleton,
                out0,
                output,
                no_suffix,
            };
            let s = st::fastx(&input, in_fmt, &p)?;
            eprintln!("{} read(s) written", s.total);
            Ok(())
        }
        SamCmd::Flagstat { input } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let fs = st::flagstat(&input, in_fmt)?;
            print!("{}", fs.report());
            Ok(())
        }
        SamCmd::Idxstats { input } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            for (name, len, mapped, unmapped) in st::idxstats(&input, in_fmt)? {
                println!("{name}\t{len}\t{mapped}\t{unmapped}");
            }
            Ok(())
        }
        SamCmd::Merge { output, format, inputs } => {
            let bam = fmt_is_bam(&format)?;
            let n = st::merge(&inputs, output.as_deref(), bam)?;
            eprintln!("merged {n} record(s) from {} file(s)", inputs.len());
            Ok(())
        }
        SamCmd::Faidx { fasta, output, regions } => {
            st::faidx(&fasta, &regions, output.as_deref())?;
            Ok(())
        }
        SamCmd::Depth { input, output, min_mapq, region } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let p = st::DepthParams { min_mapq, region, output };
            st::depth(&input, in_fmt, &p)?;
            Ok(())
        }
        SamCmd::Mpileup { input, output, reference, min_mapq, min_baseq, region } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let p = st::MpileupParams { min_mapq, min_baseq, region, reference, output };
            st::mpileup(&input, in_fmt, &p)?;
            Ok(())
        }
        SamCmd::Pileup { input, output, reference, min_mapq, min_baseq, region } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            let p = st::MpileupParams { min_mapq, min_baseq, region, reference, output };
            st::pileup(&input, in_fmt, &p)?;
            Ok(())
        }
        SamCmd::Coverage { input, output, region } => {
            let in_fmt = cleaver::formats::detect(&input)?;
            st::coverage(&input, in_fmt, &st::CoverageParams { region, output })?;
            Ok(())
        }
    }
}

fn report_split(results: &[jobs::SplitOutcome], outdir: &Path) -> Result<()> {
    let mut failures = 0;
    let mut total = 0;
    for r in results {
        if r.ok {
            total += r.chunks;
            println!("  {:<40} {:>6} {} chunks", r.input.display(), r.chunks, r.format);
        } else {
            failures += 1;
            eprintln!("  {:<40} ERROR: {}", r.input.display(), r.error);
        }
    }
    println!(
        "{} file(s) -> {} chunk(s) in {} ({} failed)",
        results.len(),
        total,
        outdir.display(),
        failures
    );
    if failures > 0 {
        anyhow::bail!("{failures} input(s) failed");
    }
    Ok(())
}

fn short(p: &Path) -> String {
    p.file_name().and_then(|s| s.to_str()).unwrap_or("?").to_string()
}

/// Read all sequences from a (possibly multi-record) FASTA file.
fn read_fasta_seqs(path: &Path) -> Result<Vec<Vec<u8>>> {
    let data = std::fs::read_to_string(path)
        .with_context(|| format!("reading adapter FASTA '{}'", path.display()))?;
    let mut seqs = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    for line in data.lines() {
        if line.starts_with('>') {
            if !cur.is_empty() {
                seqs.push(std::mem::take(&mut cur));
            }
        } else {
            cur.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    Ok(seqs)
}

/// Write the fastp JSON report and print a concise stdout summary.
fn finish_fastp(report: cleaver::fastp::Report, json: &Path) -> Result<()> {
    std::fs::write(json, report.to_json())
        .with_context(|| format!("writing JSON report '{}'", json.display()))?;
    let (b, a) = (&report.before, &report.after);
    let pct = |x: u64, t: u64| if t == 0 { 0.0 } else { x as f64 / t as f64 * 100.0 };
    let unit = if report.paired { "reads (mates)" } else { "reads" };
    println!("  {:<10} in : {:>12}  ({} bp)", unit, b.reads, b.bases);
    println!("  {:<10} out: {:>12}  ({} bp)", unit, a.reads, a.bases);
    println!("  Q20/Q30 in : {:>6.2}% / {:>6.2}%", pct(b.q20, b.bases), pct(b.q30, b.bases));
    println!("  Q20/Q30 out: {:>6.2}% / {:>6.2}%", pct(a.q20, a.bases), pct(a.q30, a.bases));
    println!(
        "  filtering : passed {}  low_qual {}  too_many_N {}  too_short {}  too_long {}  low_complexity {}",
        report.filter.passed, report.filter.low_quality, report.filter.too_many_n,
        report.filter.too_short, report.filter.too_long, report.filter.low_complexity
    );
    if report.adapter_trimmed_reads > 0 {
        println!(
            "  adapters  : trimmed from {} reads ({} bp)",
            report.adapter_trimmed_reads, report.adapter_trimmed_bases
        );
    }
    if report.corrected_bases > 0 {
        println!("  corrected : {} bases in PE overlaps", report.corrected_bases);
    }
    println!("  report    -> {}", json.display());
    Ok(())
}

fn report_stats(results: &[jobs::StatsOutcome]) -> Result<()> {
    println!(
        "  {:<22} {:>8} {:>13} {:>8} {:>8} {:>10} {:>8} {:>9} {:>5} {:>9} {:>6} {:>8}",
        "file", "seqs", "total_bp", "min", "max", "mean", "median", "N50", "L50", "N90", "GC%", "device"
    );
    let mut tot = Counts::default();
    let (mut tot_seqs, mut tot_bp) = (0u64, 0u64);
    let (mut ok, mut failures) = (0u32, 0u32);
    for r in results {
        if !r.ok {
            failures += 1;
            eprintln!("  {:<22} ERROR: {}", short(&r.input), r.error);
            continue;
        }
        ok += 1;
        tot.add(&Counts { a: r.a, c: r.c, g: r.g, t: r.t, n: r.n, other: r.other });
        tot_seqs += r.n_seqs;
        tot_bp += r.total_bp;
        println!(
            "  {:<22} {:>8} {:>13} {:>8} {:>8} {:>10.1} {:>8} {:>9} {:>5} {:>9} {:>6.2} {:>8}",
            short(&r.input),
            r.n_seqs,
            r.total_bp,
            r.min_len,
            r.max_len,
            r.mean_len,
            r.median_len,
            r.n50,
            r.l50,
            r.n90,
            r.gc_percent,
            r.device
        );
    }
    if ok > 1 {
        println!(
            "  {:<22} {:>8} {:>13} {:>8} {:>8} {:>10} {:>8} {:>9} {:>5} {:>9} {:>6.2} {:>8}",
            "TOTAL", tot_seqs, tot_bp, "", "", "", "", "", "", "", tot.gc() * 100.0, ""
        );
    }
    if failures > 0 {
        anyhow::bail!("{failures} input(s) failed");
    }
    Ok(())
}

/// Write the count matrix and `.summary`, and print a per-file assignment recap.
fn report_count(results: &[jobs::CountOutcome], output: &Path) -> Result<()> {
    let mut ok: Vec<&jobs::CountOutcome> = Vec::new();
    let mut failures = 0u32;
    for r in results {
        if r.ok {
            ok.push(r);
        } else {
            failures += 1;
            eprintln!("  {:<24} ERROR: {}", short(&r.input), r.error);
        }
    }
    if ok.is_empty() {
        anyhow::bail!("all {} input(s) failed", results.len());
    }
    let genes = &ok[0].gene_ids;
    for r in &ok {
        if r.gene_ids.len() != genes.len() {
            anyhow::bail!("inconsistent annotation across inputs (different meta-feature counts)");
        }
    }

    // Count matrix: gene_id + one column per input.
    let mut m = String::from("gene_id");
    for r in &ok {
        m.push('\t');
        m.push_str(&r.sample);
    }
    m.push('\n');
    for (i, g) in genes.iter().enumerate() {
        m.push_str(g);
        for r in &ok {
            m.push('\t');
            m.push_str(r.counts[i].to_string().as_str());
        }
        m.push('\n');
    }
    std::fs::write(output, m).with_context(|| format!("writing matrix '{}'", output.display()))?;

    // featureCounts-style summary.
    let summary_path = PathBuf::from(format!("{}.summary", output.display()));
    let mut s = String::from("Status");
    for r in &ok {
        s.push('\t');
        s.push_str(&r.sample);
    }
    s.push('\n');
    let rows: [(&str, fn(&jobs::CountOutcome) -> u64); 7] = [
        ("Assigned", |r| r.assigned),
        ("Unassigned_NoFeatures", |r| r.no_feature),
        ("Unassigned_Ambiguity", |r| r.ambiguous),
        ("Unassigned_MultiMapping", |r| r.multimapping),
        ("Unassigned_MappingQuality", |r| r.low_mapq),
        ("Unassigned_FragmentLength", |r| r.pe_filtered),
        ("Unassigned_Unmapped", |r| r.unmapped),
    ];
    for (label, f) in rows {
        s.push_str(label);
        for r in &ok {
            s.push('\t');
            s.push_str(f(r).to_string().as_str());
        }
        s.push('\n');
    }
    std::fs::write(&summary_path, s)
        .with_context(|| format!("writing summary '{}'", summary_path.display()))?;

    // Recap to stdout.
    println!(
        "  {} feature(s) -> {} meta-feature(s)",
        ok[0].n_features,
        genes.len()
    );
    for r in &ok {
        let rate = if r.total > 0 { r.assigned as f64 / r.total as f64 * 100.0 } else { 0.0 };
        println!(
            "  {:<24} {:>10} assigned ({:>5.1}%)   no_feat {:>9}  ambig {:>9}  multimap {:>9}  unmapped {:>9}",
            r.sample, r.assigned, rate, r.no_feature, r.ambiguous, r.multimapping, r.unmapped
        );
    }
    println!(
        "matrix  -> {} ({} x {})\nsummary -> {}",
        output.display(),
        genes.len(),
        ok.len(),
        summary_path.display()
    );
    if failures > 0 {
        anyhow::bail!("{failures} input(s) failed");
    }
    Ok(())
}

/// Insert `.ftype` before the extension of `base` (counts.tsv -> counts.exon.tsv).
fn typed_output(base: &Path, ftype: &str) -> PathBuf {
    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("counts");
    let ext = base.extension().and_then(|s| s.to_str()).unwrap_or("tsv");
    let parent = base.parent().filter(|p| !p.as_os_str().is_empty());
    let name = format!("{stem}.{ftype}.{ext}");
    match parent {
        Some(p) => p.join(name),
        None => PathBuf::from(name),
    }
}

/// Write per-feature-type matrices and a combined summary for a hierarchical run.
fn report_count_hier(results: &[jobs::HierOutcome], output: &Path, types: &[String]) -> Result<()> {
    let mut ok: Vec<&jobs::HierOutcome> = Vec::new();
    let mut failures = 0u32;
    for r in results {
        if r.ok {
            ok.push(r);
        } else {
            failures += 1;
            eprintln!("  {:<24} ERROR: {}", short(&r.input), r.error);
        }
    }
    if ok.is_empty() {
        anyhow::bail!("all {} input(s) failed", results.len());
    }

    for (ti, t) in types.iter().enumerate() {
        let genes = &ok[0].per_type[ti].gene_ids;
        let mut m = String::from("gene_id");
        for r in &ok {
            m.push('\t');
            m.push_str(&r.sample);
        }
        m.push('\n');
        for (gi, g) in genes.iter().enumerate() {
            m.push_str(g);
            for r in &ok {
                m.push('\t');
                m.push_str(r.per_type[ti].counts[gi].to_string().as_str());
            }
            m.push('\n');
        }
        let out = typed_output(output, t);
        std::fs::write(&out, m).with_context(|| format!("writing matrix '{}'", out.display()))?;
        println!("  [{t}] {} meta-feature(s) -> {}", genes.len(), out.display());
    }

    let summary_path = PathBuf::from(format!("{}.summary", output.display()));
    let mut s = String::from("Status");
    for r in &ok {
        s.push('\t');
        s.push_str(&r.sample);
    }
    s.push('\n');
    for (ti, t) in types.iter().enumerate() {
        s.push_str(&format!("Assigned_{t}"));
        for r in &ok {
            s.push('\t');
            s.push_str(r.per_type[ti].assigned.to_string().as_str());
        }
        s.push('\n');
    }
    let grows: [(&str, fn(&jobs::HierOutcome) -> u64); 6] = [
        ("Unassigned_NoFeatures", |r| r.no_feature),
        ("Unassigned_Ambiguity", |r| r.ambiguous),
        ("Unassigned_MultiMapping", |r| r.multimapping),
        ("Unassigned_MappingQuality", |r| r.low_mapq),
        ("Unassigned_FragmentLength", |r| r.pe_filtered),
        ("Unassigned_Unmapped", |r| r.unmapped),
    ];
    for (label, f) in grows {
        s.push_str(label);
        for r in &ok {
            s.push('\t');
            s.push_str(f(r).to_string().as_str());
        }
        s.push('\n');
    }
    std::fs::write(&summary_path, s)
        .with_context(|| format!("writing summary '{}'", summary_path.display()))?;

    for r in &ok {
        let assigned: u64 = r.per_type.iter().map(|t| t.assigned).sum();
        let rate = if r.total > 0 { assigned as f64 / r.total as f64 * 100.0 } else { 0.0 };
        print!("  {:<24} {:>10} assigned ({:>5.1}%)  [", r.sample, assigned, rate);
        for (ti, t) in types.iter().enumerate() {
            print!("{}{}={}", if ti > 0 { " " } else { "" }, t, r.per_type[ti].assigned);
        }
        println!("]  no_feat {}  ambig {}", r.no_feature, r.ambiguous);
    }
    println!("summary -> {}", summary_path.display());
    if failures > 0 {
        anyhow::bail!("{failures} input(s) failed");
    }
    Ok(())
}

fn print_version() {
    print!("{BANNER}");
    println!("cleaver {}", env!("CARGO_PKG_VERSION"));
    println!(
        "backend  {}",
        if cfg!(feature = "hydra") { "HydraMPP (multi-core / cross-node)" } else { "rayon (in-node)" }
    );
    println!(
        "gpu      {}",
        if gpu::compute_enabled() { "compiled (stats kernel)" } else { "not compiled (build --features gpu)" }
    );
}

/// Insert `.tag` before the final extension: `aln.bam` -> `aln.mapped.bam`.
fn with_tag(p: &Path, tag: &str) -> PathBuf {
    let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("out");
    let name = match p.extension().and_then(|s| s.to_str()) {
        Some(e) => format!("{stem}.{tag}.{e}"),
        None => format!("{stem}.{tag}"),
    };
    match p.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.join(name),
        _ => PathBuf::from(name),
    }
}

fn cmd_info() {
    print!("{BANNER}");
    println!("version          {}", env!("CARGO_PKG_VERSION"));
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    println!("logical cores    {cores}");
    println!(
        "scaling backend  {}",
        if cfg!(feature = "hydra") { "HydraMPP (multi-core / cross-node)" } else { "rayon (in-node)" }
    );
    println!(
        "gpu kernel       {}",
        if gpu::compute_enabled() { "compiled (stats)" } else { "not compiled (--features gpu)" }
    );
    let gpus = gpu::detect();
    if gpus.is_empty() {
        println!("gpus             none detected");
    } else {
        for (i, g) in gpus.iter().enumerate() {
            println!("gpu[{i}]           {} ({} MiB)", g.name, g.mem_mib);
        }
    }
    println!("formats          FASTA (.fasta/.fa/.fna/.ffn/.faa/.frn/.mpfa), FASTQ (.fastq/.fq), SAM, BAM (+.gz)");
    println!("convert          SAM<->BAM (--mapped-only/--unmapped-only/--partition), FASTQ->FASTA");
    println!("stats            genome stats: N50/L50/N90, lengths, GC% (GPU base-comp with --features gpu)");
    println!("count            featureCounts/htseq-style reads-per-feature from BAM/SAM + GTF/GFF");
    println!("fastp            FASTQ preprocessing: adapter/quality/polyX trim + filtering (pure-Rust)");
    println!("demux            ONT barcode demultiplexing: fit-align + trim + split (pure-Rust barbell)");
}
