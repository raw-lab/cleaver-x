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
use clap::{ArgAction, CommandFactory, Parser, Subcommand, ValueEnum};

use cleaver_core::{align, compute::Counts, parse_size, Format};

pub(crate) const BANNER: &str = r"
   ___ _
  / __\ | ___  __ ___   _____ _ __
 / /  | |/ _ \/ _` \ \ / / _ \ '__|
/ /___| |  __/ (_| |\ V /  __/ |
\____/|_|\___|\__,_| \_/ \___|_|   record-aware FASTA/FASTQ/SAM/BAM
";

/// Shown under top-level `-h` / `help`.
const EXAMPLES: &str = "\
EXAMPLES:
  # Split into record-aware chunks (never cuts a record)
  cleaver split genome.fasta -o chunks/ -c 500M      # ~500 MB FASTA chunks
  cleaver split reads.fastq.gz -o chunks/ -c 200M    # gzip in, FASTQ by size
  cleaver split aln.bam -o chunks/ -r 1000000        # SAM/BAM by record count

  # Genome statistics: N50/L50/N90, min/max/mean/median length, GC%, A/C/G/T/N
  cleaver stats genome.fasta
  cleaver stats *.fasta *.fastq.gz

  # Convert: SAM<->BAM, mapped/unmapped partition, FASTQ->FASTA
  cleaver convert aln.sam aln.bam
  cleaver convert aln.bam out.sam --partition
  cleaver convert reads.fastq reads.fasta

  cleaver info        # backends, GPUs, formats        cleaver version   (-v)
  cleaver doctor      # self-test the install          cleaver help <cmd>
";

const SPLIT_HELP: &str = "\
EXAMPLES:
  cleaver split genome.fasta -o chunks/ -c 500M       # ~500 MB FASTA chunks
  cleaver split reads.fastq.gz -o chunks/ -c 200M     # gzip in, FASTQ by size
  cleaver split aln.bam -o chunks/ -r 500000          # 500k records / chunk
  cleaver split a.fasta b.fasta -o chunks/            # many inputs, one pool
";

const STATS_HELP: &str = "\
EXAMPLES:
  cleaver stats genome.fasta              # N50/L50/N90, min/max/mean, GC%, A/C/G/T/N
  cleaver stats *.fasta                   # one row per file + a TOTAL row
  cleaver stats reads.fastq.gz --format fastq

The N/L assembly metrics are computed by rustyomestats' compute_nl kernel. Build
cleaver --features full (rustc 1.85+) to add rustyomestats' codon density, Castro
U50 metrics, and FragGeneScan ORF density.
";

const CONVERT_HELP: &str = "\
EXAMPLES:
  cleaver convert aln.sam aln.bam              # SAM -> BAM (BGZF-compressed)
  cleaver convert aln.bam aln.sam              # BAM -> SAM
  cleaver convert aln.bam kept.bam --mapped-only
  cleaver convert aln.sam out.sam --partition  # -> out.mapped.sam + out.unmapped.sam
  cleaver convert reads.fastq reads.fasta      # drop quality lines
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
    before_help = BANNER,
    about = "Split & convert FASTA/FASTQ/SAM/BAM without cutting records.",
    after_help = EXAMPLES,
    disable_help_flag = true,
    disable_version_flag = true
)]
struct Cli {
    /// Print help (use `cleaver help <command>` for one command).
    #[arg(short = 'h', global = true, action = ArgAction::Help)]
    help: Option<bool>,
    /// Print version and build configuration.
    #[arg(short = 'v', action = ArgAction::SetTrue)]
    version: bool,
    #[cfg(feature = "hydra")]
    #[command(flatten)]
    cluster: Cluster,
    #[command(subcommand)]
    cmd: Option<Cmd>,
}

#[derive(Subcommand)]
enum Cmd {
    /// Split one or more files into record-aware chunks (format auto-detected).
    #[command(before_help = BANNER, after_help = SPLIT_HELP)]
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
    #[command(before_help = BANNER, after_help = CONVERT_HELP)]
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
    /// Genome statistics: N50/L50/N90, lengths, GC%, composition (via rustyomestats).
    #[command(before_help = BANNER, after_help = STATS_HELP)]
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
    /// Print version and build configuration.
    Version,
    /// Self-test: verify the install and that every code path works.
    Doctor,
    /// Show version, backends, detected GPUs, and supported formats.
    Info,
}

#[derive(Clone, Copy, ValueEnum)]
enum FmtArg {
    Fasta,
    Fastq,
    Sam,
    Bam,
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

    if cli.version {
        print_version();
        return Ok(());
    }

    #[cfg(feature = "hydra")]
    let cluster = cli.cluster;

    let cmd = match cli.cmd {
        Some(c) => c,
        None => {
            // bare `cleaver`: show the full help (with examples).
            Cli::command().print_help()?;
            println!();
            return Ok(());
        }
    };

    match cmd {
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

        Cmd::Version => {
            print_version();
            Ok(())
        }

        Cmd::Doctor => doctor::run(),

        Cmd::Info => {
            cmd_info();
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

fn report_stats(results: &[jobs::StatsOutcome]) -> Result<()> {
    println!(
        "  {:<22} {:>8} {:>13} {:>8} {:>9} {:>10} {:>8} {:>10} {:>5} {:>9} {:>6} {:>11} {:>11} {:>11} {:>11} {:>9} {:>10}",
        "file", "seqs", "total_bp", "min", "max", "mean", "median", "N50", "L50", "N90", "GC%", "A", "C", "G", "T", "N", "device"
    );
    let mut tot = Counts::default();
    let mut tot_seqs = 0u64;
    let mut tot_bp = 0u64;
    let mut gmin = u64::MAX;
    let mut gmax = 0u64;
    let mut ok = 0;
    let mut failures = 0;

    for r in results {
        let name = r.input.file_name().and_then(|s| s.to_str()).unwrap_or("?");
        if !r.ok {
            failures += 1;
            eprintln!("  {name:<22} ERROR: {}", r.error);
            continue;
        }
        ok += 1;
        tot.add(&Counts { a: r.a, c: r.c, g: r.g, t: r.t, n: r.n, other: r.other });
        tot_seqs += r.n_seqs;
        tot_bp += r.total_bp;
        if r.n_seqs > 0 {
            gmin = gmin.min(r.min_len);
            gmax = gmax.max(r.max_len);
        }
        println!(
            "  {:<22} {:>8} {:>13} {:>8} {:>9} {:>10.1} {:>8} {:>10} {:>5} {:>9} {:>6.2} {:>11} {:>11} {:>11} {:>11} {:>9} {:>10}",
            name, r.n_seqs, r.total_bp, r.min_len, r.max_len, r.mean_len, r.median_len,
            r.n50, r.l50, r.n90, r.gc_percent, r.a, r.c, r.g, r.t, r.n, r.device
        );
    }

    if ok > 1 {
        let mean = if tot_seqs > 0 { tot_bp as f64 / tot_seqs as f64 } else { 0.0 };
        println!(
            "  {:<22} {:>8} {:>13} {:>8} {:>9} {:>10.1} {:>8} {:>10} {:>5} {:>9} {:>6.2} {:>11} {:>11} {:>11} {:>11} {:>9} {:>10}",
            "TOTAL", tot_seqs, tot_bp, if gmin == u64::MAX { 0 } else { gmin }, gmax, mean,
            "-", "-", "-", "-", 100.0 * tot.gc(), tot.a, tot.c, tot.g, tot.t, tot.n, ""
        );
    }
    if failures > 0 {
        anyhow::bail!("{failures} input(s) failed");
    }
    Ok(())
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
    println!("stats            genome stats: N50/L50/N90, lengths, GC%, A/C/G/T/N (rustyomestats compute_nl)");
}

/// Print the ASCII banner, version, and build configuration (`version` / `-v`).
fn print_version() {
    print!("{BANNER}");
    println!("cleaver {}", env!("CARGO_PKG_VERSION"));
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    println!(
        "  backend {} | gpu kernel {} | {} logical core(s)",
        if cfg!(feature = "hydra") {
            "HydraMPP (multi-core/cross-node)"
        } else {
            "rayon (in-node)"
        },
        if gpu::compute_enabled() { "compiled" } else { "not compiled" },
        cores
    );
}
