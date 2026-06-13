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

use cleaver_core::{align, compute::Counts, parse_size, Format};

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
        /// Feature type in column 3 to count (e.g. exon, CDS).
        #[arg(short = 't', long, default_value = "exon")]
        feature_type: String,
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
        /// Count reads overlapping >1 meta-feature against all of them.
        #[arg(short = 'O', long)]
        allow_multi_overlap: bool,
        /// Overlap-resolution mode.
        #[arg(long, value_enum, default_value_t = ModeArg::Union)]
        mode: ModeArg,
        /// Worker threads (one file per task) for the in-node backend. 0 = all cores.
        #[arg(short = 'T', long, default_value_t = 0)]
        threads: usize,
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

impl ModeArg {
    fn as_u8(self) -> u8 {
        match self {
            ModeArg::Union => 0,
            ModeArg::Strict => 1,
            ModeArg::Nonempty => 2,
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
            group_by,
            stranded,
            min_mapq,
            count_multimappers,
            allow_multi_overlap,
            mode,
            threads,
        } => {
            let work: Vec<jobs::CountJob> = inputs
                .iter()
                .map(|i| jobs::CountJob {
                    input: i.clone(),
                    annotation: annotation.clone(),
                    feature_type: feature_type.clone(),
                    group_by: group_by.clone(),
                    stranded,
                    min_mapq,
                    count_multimappers,
                    allow_multi_overlap,
                    mode: mode.as_u8(),
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
                hydra::run("cleaver_count", jobs::do_count, work, 1, 0, cfg)?
            };
            #[cfg(not(feature = "hydra"))]
            let results = parallel::run(&work, threads, |j| jobs::do_count(j.clone()))?;

            report_count(&results, &output)
        }

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
    let rows: [(&str, fn(&jobs::CountOutcome) -> u64); 6] = [
        ("Assigned", |r| r.assigned),
        ("Unassigned_NoFeatures", |r| r.no_feature),
        ("Unassigned_Ambiguity", |r| r.ambiguous),
        ("Unassigned_MultiMapping", |r| r.multimapping),
        ("Unassigned_MappingQuality", |r| r.low_mapq),
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
}
