//! rustyomestats CLI.
//!
//! Two subcommands:
//!   * `genome` — length / GC / N-L / 6-frame codon density / FragGeneScan
//!   * `u50`    — Castro et al. (2016) N50/NG50/L50/LG50/U50/UL50/UG50/
//!                ULG50/UG50% from a reference FASTA + sorted mapped-contigs BED

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use rustyomestats::{codon, fgs, io_utils, stats, u50};

#[derive(Parser, Debug)]
#[command(
    name    = "rustyomestats",
    version,
    about   = "Genome statistics, codon density, and Castro U50 assembly metrics.",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Length, GC, N/L stats, 6-frame + FragGeneScan codon density from FASTA.
    Genome(GenomeArgs),
    /// Castro et al. (2016) U50-family assembly metrics from reference + BED.
    U50(U50Args),
}

// --------------------------------------------------------------------------- //
// `genome` subcommand                                                          //
// --------------------------------------------------------------------------- //

#[derive(Parser, Debug)]
struct GenomeArgs {
    /// FASTA file or folder of FASTA files
    #[arg(short, long)]
    fasta: PathBuf,

    /// Interval size in bp for the length histogram
    #[arg(short, long, default_value_t = 1000)]
    interval: usize,

    /// Output directory
    #[arg(short, long, default_value = "rustyomestats_out")]
    outdir: PathBuf,

    /// Skip 6-frame absolute codon density
    #[arg(long)]
    skip_absolute: bool,

    /// Skip FragGeneScan-predicted codon density
    #[arg(long)]
    skip_predicted: bool,

    /// Path or name of the FragGeneScanRs binary
    #[arg(long, default_value = "FragGeneScanRs")]
    fgs_bin: String,

    /// FragGeneScan training model. `complete` = assembled genome;
    /// `illumina_5`, `illumina_10`, `454_30`, … for short reads.
    #[arg(long, default_value = "complete")]
    fgs_model: String,

    /// Threads (0 = all available)
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
}

fn run_genome(a: GenomeArgs) -> Result<()> {
    if a.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(a.threads)
            .build_global()
            .ok();
    }

    std::fs::create_dir_all(&a.outdir)
        .with_context(|| format!("creating output dir {:?}", a.outdir))?;

    let files = io_utils::collect_fasta_files(&a.fasta)?;
    eprintln!("[genome] found {} FASTA file(s)", files.len());
    let records = io_utils::load_all_records(&files)?;
    eprintln!("[genome] loaded {} sequences", records.len());

    let summary = stats::compute_length_stats(&records, a.interval, &a.outdir)?;
    println!("\n=== genome summary ===\n{summary}\n");

    if !a.skip_absolute {
        codon::compute_absolute_codon_density(&records, &a.outdir)?;
        eprintln!("[genome] 6-frame absolute codon density written");
    }

    if !a.skip_predicted {
        let fgs_input: PathBuf = if a.fasta.is_file() {
            a.fasta.clone()
        } else {
            let concat = a.outdir.join("fgs_input_concat.fa");
            use std::io::{BufWriter, Write};
            let mut out = BufWriter::new(std::fs::File::create(&concat)?);
            for f in &files {
                let src = std::fs::read(f)
                    .with_context(|| format!("reading {:?}", f))?;
                out.write_all(&src)?;
                if !src.ends_with(b"\n") {
                    out.write_all(b"\n")?;
                }
            }
            out.flush()?;
            concat
        };
        match fgs::compute_predicted_codon_density(
            &fgs_input, &a.outdir, &a.fgs_bin, &a.fgs_model, a.threads.max(1),
        ) {
            Ok(()) => eprintln!("[genome] FGS predicted codon density written"),
            Err(e) => eprintln!("[genome] WARN: FGS step skipped: {e:#}"),
        }
    }

    if !a.skip_absolute && !a.skip_predicted {
        if let Err(e) = codon::write_codon_comparison(&a.outdir) {
            eprintln!("[genome] WARN: could not write comparison: {e:#}");
        }
    }

    eprintln!("[genome] done. outputs in {:?}", a.outdir);
    eprintln!("[genome] plot:  python scripts/plot_stats.py -d {:?}", a.outdir);
    Ok(())
}

// --------------------------------------------------------------------------- //
// `u50` subcommand                                                             //
// --------------------------------------------------------------------------- //

#[derive(Parser, Debug)]
struct U50Args {
    /// Reference FASTA (single record — first record is used as the baseline).
    #[arg(short, long)]
    reference: PathBuf,

    /// Sorted BED file of contigs mapped to the reference.
    /// 0-based half-open (standard UCSC BED); column 1 is chrom, 2 is start,
    /// 3 is stop. Any extra columns are ignored.
    #[arg(short, long)]
    bed: PathBuf,

    /// Output directory
    #[arg(short, long, default_value = "rustyomestats_out")]
    outdir: PathBuf,
}

fn run_u50(a: U50Args) -> Result<()> {
    let res = u50::compute_u50(&a.reference, &a.bed, &a.outdir)?;

    println!("\n=== Castro U50 assembly summary ===");
    println!("reference                 : {} ({} bp)", res.ref_name, res.ref_length);
    println!("contigs in BED            : {}", res.num_contigs_in_bed);
    println!("contigs after dedup       : {}", res.num_contigs_kept);
    println!("contigs with unique bp    : {}", res.num_contigs_unique);
    println!("total original length     : {}", res.total_orig_length);
    println!("total unique length       : {}", res.total_unique_length);
    println!("gap bp                    : {}", res.gap_bp);
    println!("overlap bp                : {}", res.overlap_bp);
    println!();
    println!("  N50 = {:>10}   L50  = {:>6}", res.n50,  res.l50);
    println!("  NG50= {:>10}   LG50 = {:>6}", res.ng50, res.lg50);
    println!("  U50 = {:>10}   UL50 = {:>6}", res.u50,  res.ul50);
    println!("  UG50= {:>10}   ULG50= {:>6}", res.ug50, res.ulg50);
    println!("  UG50% = {:.4} %", res.ug50_pct);

    eprintln!("\n[u50] outputs written to {:?}", a.outdir);
    eprintln!("[u50] plot:  python scripts/plot_stats.py -d {:?}", a.outdir);
    Ok(())
}

// --------------------------------------------------------------------------- //
// entry point                                                                  //
// --------------------------------------------------------------------------- //

fn main() -> Result<()> {
    match Cli::parse().cmd {
        Cmd::Genome(a) => run_genome(a),
        Cmd::U50(a)    => run_u50(a),
    }
}
