//! Embedding cleaver as a library: genome stats, flag tallies, size parsing.
//!
//! Run with:  cargo run --example library_usage -- genome.fna aln.bam

use anyhow::Result;
use cleaver::{count, genome, parse_size, samtools, Format};
use std::path::Path;

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let fasta = args.next().unwrap_or_else(|| "genome.fna".into());
    let bam = args.next();

    // ---- genome / assembly statistics (N50, L50, GC%, lengths) --------------
    let stats = genome::genome_stats_file(Path::new(&fasta), Format::Fasta, None)?;
    println!(
        "{}: {} sequences, {} bp, N50 {} (L50 {}), GC {:.2}%",
        fasta,
        stats.n_seqs,
        stats.total_bp,
        stats.n50(),
        stats.l50(),
        stats.gc_percent
    );

    // ---- chunk-size strings ("1G", "256M", "50Mi") -------------------------
    println!("256M = {} bytes", parse_size("256M")?);

    // ---- alignment summaries, straight from the samtools engine ------------
    if let Some(bam) = bam {
        let fs = samtools::flagstat(Path::new(&bam), Format::Bam)?;
        print!("{}", fs.report());

        // reads-per-feature needs an annotation; shown here with defaults
        let _params = count::CountParams::default();
    }

    Ok(())
}
