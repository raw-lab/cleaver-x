//! FragGeneScanRs subprocess wrapper.
//!
//! FragGeneScanRs is the Rust port of FragGeneScan by unipept
//! (<https://github.com/unipept/FragGeneScanRs>). Install via either:
//!
//! ```text
//! cargo install fraggenescanrs
//! # or
//! conda install -c bioconda fraggenescanrs
//! ```
//!
//! We shell out rather than linking it as a library, so users can swap in a
//! specific FGS version or use the classic C FragGeneScan binary.
//! Expected invocation:
//!
//! ```text
//! FragGeneScanRs -s <input.fa> -o <out_prefix> -w -t <model> -p <threads>
//! ```
//!
//! The `-w` flag selects whole-genome mode (the default here). For short reads
//! pass `--model illumina_10` etc. — that model name is forwarded as `-t` and
//! `-w` is omitted.

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Runs FragGeneScanRs on `input_fasta` and returns the path to the predicted
/// nucleotide ORF file (`<prefix>.ffn`).
pub fn run_fgs(
    input_fasta: &Path,
    outdir: &Path,
    fgs_bin: &str,
    model: &str,
    threads: usize,
) -> Result<PathBuf> {
    let out_prefix = outdir.join("fgs_predicted");
    let whole_genome = model.eq_ignore_ascii_case("complete");

    let mut cmd = Command::new(fgs_bin);
    cmd.arg("-s").arg(input_fasta)
        .arg("-o").arg(&out_prefix)
        .arg("-t").arg(model)
        .arg("-p").arg(threads.to_string());
    if whole_genome {
        cmd.arg("-w");
    }

    eprintln!("[rustyomestats] running: {:?}", cmd);
    let status = cmd
        .status()
        .with_context(|| format!(
            "failed to launch FragGeneScanRs binary `{}`. \
             Install with `cargo install fraggenescanrs` or \
             `conda install -c bioconda fraggenescanrs`, or pass --fgs-bin <path>.",
            fgs_bin
        ))?;

    if !status.success() {
        anyhow::bail!("FragGeneScanRs exited non-zero: {:?}", status);
    }

    let ffn = out_prefix.with_extension("ffn");
    if !ffn.exists() {
        anyhow::bail!(
            "expected FragGeneScanRs output {:?} not found. \
             Some FGS versions need `-n` to emit nucleotide ORFs — \
             check the installed FGS CLI and adjust if needed.",
            ffn
        );
    }
    Ok(ffn)
}

/// Run FragGeneScanRs, then count codon density on its predicted ORFs.
pub fn compute_predicted_codon_density(
    input_fasta: &Path,
    outdir: &Path,
    fgs_bin: &str,
    model: &str,
    threads: usize,
) -> Result<()> {
    let ffn = run_fgs(input_fasta, outdir, fgs_bin, model, threads)?;
    crate::codon::compute_predicted_codon_density(&ffn, outdir)
}
