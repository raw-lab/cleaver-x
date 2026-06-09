//! FASTA discovery and loading helpers.

use anyhow::{Context, Result};
use bio::io::fasta;
use std::path::{Path, PathBuf};

/// Accepted FASTA extensions (case-insensitive).
const FASTA_EXTS: &[&str] = &["fasta", "fa", "fna", "ffn", "faa"];

/// Return every FASTA file under `input`. If `input` is itself a file, return it.
pub fn collect_fasta_files(input: &Path) -> Result<Vec<PathBuf>> {
    if input.is_file() {
        return Ok(vec![input.to_path_buf()]);
    }
    if !input.is_dir() {
        anyhow::bail!("input {:?} is neither a file nor a directory", input);
    }

    let mut out = Vec::new();
    for entry in std::fs::read_dir(input)
        .with_context(|| format!("reading dir {:?}", input))?
    {
        let p = entry?.path();
        if !p.is_file() {
            continue;
        }
        if let Some(ext) = p.extension().and_then(|e| e.to_str()) {
            let ext_lower = ext.to_ascii_lowercase();
            if FASTA_EXTS.contains(&ext_lower.as_str()) {
                out.push(p);
            }
        }
    }
    out.sort();
    if out.is_empty() {
        anyhow::bail!("no FASTA files found in {:?}", input);
    }
    Ok(out)
}

/// Load every record from every file into memory.
///
/// For very large genomes, consider switching this to a streaming pass —
/// the codon + stats code already operates per-record so it is friendly to that.
pub fn load_all_records(files: &[PathBuf]) -> Result<Vec<fasta::Record>> {
    let mut all = Vec::new();
    for f in files {
        let reader = fasta::Reader::from_file(f)
            .with_context(|| format!("opening FASTA {:?}", f))?;
        for r in reader.records() {
            all.push(r.with_context(|| format!("parsing record in {:?}", f))?);
        }
    }
    Ok(all)
}
