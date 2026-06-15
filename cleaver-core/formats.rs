//! Format detection for the sequence/alignment files Cleaver handles.
//!
//! Resolution order: file extension first (after stripping a `.gz` suffix),
//! then a content sniff of the first bytes if the extension is unknown.

use std::fs::File;
use std::io::Read;
use std::path::Path;

use anyhow::{Context, Result};

/// A file format Cleaver can split and/or convert.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// FASTA nucleotide/protein (.fasta/.fa/.fna/.ffn/.faa/.frn/.mpfa).
    Fasta,
    /// FASTQ reads (.fastq/.fq).
    Fastq,
    /// SAM alignments (text).
    Sam,
    /// BAM alignments (BGZF-compressed binary).
    Bam,
}

impl Format {
    /// The chunk file extension this format writes (without the dot).
    pub fn ext(self) -> &'static str {
        match self {
            Format::Fasta => "fasta",
            Format::Fastq => "fastq",
            Format::Sam => "sam",
            Format::Bam => "bam",
        }
    }
    /// True for the line-oriented sequence formats (split by size).
    pub fn is_sequence(self) -> bool {
        matches!(self, Format::Fasta | Format::Fastq)
    }
    /// True for the alignment formats (split by record count).
    pub fn is_alignment(self) -> bool {
        matches!(self, Format::Sam | Format::Bam)
    }
}

/// Lowercased final extension of `name`, ignoring a trailing `.gz`.
fn inner_ext(name: &str) -> Option<String> {
    let base = name.strip_suffix(".gz").or_else(|| name.strip_suffix(".bgz")).unwrap_or(name);
    base.rsplit_once('.').map(|(_, e)| e.to_ascii_lowercase())
}

/// Detect by extension only (works for paths that need not yet exist, e.g. an
/// output target). Returns `None` if the extension is not recognised.
pub fn detect_by_ext(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?;
    match inner_ext(name)?.as_str() {
        "fasta" | "fa" | "fna" | "ffn" | "faa" | "frn" | "mpfa" | "fas" => Some(Format::Fasta),
        "fastq" | "fq" => Some(Format::Fastq),
        "sam" => Some(Format::Sam),
        "bam" => Some(Format::Bam),
        _ => None,
    }
}

/// Detect the format of an existing input: extension first, then content.
pub fn detect(path: &Path) -> Result<Format> {
    if let Some(f) = detect_by_ext(path) {
        return Ok(f);
    }
    sniff(path).with_context(|| format!("could not determine format of '{}'", path.display()))
}

/// Content-based detection from the first bytes (gzip-aware).
fn sniff(path: &Path) -> Result<Format> {
    let mut file = File::open(path).with_context(|| format!("opening '{}'", path.display()))?;
    let mut magic = [0u8; 4];
    let n = file.read(&mut magic)?;
    let head = &magic[..n];

    // gzip / BGZF member starts with 1f 8b. BGZF additionally sets FLG.FEXTRA
    // (byte 3 == 0x04) and carries a "BAM\1" payload — treat as BAM.
    if head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        // Re-open and decompress the first block to read the inner magic.
        let f2 = File::open(path)?;
        #[cfg(feature = "gzip")]
        {
            let mut dec = flate2::read::MultiGzDecoder::new(f2);
            let mut inner = [0u8; 16];
            let m = dec.read(&mut inner).unwrap_or(0);
            return Ok(classify_text(&inner[..m]));
        }
        #[cfg(not(feature = "gzip"))]
        {
            let _ = f2;
            // BGZF FEXTRA flag strongly implies BAM even without decompression.
            if head.len() >= 4 && head[3] == 0x04 {
                return Ok(Format::Bam);
            }
            anyhow::bail!("gzip input but built without the `gzip` feature");
        }
    }
    Ok(classify_text(head))
}

/// Classify decompressed/plain leading bytes.
fn classify_text(bytes: &[u8]) -> Format {
    // skip leading whitespace
    let start = bytes.iter().position(|b| !b.is_ascii_whitespace()).unwrap_or(0);
    let b = &bytes[start..];
    if b.starts_with(b"BAM\x01") {
        return Format::Bam;
    }
    if b.first() == Some(&b'>') {
        return Format::Fasta;
    }
    if b.starts_with(b"@HD") || b.starts_with(b"@SQ") || b.starts_with(b"@RG") || b.starts_with(b"@PG") || b.starts_with(b"@CO") {
        return Format::Sam;
    }
    if b.first() == Some(&b'@') {
        return Format::Fastq;
    }
    // default to FASTA (most permissive line-oriented sequence format)
    Format::Fasta
}
