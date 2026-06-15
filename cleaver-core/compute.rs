//! Base-composition counting (A/C/G/T/N/other + GC%).
//!
//! This is the CPU reference implementation: pure, `forbid(unsafe)`-clean, and
//! streamed so memory stays flat. The CLI's `--features gpu` path runs the same
//! arithmetic as a CUDA kernel over the sequence bytes and falls back here when
//! no device is available, so results are identical either way.

use std::path::Path;

use anyhow::Result;

use crate::formats::Format;

/// Nucleotide tallies over sequence bytes (case-insensitive).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Counts {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub other: u64,
}

impl Counts {
    pub fn total(&self) -> u64 {
        self.a + self.c + self.g + self.t + self.n + self.other
    }
    /// GC fraction over A/C/G/T (ignores N/other), 0.0 if no ACGT.
    pub fn gc(&self) -> f64 {
        let gc = self.g + self.c;
        let acgt = self.a + self.c + self.g + self.t;
        if acgt == 0 {
            0.0
        } else {
            gc as f64 / acgt as f64
        }
    }
    pub fn add(&mut self, o: &Counts) {
        self.a += o.a;
        self.c += o.c;
        self.g += o.g;
        self.t += o.t;
        self.n += o.n;
        self.other += o.other;
    }
}

/// Count bases over a raw byte buffer. This is exactly the per-byte mapping the
/// GPU kernel performs, kept here as the single source of truth.
pub fn count_bytes(buf: &[u8]) -> Counts {
    let mut c = Counts::default();
    for &b in buf {
        match b {
            b'A' | b'a' => c.a += 1,
            b'C' | b'c' => c.c += 1,
            b'G' | b'g' => c.g += 1,
            b'T' | b't' => c.t += 1,
            b'N' | b'n' => c.n += 1,
            b'\n' | b'\r' => {} // line breaks are not bases
            _ => c.other += 1,
        }
    }
    c
}

/// Stream a FASTA/FASTQ file and tally base composition over *sequence* bytes
/// only (headers, `+` separators and quality lines are skipped).
pub fn count_file(path: &Path, fmt: Format) -> Result<Counts> {
    use std::io::BufRead;
    let (mut reader, _gz) = crate::open_reader(path)?;
    let mut total = Counts::default();
    let mut line = Vec::with_capacity(256);
    let mut lineno: u64 = 0;
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        let is_seq = match fmt {
            Format::Fasta => line.first() != Some(&b'>'),
            Format::Fastq => lineno % 4 == 1, // 0:@id 1:seq 2:+ 3:qual
            _ => true,
        };
        if is_seq {
            total.add(&count_bytes(&line));
        }
        lineno += 1;
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn counts_fasta_sequence_only() {
        let d = std::env::temp_dir().join(format!("cleaver-compute-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        let p = d.join("s.fasta");
        File::create(&p).unwrap().write_all(b">h1\nACGT\nAACC\n>h2\nGGTTNN\n").unwrap();
        let c = count_file(&p, Format::Fasta).unwrap();
        // seq bytes: ACGT + AACC + GGTTNN  =>  A:3 C:3 G:3 T:3 N:2 (headers excluded)
        assert_eq!((c.a, c.c, c.g, c.t, c.n), (3, 3, 3, 3, 2));
        assert_eq!(c.total(), 14);
        assert!((c.gc() - 0.5).abs() < 1e-9); // GC=(C+G)=6 over ACGT=12
    }

    #[test]
    fn counts_fastq_sequence_lines_only() {
        let d = std::env::temp_dir().join(format!("cleaver-compute2-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        let p = d.join("r.fastq");
        File::create(&p).unwrap().write_all(b"@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nJJJJ\n").unwrap();
        let c = count_file(&p, Format::Fastq).unwrap();
        // sequences ACGT + GGCC; quality (IIII/JJJJ) and '+' skipped
        assert_eq!((c.a, c.c, c.g, c.t), (1, 3, 3, 1));
        assert_eq!(c.other, 0);
    }
}
