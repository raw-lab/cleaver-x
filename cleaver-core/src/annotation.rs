//! Annotation model for `cleaver count`: a GTF/GFF3 parser plus a binned
//! interval index and overlap query supporting htseq-style modes.
//!
//! Features of the chosen type (e.g. `exon`) are grouped into meta-features by
//! an attribute key (e.g. `gene_id`). A read's covered reference blocks are
//! tested against the index; the resulting set of meta-feature ids drives the
//! assignment (unique → counted, none → no-feature, many → ambiguous).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

/// Binning granularity for the per-reference interval index (bases).
const BIN: u64 = 16_384;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Strand {
    Fwd,
    Rev,
    Unknown,
}

/// htseq overlap-resolution mode.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    /// Any feature overlapping any base of the read (htseq `union`,
    /// featureCounts default).
    Union,
    /// Meta-features whose features cover *every* base of the read.
    IntersectionStrict,
    /// As strict, but bases covered by no feature are ignored.
    IntersectionNonempty,
}

struct Feat {
    start: u64, // 1-based inclusive
    end: u64,
    strand: Strand,
    gene: u32,
}

#[derive(Default)]
struct Chrom {
    feats: Vec<Feat>,
    bins: HashMap<u64, Vec<u32>>,
}

pub struct Annotation {
    pub feature_type: String,
    pub group_key: String,
    gene_ids: Vec<String>,
    gene_index: HashMap<String, u32>,
    chroms: HashMap<String, Chrom>,
    /// Feature lines of the requested type that were indexed.
    pub n_features: u64,
}

/// Extract the value for `key` from a GTF (`key "value"`) or GFF3 (`key=value`)
/// attribute column.
fn attr_value(attrs: &str, key: &str) -> Option<String> {
    for field in attrs.split(';') {
        let f = field.trim();
        if !f.starts_with(key) {
            continue;
        }
        let after = f.as_bytes().get(key.len()).copied();
        let rest = f[key.len()..].trim_start();
        // GFF3: key=value
        if let Some(v) = rest.strip_prefix('=') {
            return Some(v.trim().trim_matches('"').to_string());
        }
        // GTF: key "value" (a space separated the key from the value)
        if after == Some(b' ') || after == Some(b'\t') {
            let v = rest.trim().trim_matches('"');
            if !v.is_empty() {
                return Some(v.to_string());
            }
        }
        // otherwise this was a different key sharing a prefix; keep scanning
    }
    None
}

fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Fwd,
        "-" => Strand::Rev,
        _ => Strand::Unknown,
    }
}

impl Annotation {
    pub fn n_genes(&self) -> usize {
        self.gene_ids.len()
    }
    pub fn gene_id(&self, i: u32) -> &str {
        &self.gene_ids[i as usize]
    }
    pub fn gene_ids(&self) -> &[String] {
        &self.gene_ids
    }

    /// Parse a GTF/GFF file, indexing only rows whose column 3 equals
    /// `feature_type` and that carry the `group_key` attribute.
    pub fn from_path(path: &Path, feature_type: &str, group_key: &str) -> Result<Annotation> {
        let f = File::open(path).with_context(|| format!("opening annotation '{}'", path.display()))?;
        let mut reader = BufReader::new(f);

        let mut ann = Annotation {
            feature_type: feature_type.to_string(),
            group_key: group_key.to_string(),
            gene_ids: Vec::new(),
            gene_index: HashMap::new(),
            chroms: HashMap::new(),
            n_features: 0,
        };

        let mut line = String::new();
        loop {
            line.clear();
            let n = reader.read_line(&mut line).context("reading annotation")?;
            if n == 0 {
                break;
            }
            let l = line.trim_end_matches(['\n', '\r']);
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let mut cols = l.split('\t');
            let seqname = match cols.next() {
                Some(s) => s,
                None => continue,
            };
            let _source = cols.next();
            let ftype = match cols.next() {
                Some(s) => s,
                None => continue,
            };
            if ftype != feature_type {
                continue;
            }
            let start: u64 = match cols.next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let end: u64 = match cols.next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let _score = cols.next();
            let strand = parse_strand(cols.next().unwrap_or("."));
            let _frame = cols.next();
            let attrs = cols.next().unwrap_or("");
            let gene_name = match attr_value(attrs, group_key) {
                Some(g) => g,
                None => continue,
            };
            if start == 0 || end < start {
                continue;
            }

            let gene = match ann.gene_index.get(&gene_name) {
                Some(&i) => i,
                None => {
                    let i = ann.gene_ids.len() as u32;
                    ann.gene_ids.push(gene_name.clone());
                    ann.gene_index.insert(gene_name, i);
                    i
                }
            };

            let chrom = ann.chroms.entry(seqname.to_string()).or_default();
            let idx = chrom.feats.len() as u32;
            chrom.feats.push(Feat { start, end, strand, gene });
            for b in (start / BIN)..=(end / BIN) {
                chrom.bins.entry(b).or_default().push(idx);
            }
            ann.n_features += 1;
        }
        Ok(ann)
    }

    /// Distinct meta-feature ids assigned to a read, given its covered reference
    /// blocks (1-based inclusive, sorted), read strand, strandedness (0/1/2) and
    /// overlap mode. Empty → no feature; len 1 → unique; len > 1 → ambiguous.
    pub fn assign(
        &self,
        seqname: &str,
        blocks: &[(u64, u64)],
        read_strand: Strand,
        stranded: u8,
        mode: Mode,
    ) -> Vec<u32> {
        let chrom = match self.chroms.get(seqname) {
            Some(c) => c,
            None => return Vec::new(),
        };

        // Candidate feature indices from the bins the read touches.
        let mut cidx: Vec<u32> = Vec::new();
        for &(bs, be) in blocks {
            for b in (bs / BIN)..=(be / BIN) {
                if let Some(v) = chrom.bins.get(&b) {
                    cidx.extend_from_slice(v);
                }
            }
        }
        cidx.sort_unstable();
        cidx.dedup();

        // Strand-filter and keep those overlapping at least one block.
        let mut cand: Vec<(u64, u64, u32)> = Vec::new();
        for &i in &cidx {
            let f = &chrom.feats[i as usize];
            if !strand_ok(f.strand, read_strand, stranded) {
                continue;
            }
            if blocks.iter().any(|&(bs, be)| f.start <= be && f.end >= bs) {
                cand.push((f.start, f.end, f.gene));
            }
        }
        if cand.is_empty() {
            return Vec::new();
        }

        match mode {
            Mode::Union => {
                let mut g: Vec<u32> = cand.iter().map(|c| c.2).collect();
                g.sort_unstable();
                g.dedup();
                g
            }
            Mode::IntersectionStrict | Mode::IntersectionNonempty => sweep(blocks, &cand, mode),
        }
    }
}

fn strand_ok(feat: Strand, read: Strand, stranded: u8) -> bool {
    match stranded {
        0 => true,
        1 => feat == Strand::Unknown || feat == read,
        2 => feat == Strand::Unknown || feat == opposite(read),
        _ => true,
    }
}

fn opposite(s: Strand) -> Strand {
    match s {
        Strand::Fwd => Strand::Rev,
        Strand::Rev => Strand::Fwd,
        Strand::Unknown => Strand::Unknown,
    }
}

fn intersect(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

/// Per-base intersection sweep for the strict / nonempty modes. Segments tile
/// each read block at feature boundaries; the gene set of a segment is the set
/// of meta-features whose feature fully covers it.
fn sweep(blocks: &[(u64, u64)], cand: &[(u64, u64, u32)], mode: Mode) -> Vec<u32> {
    let mut result: Option<Vec<u32>> = None;
    let mut any_nonempty = false;

    for &(bs, be) in blocks {
        let mut pts: Vec<u64> = vec![bs, be + 1];
        for &(fs, fe, _) in cand {
            let s = fs.max(bs);
            let e = fe.min(be);
            if s <= e {
                pts.push(s);
                pts.push(e + 1);
            }
        }
        pts.sort_unstable();
        pts.dedup();

        for w in pts.windows(2) {
            let (x, y) = (w[0], w[1]);
            if y <= x {
                continue;
            }
            let last = y - 1; // segment [x, last]
            if x < bs || last > be {
                continue;
            }
            let mut gs: Vec<u32> = Vec::new();
            for &(fs, fe, g) in cand {
                if fs <= x && fe >= last {
                    gs.push(g);
                }
            }
            gs.sort_unstable();
            gs.dedup();

            if gs.is_empty() {
                if mode == Mode::IntersectionStrict {
                    return Vec::new(); // a base with no feature kills strict
                }
                continue; // nonempty: ignore uncovered bases
            }
            any_nonempty = true;
            result = Some(match result.take() {
                None => gs,
                Some(prev) => intersect(&prev, &gs),
            });
            if result.as_ref().map(|r| r.is_empty()).unwrap_or(false) {
                return Vec::new();
            }
        }
    }
    match result {
        Some(r) if any_nonempty => r,
        _ => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_gtf() -> std::path::PathBuf {
        // geneA exons on '+' at 100-200 and 300-400; geneB on '-' at 1000-1100.
        let d = std::env::temp_dir().join(format!("cleaver-ann-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        let p = d.join("a.gtf");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "# comment").unwrap();
        writeln!(f, "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"geneA\"; transcript_id \"tA\";").unwrap();
        writeln!(f, "chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id \"geneA\";").unwrap();
        writeln!(f, "chr1\tsrc\texon\t1000\t1100\t.\t-\t.\tgene_id \"geneB\";").unwrap();
        writeln!(f, "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tgene_id \"geneA\";").unwrap(); // wrong type, ignored
        p
    }

    #[test]
    fn parse_and_index() {
        let p = write_gtf();
        let a = Annotation::from_path(&p, "exon", "gene_id").unwrap();
        assert_eq!(a.n_genes(), 2);
        assert_eq!(a.n_features, 3); // gene row excluded
        assert_eq!(a.gene_id(0), "geneA");
        assert_eq!(a.gene_id(1), "geneB");
    }

    #[test]
    fn union_overlap_unique_ambiguous_none() {
        let p = write_gtf();
        let a = Annotation::from_path(&p, "exon", "gene_id").unwrap();
        // read inside geneA exon1
        assert_eq!(a.assign("chr1", &[(120, 180)], Strand::Fwd, 0, Mode::Union), vec![0]);
        // read in intergenic space
        assert!(a.assign("chr1", &[(500, 600)], Strand::Fwd, 0, Mode::Union).is_empty());
        // unknown chromosome
        assert!(a.assign("chrX", &[(120, 180)], Strand::Fwd, 0, Mode::Union).is_empty());
        // spliced read across geneA exon1 and exon2 (gap = intron) -> still geneA
        assert_eq!(a.assign("chr1", &[(150, 200), (300, 350)], Strand::Fwd, 0, Mode::Union), vec![0]);
    }

    #[test]
    fn strandedness() {
        let p = write_gtf();
        let a = Annotation::from_path(&p, "exon", "gene_id").unwrap();
        // geneB is on '-'; a forward read with -s 1 should NOT match it
        assert!(a.assign("chr1", &[(1020, 1080)], Strand::Fwd, 1, Mode::Union).is_empty());
        // a reverse read with -s 1 matches geneB
        assert_eq!(a.assign("chr1", &[(1020, 1080)], Strand::Rev, 1, Mode::Union), vec![1]);
        // unstranded matches regardless
        assert_eq!(a.assign("chr1", &[(1020, 1080)], Strand::Fwd, 0, Mode::Union), vec![1]);
    }

    #[test]
    fn strict_vs_nonempty() {
        let p = write_gtf();
        let a = Annotation::from_path(&p, "exon", "gene_id").unwrap();
        // read 150-250: 150-200 in geneA exon1, 201-250 intergenic
        // union -> {geneA}; strict -> {} (uncovered bases); nonempty -> {geneA}
        assert_eq!(a.assign("chr1", &[(150, 250)], Strand::Fwd, 0, Mode::Union), vec![0]);
        assert!(a.assign("chr1", &[(150, 250)], Strand::Fwd, 0, Mode::IntersectionStrict).is_empty());
        assert_eq!(
            a.assign("chr1", &[(150, 250)], Strand::Fwd, 0, Mode::IntersectionNonempty),
            vec![0]
        );
        // fully inside exon1 -> strict also assigns
        assert_eq!(a.assign("chr1", &[(120, 180)], Strand::Fwd, 0, Mode::IntersectionStrict), vec![0]);
    }
}
