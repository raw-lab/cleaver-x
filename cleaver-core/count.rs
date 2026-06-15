//! featureCounts / htseq-count-style read counting.
//!
//! For one alignment file, every primary mapped record passing the quality and
//! multi-mapping filters is reduced to its covered reference blocks (split at
//! CIGAR `N` introns), those blocks are tested against the [`Annotation`] index,
//! and the read is assigned to a meta-feature when exactly one is hit. The
//! summary categories mirror featureCounts (`Assigned`,
//! `Unassigned_NoFeatures`, `Unassigned_Ambiguity`, …).

use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::Data;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::Header;

use crate::align::{self, AlignReader};
use crate::annotation::{Annotation, Mode, Strand};
use crate::formats::Format;

/// Counting knobs (featureCounts/htseq flag semantics).
#[derive(Debug, Clone)]
pub struct CountParams {
    /// 0 = unstranded, 1 = stranded, 2 = reverse-stranded.
    pub stranded: u8,
    /// Drop alignments with MAPQ below this (0 = keep all).
    pub min_mapq: u8,
    /// Count multi-mapping reads (NH > 1) instead of discarding them.
    pub count_multimappers: bool,
    /// Count reads overlapping >1 meta-feature against all of them.
    pub allow_multi_overlap: bool,
    /// Overlap-resolution mode.
    pub mode: Mode,
    /// Skip secondary/supplementary alignments (recommended).
    pub primary_only: bool,
    /// PE: require both ends mapped — skip a read whose mate is unmapped (`-B`).
    pub require_both_ends: bool,
    /// PE: exclude chimeric reads (supplementary, or mate on a different
    /// reference) (`-C`).
    pub exclude_chimeric: bool,
    /// PE: only count fragments whose |TLEN| is within [min,max] (`-P`).
    pub check_pe_dist: bool,
    pub min_frag_len: u32, // -d
    pub max_frag_len: u32, // -D
}

impl Default for CountParams {
    fn default() -> Self {
        CountParams {
            stranded: 0,
            min_mapq: 0,
            count_multimappers: false,
            allow_multi_overlap: false,
            mode: Mode::Union,
            primary_only: true,
            require_both_ends: false,
            exclude_chimeric: false,
            check_pe_dist: false,
            min_frag_len: 50,
            max_frag_len: 600,
        }
    }
}

/// Per-file result: a count per meta-feature plus featureCounts-style summary.
#[derive(Debug, Clone)]
pub struct FileCounts {
    pub counts: Vec<u64>,
    pub assigned: u64,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub unmapped: u64,
    pub low_mapq: u64,
    pub multimapping: u64,
    pub pe_filtered: u64,
    pub total: u64,
}

impl FileCounts {
    fn new(n_genes: usize) -> Self {
        FileCounts {
            counts: vec![0; n_genes],
            assigned: 0,
            no_feature: 0,
            ambiguous: 0,
            unmapped: 0,
            low_mapq: 0,
            multimapping: 0,
            pe_filtered: 0,
            total: 0,
        }
    }
    pub fn assignment_rate(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.assigned as f64 / self.total as f64 * 100.0
        }
    }
}

enum Pre {
    Unmapped,
    LowMapq,
    Multi,
    PeFiltered,
    /// Passed all read-level filters; carries the data needed to assign.
    Ok { name: String, blocks: Vec<(u64, u64)>, strand: Strand },
}

/// Covered reference intervals (1-based inclusive). `M/=/X/D` extend the current
/// block; `N` (intron) closes it and starts a new one after the gap; `I/S/H/P`
/// do not consume the reference.
fn ref_blocks(start: u64, ops: &[Op]) -> Vec<(u64, u64)> {
    let mut blocks = Vec::new();
    let mut pos = start;
    let mut cur_start = start;
    let mut in_block = false;
    for op in ops {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                if !in_block {
                    cur_start = pos;
                    in_block = true;
                }
                pos += op.len() as u64;
            }
            Kind::Skip => {
                if in_block {
                    blocks.push((cur_start, pos - 1));
                    in_block = false;
                }
                pos += op.len() as u64;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if in_block {
        blocks.push((cur_start, pos - 1));
    }
    blocks
}

fn nh_tag(d: &Data) -> Option<i64> {
    match d.get(&Tag::ALIGNMENT_HIT_COUNT) {
        Some(Value::Int8(v)) => Some(*v as i64),
        Some(Value::UInt8(v)) => Some(*v as i64),
        Some(Value::Int16(v)) => Some(*v as i64),
        Some(Value::UInt16(v)) => Some(*v as i64),
        Some(Value::Int32(v)) => Some(*v as i64),
        Some(Value::UInt32(v)) => Some(*v as i64),
        _ => None,
    }
}

/// Read-level filtering, independent of any annotation. Returns the read's
/// reference name, covered blocks, and strand when it passes.
fn prefilter(rec: &RecordBuf, header: &Header, p: &CountParams) -> Pre {
    let flags = rec.flags();
    if p.primary_only && (flags.is_secondary() || flags.is_supplementary()) {
        return Pre::Multi;
    }
    if p.exclude_chimeric && flags.is_supplementary() {
        return Pre::PeFiltered;
    }
    if flags.is_unmapped() {
        return Pre::Unmapped;
    }
    // Paired-end filters (only meaningful for segmented reads).
    if flags.is_segmented() {
        if p.require_both_ends && flags.is_mate_unmapped() {
            return Pre::PeFiltered;
        }
        if p.exclude_chimeric {
            if let (Some(a), Some(b)) = (rec.reference_sequence_id(), rec.mate_reference_sequence_id()) {
                if a != b {
                    return Pre::PeFiltered;
                }
            }
        }
        if p.check_pe_dist {
            let tlen = rec.template_length().unsigned_abs();
            if tlen != 0 && (tlen < p.min_frag_len || tlen > p.max_frag_len) {
                return Pre::PeFiltered;
            }
        }
    }
    let mapq = rec.mapping_quality().map(|m| m.get()).unwrap_or(0);
    if (mapq as u16) < p.min_mapq as u16 {
        return Pre::LowMapq;
    }
    if !p.count_multimappers {
        if let Some(nh) = nh_tag(rec.data()) {
            if nh > 1 {
                return Pre::Multi;
            }
        }
    }
    let rid = match rec.reference_sequence_id() {
        Some(i) => i,
        None => return Pre::Unmapped,
    };
    let name: &[u8] = match header.reference_sequences().get_index(rid) {
        Some((n, _)) => n.as_ref(),
        None => return Pre::Ok { name: String::new(), blocks: Vec::new(), strand: Strand::Fwd },
    };
    let name = match std::str::from_utf8(name) {
        Ok(s) => s.to_string(),
        Err(_) => return Pre::Ok { name: String::new(), blocks: Vec::new(), strand: Strand::Fwd },
    };
    let start = match rec.alignment_start() {
        Some(pos) => pos.get() as u64,
        None => return Pre::Unmapped,
    };
    let blocks = ref_blocks(start, rec.cigar().as_ref());
    let strand = if flags.is_reverse_complemented() { Strand::Rev } else { Strand::Fwd };
    Pre::Ok { name, blocks, strand }
}

/// Fold a distinct-gene assignment into a [`FileCounts`].
fn tally_assignment(fc: &mut FileCounts, genes: &[u32], allow_multi: bool) {
    match genes.len() {
        0 => fc.no_feature += 1,
        1 => {
            fc.assigned += 1;
            fc.counts[genes[0] as usize] += 1;
        }
        _ => {
            if allow_multi {
                fc.assigned += 1;
                for &g in genes {
                    fc.counts[g as usize] += 1;
                }
            } else {
                fc.ambiguous += 1;
            }
        }
    }
}

/// Count one SAM/BAM file against `ann`.
pub fn count_file(path: &Path, fmt: Format, ann: &Annotation, params: &CountParams) -> Result<FileCounts> {
    let mut reader = align::open_align_reader(path, fmt)?;
    let header = align::read_align_header(&mut reader)?;
    let mut fc = FileCounts::new(ann.n_genes());

    macro_rules! go {
        ($r:expr) => {{
            for result in $r.record_bufs(&header) {
                let rec = result.context("reading alignment record")?;
                fc.total += 1;
                match prefilter(&rec, &header, params) {
                    Pre::Unmapped => fc.unmapped += 1,
                    Pre::LowMapq => fc.low_mapq += 1,
                    Pre::Multi => fc.multimapping += 1,
                    Pre::PeFiltered => fc.pe_filtered += 1,
                    Pre::Ok { name, blocks, strand } => {
                        if blocks.is_empty() {
                            fc.no_feature += 1;
                        } else {
                            let genes = ann.assign(&name, &blocks, strand, params.stranded, params.mode);
                            tally_assignment(&mut fc, &genes, params.allow_multi_overlap);
                        }
                    }
                }
            }
        }};
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    Ok(fc)
}

/// Hierarchical multi-feature-type counting (VERSE `-t a;b;c` priority order):
/// each read is assigned to the first feature type that yields a unique
/// meta-feature. Returns per-type gene counts plus a shared summary.
#[derive(Debug, Clone)]
pub struct HierCounts {
    /// One gene-count vector per annotation, in priority order.
    pub per_type: Vec<Vec<u64>>,
    /// Reads uniquely assigned at each level.
    pub assigned_per_type: Vec<u64>,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub unmapped: u64,
    pub low_mapq: u64,
    pub multimapping: u64,
    pub pe_filtered: u64,
    pub total: u64,
}

pub fn count_file_hier(
    path: &Path,
    fmt: Format,
    anns: &[Annotation],
    params: &CountParams,
) -> Result<HierCounts> {
    let mut reader = align::open_align_reader(path, fmt)?;
    let header = align::read_align_header(&mut reader)?;
    let mut hc = HierCounts {
        per_type: anns.iter().map(|a| vec![0u64; a.n_genes()]).collect(),
        assigned_per_type: vec![0u64; anns.len()],
        no_feature: 0,
        ambiguous: 0,
        unmapped: 0,
        low_mapq: 0,
        multimapping: 0,
        pe_filtered: 0,
        total: 0,
    };

    macro_rules! go {
        ($r:expr) => {{
            for result in $r.record_bufs(&header) {
                let rec = result.context("reading alignment record")?;
                hc.total += 1;
                match prefilter(&rec, &header, params) {
                    Pre::Unmapped => hc.unmapped += 1,
                    Pre::LowMapq => hc.low_mapq += 1,
                    Pre::Multi => hc.multimapping += 1,
                    Pre::PeFiltered => hc.pe_filtered += 1,
                    Pre::Ok { name, blocks, strand } => {
                        if blocks.is_empty() {
                            hc.no_feature += 1;
                            continue;
                        }
                        // try each feature type in priority order
                        let mut placed = false;
                        let mut saw_ambiguous = false;
                        for (i, ann) in anns.iter().enumerate() {
                            let genes = ann.assign(&name, &blocks, strand, params.stranded, params.mode);
                            match genes.len() {
                                0 => {}
                                1 => {
                                    hc.assigned_per_type[i] += 1;
                                    hc.per_type[i][genes[0] as usize] += 1;
                                    placed = true;
                                    break;
                                }
                                _ => {
                                    if params.allow_multi_overlap {
                                        hc.assigned_per_type[i] += 1;
                                        for g in genes {
                                            hc.per_type[i][g as usize] += 1;
                                        }
                                        placed = true;
                                        break;
                                    } else {
                                        saw_ambiguous = true;
                                    }
                                }
                            }
                        }
                        if !placed {
                            if saw_ambiguous {
                                hc.ambiguous += 1;
                            } else {
                                hc.no_feature += 1;
                            }
                        }
                    }
                }
            }
        }};
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    Ok(hc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn workdir() -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("cleaver-count-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        d
    }

    fn write_gtf(dir: &Path) -> std::path::PathBuf {
        // geneX 100-300, geneY 200-400 (overlap 200-300), geneZ 1000-1100 — all '+'
        let p = dir.join("c.gtf");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "chr1\ts\texon\t100\t300\t.\t+\t.\tgene_id \"geneX\";").unwrap();
        writeln!(f, "chr1\ts\texon\t200\t400\t.\t+\t.\tgene_id \"geneY\";").unwrap();
        writeln!(f, "chr1\ts\texon\t1000\t1100\t.\t+\t.\tgene_id \"geneZ\";").unwrap();
        p
    }

    fn write_sam(dir: &Path) -> std::path::PathBuf {
        let p = dir.join("c.sam");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "@HD\tVN:1.6\tSO:coordinate").unwrap();
        writeln!(f, "@SQ\tSN:chr1\tLN:2000").unwrap();
        // rA 120-169 -> geneX only
        writeln!(f, "rA\t0\tchr1\t120\t60\t50M\t*\t0\t0\t*\t*").unwrap();
        // rB 250-269 -> geneX + geneY -> ambiguous
        writeln!(f, "rB\t0\tchr1\t250\t60\t20M\t*\t0\t0\t*\t*").unwrap();
        // rC 500-529 -> intergenic -> no_feature
        writeln!(f, "rC\t0\tchr1\t500\t60\t30M\t*\t0\t0\t*\t*").unwrap();
        // rD 350-369 -> geneY only
        writeln!(f, "rD\t0\tchr1\t350\t60\t20M\t*\t0\t0\t*\t*").unwrap();
        // rE unmapped
        writeln!(f, "rE\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*").unwrap();
        // rF multimapper (NH:i:3)
        writeln!(f, "rF\t0\tchr1\t1000\t60\t20M\t*\t0\t0\t*\t*\tNH:i:3").unwrap();
        p
    }

    #[test]
    fn counts_and_summary_sam() {
        let d = workdir();
        let gtf = write_gtf(&d);
        let sam = write_sam(&d);
        let ann = Annotation::from_path(&gtf, "exon", "gene_id").unwrap();
        let fc = count_file(&sam, Format::Sam, &ann, &CountParams::default()).unwrap();

        assert_eq!(fc.counts[0], 1, "geneX"); // rA
        assert_eq!(fc.counts[1], 1, "geneY"); // rD
        assert_eq!(fc.counts[2], 0, "geneZ");
        assert_eq!(fc.assigned, 2);
        assert_eq!(fc.ambiguous, 1); // rB
        assert_eq!(fc.no_feature, 1); // rC
        assert_eq!(fc.unmapped, 1); // rE
        assert_eq!(fc.multimapping, 1); // rF
        assert_eq!(fc.total, 6);
    }

    #[test]
    fn allow_multi_overlap_assigns_both() {
        let d = workdir();
        let gtf = write_gtf(&d);
        let sam = write_sam(&d);
        let ann = Annotation::from_path(&gtf, "exon", "gene_id").unwrap();
        let p = CountParams { allow_multi_overlap: true, ..Default::default() };
        let fc = count_file(&sam, Format::Sam, &ann, &p).unwrap();
        // rB now counts for BOTH geneX and geneY
        assert_eq!(fc.counts[0], 2, "geneX gets rA + rB");
        assert_eq!(fc.counts[1], 2, "geneY gets rD + rB");
        assert_eq!(fc.ambiguous, 0);
        assert_eq!(fc.assigned, 3);
    }

    #[test]
    fn min_mapq_filters() {
        let d = workdir();
        let gtf = write_gtf(&d);
        let p = d.join("lowq.sam");
        {
            let mut f = std::fs::File::create(&p).unwrap();
            writeln!(f, "@HD\tVN:1.6").unwrap();
            writeln!(f, "@SQ\tSN:chr1\tLN:2000").unwrap();
            writeln!(f, "q1\t0\tchr1\t120\t5\t30M\t*\t0\t0\t*\t*").unwrap(); // mapq 5
        }
        let ann = Annotation::from_path(&gtf, "exon", "gene_id").unwrap();
        let p10 = CountParams { min_mapq: 10, ..Default::default() };
        let fc = count_file(&p, Format::Sam, &ann, &p10).unwrap();
        assert_eq!(fc.low_mapq, 1);
        assert_eq!(fc.assigned, 0);
    }

    #[test]
    fn ref_blocks_splits_on_intron() {
        let ops = [
            Op::new(Kind::Match, 51),
            Op::new(Kind::Skip, 99),
            Op::new(Kind::Match, 51),
        ];
        let b = ref_blocks(150, &ops);
        assert_eq!(b, vec![(150, 200), (300, 350)]);
    }
}
