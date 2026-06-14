/*!
Render plots from the CSV files written by `rustyomestats`, and bundle them
into a single self-contained HTML report.

Inputs : CSV files produced by `rustyomestats genome` and/or `rustyomestats u50`.
Outputs: one PNG per figure, plus `report.html` with every PNG embedded
         (base64) so the HTML can be shared as a single file.

Stack  : csv (read CSV), plotters (plots), base64, clap, chrono.
*/

use anyhow::{Context, Result};
use base64::{engine::general_purpose::STANDARD as B64, Engine};
use chrono::Local;
use clap::Parser;
use plotters::prelude::*;
use serde::Deserialize;
use std::{
    collections::HashMap,
    fs,
//    io::Write,
    path::{Path, PathBuf},
};

// ── frame ordering ────────────────────────────────────────────────────────────
const FRAME_ORDER: [&str; 6] = ["+1", "+2", "+3", "-1", "-2", "-3"];

// ── colour palette ────────────────────────────────────────────────────────────
const STEEL_BLUE: RGBColor = RGBColor(70, 130, 180);
const SEA_GREEN: RGBColor = RGBColor(46, 139, 87);
const FIREBRICK: RGBColor = RGBColor(178, 34, 34);
const PURPLE: RGBColor = RGBColor(129, 114, 178);
const NAVY: RGBColor = RGBColor(76, 114, 176);
const GREEN: RGBColor = RGBColor(85, 168, 104);

// ── PNG dimensions ────────────────────────────────────────────────────────────
const W_WIDE: u32 = 1400;
const W_STD: u32 = 1000;
const H_STD: u32 = 600;
const H_SHORT: u32 = 380;

// ══════════════════════════════════════════════════════════════════════════════
// CSV row types
// ══════════════════════════════════════════════════════════════════════════════

#[derive(Debug, Deserialize)]
struct PerSequenceRow {
    length: f64,
    gc: f64,
}

#[derive(Debug, Deserialize)]
struct CodonAbsRow {
    frame: String,
    codon: String,
    count: f64,
}

#[derive(Debug, Deserialize)]
struct CodonAggRow {
    codon: String,
    density: f64,
}

#[derive(Debug, Deserialize)]
struct CodonCompRow {
    codon: String,
    enrichment_pred_over_abs: f64,
}

#[derive(Debug, Deserialize)]
struct U50ContigRow {
    orig_length: f64,
    unique_bp: f64,
}

#[derive(Debug, Deserialize)]
struct GapIntervalRow {
    start: usize,
    end: usize,
}

// ══════════════════════════════════════════════════════════════════════════════
// generic CSV helpers
// ══════════════════════════════════════════════════════════════════════════════

fn read_csv<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<Vec<T>> {
    let mut rdr = csv::Reader::from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    rdr.deserialize()
        .map(|r| r.context("deserialising CSV row"))
        .collect()
}

fn read_summary_map(path: &Path) -> Option<HashMap<String, String>> {
    #[derive(Deserialize)]
    struct Row {
        metric: String,
        value: String,
    }
    let rows: Vec<Row> = read_csv(path).ok()?;
    Some(rows.into_iter().map(|r| (r.metric, r.value)).collect())
}

// ══════════════════════════════════════════════════════════════════════════════
// histogram helpers (bin the data ourselves, plotters draws rectangles)
// ══════════════════════════════════════════════════════════════════════════════

/// Build `n_bins` equal-width histogram bins over `values`.
/// Returns `(bin_edges, counts)` where `bin_edges` has length `n_bins + 1`.
fn histogram(values: &[f64], n_bins: usize) -> (Vec<f64>, Vec<u64>) {
    if values.is_empty() || n_bins == 0 {
        return (vec![], vec![]);
    }
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = (max - min).max(f64::EPSILON);
    let width = range / n_bins as f64;

    let mut counts = vec![0u64; n_bins];
    for &v in values {
        let idx = ((v - min) / width) as usize;
        let idx = idx.min(n_bins - 1);
        counts[idx] += 1;
    }
    let edges: Vec<f64> = (0..=n_bins).map(|i| min + i as f64 * width).collect();
    (edges, counts)
}

/// Same as `histogram` but on `log10(v)` for positive values, returns
/// linear-scale edges (i.e. powers of 10) and counts.
fn log_histogram(values: &[f64], n_bins: usize) -> (Vec<f64>, Vec<u64>) {
    let logs: Vec<f64> = values.iter().filter(|&&v| v > 0.0).map(|v| v.log10()).collect();
    let (log_edges, counts) = histogram(&logs, n_bins);
    let lin_edges: Vec<f64> = log_edges.iter().map(|&e| 10f64.powf(e)).collect();
    (lin_edges, counts)
}

// ══════════════════════════════════════════════════════════════════════════════
// save helper
// ══════════════════════════════════════════════════════════════════════════════

fn save_path(outdir: &Path, name: &str) -> PathBuf {
    outdir.join(format!("{name}.png"))
}

// ══════════════════════════════════════════════════════════════════════════════
// superscript helper
// ══════════════════════════════════════════════════════════════════════════════

fn superscript(n: i32) -> String {
    n.to_string()
        .chars()
        .map(|c| match c {
            '0' => '⁰',
            '1' => '¹',
            '2' => '²',
            '3' => '³',
            '4' => '⁴',
            '5' => '⁵',
            '6' => '⁶',
            '7' => '⁷',
            '8' => '⁸',
            '9' => '⁹',
            '-' => '⁻',
            _ => c,
        })
        .collect()
}

// ══════════════════════════════════════════════════════════════════════════════
// per-sequence plots
// ══════════════════════════════════════════════════════════════════════════════

pub fn plot_length_hist(outdir: &Path) -> Option<PathBuf> {
    let p = outdir.join("per_sequence.csv");
    if !p.exists() {
        return None;
    }
    let rows: Vec<PerSequenceRow> = read_csv(&p).ok()?;
    let lengths: Vec<f64> = rows.iter().map(|r| r.length).collect();

    let out = save_path(outdir, "plot_length_histogram");
    let root = BitMapBackend::new(&out, (W_STD, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let (log_edges, counts) = log_histogram(&lengths, 50);
    if log_edges.len() < 2 {
        return None;
    }
    let x_min = log_edges[0];
    let x_max = *log_edges.last().unwrap();
    let y_max = *counts.iter().max().unwrap_or(&1) as f64 * 1.1;

    let mut chart = ChartBuilder::on(&root)
        .caption("Sequence length distribution", ("sans-serif", 22))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d((x_min..x_max).log_scale(), 0f64..y_max)
        .ok()?;

    chart
        .configure_mesh()
        .max_light_lines(0)
        .x_desc("length (bp, log10)")
        .y_desc("# sequences")
        .x_label_formatter(&|v| {
            let exp = v.log10().round() as i32;
            format!("10{}", superscript(exp))})
        .draw()
        .ok()?;

    for i in 0..counts.len() {
        let x0 = log_edges[i];
        let x1 = log_edges[i + 1];
        let h = counts[i] as f64;
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(x0, 0f64), (x1, h)],
                STEEL_BLUE.mix(0.8).filled(),
            )))
            .ok()?;
    }

    root.present().ok()?;
    Some(out.clone())
}

pub fn plot_gc_distribution(outdir: &Path) -> Option<PathBuf> {
    let p = outdir.join("per_sequence.csv");
    if !p.exists() {
        return None;
    }
    let rows: Vec<PerSequenceRow> = read_csv(&p).ok()?;
    let gc_vals: Vec<f64> = rows.iter().map(|r| r.gc).collect();

    let out = save_path(outdir, "plot_gc_distribution");
    let root = BitMapBackend::new(&out, (W_STD, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let (edges, counts) = histogram(&gc_vals, 40);
    if edges.len() < 2 {
        return None;
    }
    let y_max = *counts.iter().max().unwrap_or(&1) as f64 * 1.15;

    let mut chart = ChartBuilder::on(&root)
        .caption("GC% per sequence", ("sans-serif", 22))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(edges[0]..(*edges.last().unwrap()), 0f64..y_max)
        .ok()?;

    chart
        .configure_mesh()
        .x_desc("GC %")
        .y_desc("count")
		.max_light_lines(0)
        .draw()
        .ok()?;

    for i in 0..counts.len() {
        let x0 = edges[i];
        let x1 = edges[i + 1];
        let h = counts[i] as f64;
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(x0, 0f64), (x1, h)],
                SEA_GREEN.mix(0.8).filled(),
            )))
            .ok()?;
    }

    root.present().ok()?;
    Some(out.clone())
}

pub fn plot_gc_vs_length(outdir: &Path) -> Option<PathBuf> {
    let p = outdir.join("per_sequence.csv");
    if !p.exists() {
        return None;
    }
    let rows: Vec<PerSequenceRow> = read_csv(&p).ok()?;
    if rows.is_empty() {
        return None;
    }

    let len_min = rows.iter().map(|r| r.length).fold(f64::INFINITY, f64::min);
    let len_max = rows.iter().map(|r| r.length).fold(f64::NEG_INFINITY, f64::max);
    let gc_min = rows.iter().map(|r| r.gc).fold(f64::INFINITY, f64::min);
    let gc_max = rows.iter().map(|r| r.gc).fold(f64::NEG_INFINITY, f64::max);

    let x_min = len_min / 1.2;
    let x_max = len_max * 1.2;

    let out = save_path(outdir, "plot_gc_vs_length");
    let root = BitMapBackend::new(&out, (W_STD, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let mut chart = ChartBuilder::on(&root)
        .caption("GC% vs sequence length", ("sans-serif", 22))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (x_min..x_max).log_scale(),
            (gc_min - 1.0)..(gc_max + 1.0),
        )
        .ok()?;

    chart
        .configure_mesh()
        .x_desc("length (bp)")
        .y_desc("GC %")
        .x_label_formatter(&|v| {
            let exp = v.log10().round() as i32;
            format!("10{}", superscript(exp))})
		.max_light_lines(0)
        .draw()
        .ok()?;

    chart
        .draw_series(rows.iter().map(|r| {
            Circle::new((r.length, r.gc), 2, NAVY.mix(0.45).filled())
        }))
        .ok()?;

    root.present().ok()?;
    Some(out.clone())
}

// ══════════════════════════════════════════════════════════════════════════════
// codon plots
// ══════════════════════════════════════════════════════════════════════════════

pub fn plot_codon_usage_bar(outdir: &Path) -> Option<PathBuf> {
    let abs_p = outdir.join("codon_absolute_aggregate.csv");
    let pred_p = outdir.join("codon_predicted_aggregate.csv");

    let mut abs_map: HashMap<String, f64> = HashMap::new();
    let mut pred_map: HashMap<String, f64> = HashMap::new();
    let mut codons: Vec<String> = Vec::new();

    if abs_p.exists() {
        let rows: Vec<CodonAggRow> = read_csv(&abs_p).ok()?;
        for r in rows {
            if !codons.contains(&r.codon) {
                codons.push(r.codon.clone());
            }
            abs_map.insert(r.codon, r.density);
        }
    }
    if pred_p.exists() {
        let rows: Vec<CodonAggRow> = read_csv(&pred_p).ok()?;
        for r in rows {
            if !codons.contains(&r.codon) {
                codons.push(r.codon.clone());
            }
            pred_map.insert(r.codon, r.density);
        }
    }
    if codons.is_empty() {
        return None;
    }
    codons.sort();

    let y_max = abs_map
        .values()
        .chain(pred_map.values())
        .cloned()
        .fold(0f64, f64::max)
        * 1.15;

    let n = codons.len();
    let out = save_path(outdir, "plot_codon_usage_bar");
    let root = BitMapBackend::new(&out, (W_WIDE, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Codon density: 6-frame vs FragGeneScan",
            ("sans-serif", 22),
        )
        .margin(10)
        .x_label_area_size(90)
        .y_label_area_size(60)
        .build_cartesian_2d(0f64..(n as f64 * 2.2), 0f64..y_max)
        .ok()?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .x_labels(0)
        .x_desc("codon")
        .y_desc("density (fraction)")
		.max_light_lines(0)
        .draw()
        .ok()?;

    let bar_w = 0.9f64;
    for (i, codon) in codons.iter().enumerate() {
        let x_base = i as f64 * 2.2;

        // absolute bar
        if let Some(&d) = abs_map.get(codon) {
            chart
                .draw_series(std::iter::once(Rectangle::new(
                    [(x_base, 0f64), (x_base + bar_w, d)],
                    NAVY.mix(0.8).filled(),
                )))
                .ok()?;
        }
        // predicted bar
        if let Some(&d) = pred_map.get(codon) {
            chart
                .draw_series(std::iter::once(Rectangle::new(
                    [(x_base + bar_w, 0f64), (x_base + bar_w * 2.0, d)],
                    GREEN.mix(0.8).filled(),
                )))
                .ok()?;
        }
    }
    // Draw codon labels
    let label_y = -y_max * 0.05;
    let label_style = ("sans-serif", 12)
        .into_font()
        .transform(FontTransform::Rotate90)
        .color(&BLACK);

    for (i, codon) in codons.iter().enumerate() {
        let x_base = i as f64 * 2.2;

        // Center of the pair of bars
        let x_center = x_base + bar_w;

        chart
            .draw_series(std::iter::once(Text::new(
                codon.clone(),
                (x_center, label_y),
                label_style.clone(),
            )))
            .ok()?;
    }


    // Draw legend
    let legend_x = n as f64 * 2.2 * 0.78;
    let legend_y = y_max * 0.98;
    let legend_w = 25.0;
    let legend_h = y_max * 0.12;
    let box_w = 0.8;
    let box_h = y_max * 0.03;
    chart
        .draw_series(std::iter::once(Rectangle::new(
            [
                (legend_x, legend_y - legend_h),
                (legend_x + legend_w, legend_y),
            ],
            WHITE.filled(),
        )))
        .ok()?;
    chart
        .draw_series(std::iter::once(Rectangle::new(
            [
                (legend_x, legend_y - legend_h),
                (legend_x + legend_w, legend_y),
            ],
            BLACK.stroke_width(1),
        )))
        .ok()?;
    chart
        .draw_series(std::iter::once(Rectangle::new(
            [
                (legend_x+1.0, legend_y - box_h),
                (legend_x+1.0 + box_w, legend_y),
            ],
            NAVY.mix(0.8).filled(),
        )))
        .ok()?;

    chart
        .draw_series(std::iter::once(Text::new(
            "6-frame background",
            (legend_x + 2.0, legend_y - box_h / 2.0),
            ("sans-serif", 14).into_font(),
        )))
        .ok()?;
    chart
        .draw_series(std::iter::once(Rectangle::new(
            [
                (legend_x+1.0, legend_y - box_h * 2.5),
                (legend_x+1.0 + box_w, legend_y - box_h * 1.5),
            ],
            GREEN.mix(0.8).filled(),
        )))
        .ok()?;

    chart
        .draw_series(std::iter::once(Text::new(
            "FragGeneScan ORFs",
            (legend_x + 2.0, legend_y - box_h * 2.0),
            ("sans-serif", 14).into_font(),
        )))
        .ok()?;

    root.present().ok()?;
    Some(out.clone())
}

pub fn plot_codon_heatmap_by_frame(outdir: &Path) -> Option<PathBuf> {
    let p = outdir.join("codon_absolute.csv");
    if !p.exists() {
        return None;
    }
    let rows: Vec<CodonAbsRow> = read_csv(&p).ok()?;

    // aggregate: (frame, codon) -> total count
    let mut agg: HashMap<(String, String), f64> = HashMap::new();
    let mut frame_totals: HashMap<String, f64> = HashMap::new();
    let mut codon_set: std::collections::BTreeSet<String> = Default::default();

    for r in &rows {
        *agg.entry((r.frame.clone(), r.codon.clone())).or_insert(0.0) += r.count;
        *frame_totals.entry(r.frame.clone()).or_insert(0.0) += r.count;
        codon_set.insert(r.codon.clone());
    }

    let codons: Vec<String> = codon_set.into_iter().collect();
    let n_codons = codons.len();
    let n_frames = FRAME_ORDER.len();

    // density matrix [frame][codon]
    let mut matrix = vec![vec![0f64; n_codons]; n_frames];
    for (fi, frame) in FRAME_ORDER.iter().enumerate() {
        let total = frame_totals.get(*frame).cloned().unwrap_or(1.0).max(1.0);
        for (ci, codon) in codons.iter().enumerate() {
            let count = agg
                .get(&(frame.to_string(), codon.clone()))
                .cloned()
                .unwrap_or(0.0);
            matrix[fi][ci] = count / total;
        }
    }

    let max_density = matrix
        .iter()
        .flat_map(|r| r.iter())
        .cloned()
        .fold(0f64, f64::max);

    let out = save_path(outdir, "plot_codon_heatmap_by_frame");
    let root = BitMapBackend::new(&out, (W_WIDE, 500)).into_drawing_area();
    root.fill(&WHITE).ok()?;
    let (heatmap_area, legend_area) = root.split_horizontally(W_WIDE as i32 - 90);

    let mut chart = ChartBuilder::on(&heatmap_area)
        .caption("6-frame codon density heatmap", ("sans-serif", 22))
        .margin(20)
        .x_label_area_size(80)
        .y_label_area_size(60)
        .build_cartesian_2d(0..n_codons, 0..n_frames)
        .ok()?;

    chart
        .configure_mesh()
        .x_labels(n_codons.min(20))
        .y_labels(n_frames)
        .x_label_formatter(&|v| {
            codons.get(*v).cloned().unwrap_or_default()
        })
        .y_label_formatter(&|v| {
            FRAME_ORDER.get(*v).map(|s| s.to_string()).unwrap_or_default()
        })
        .x_desc("codon")
        .y_desc("frame")
        .draw()
        .ok()?;

    // draw heatmap
    for fi in 0..n_frames {
        for ci in 0..n_codons {
            let v = matrix[fi][ci];
            let intensity = if max_density > 0.0 {
                (v / max_density).clamp(0.0, 1.0)
            } else {
                0.0
            };
            // viridis-like: dark purple → teal → yellow
            let r = (intensity * 253.0) as u8;
            let g = ((0.5 - (intensity - 0.5).abs()) * 2.0 * 231.0) as u8;
            let b = ((1.0 - intensity) * 220.0) as u8;
            let colour = RGBColor(r, g, b).filled();
            chart
                .draw_series(std::iter::once(Rectangle::new(
                    [(ci, fi), (ci + 1, fi + 1)],
                    colour,
                )))
                .ok()?;
        }
    }

    // draw legend
    let (_lw, lh) = legend_area.dim_in_pixel();
    let title_y = 15;
    let bar_top = 35;          // below title
    let bar_bottom = lh as i32 - 20; // leave padding at bottom
    let bar_left = 10;
    let bar_right = 30;

    let bar_height = bar_bottom - bar_top;
    for y in bar_top..bar_bottom {
        let intensity =
            1.0 - ((y - bar_top) as f64 / bar_height as f64);

        let r = (intensity * 253.0) as u8;
        let g = ((0.5 - (intensity - 0.5).abs()) * 2.0 * 231.0) as u8;
        let b = ((1.0 - intensity) * 220.0) as u8;

        legend_area.draw(&Rectangle::new(
            [(bar_left, y), (bar_right, y + 1)],
            RGBColor(r, g, b).filled(),
        )).ok()?;
    }
    for i in 0..=4 {
        let frac = i as f64 / 4.0;

        let y = bar_bottom
            - ((bar_height as f64 * frac) as i32);

        legend_area.draw(&Text::new(
            format!("{:.4}", frac * max_density),
            (bar_right + 8, y),
            ("sans-serif", 10).into_font(),
        )).ok()?;
    }
    legend_area.draw(&Text::new(
        "0",
        (35, lh as i32 - 5),
        ("sans-serif", 12).into_font(),
    )).ok()?;
    legend_area.draw(&Text::new(
        "Density",
        (bar_left, title_y),
        ("sans-serif", 14).into_font(),
    )).ok()?;

    root.present().ok()?;
    Some(out.clone())
}

pub fn plot_codon_enrichment(outdir: &Path) -> Option<PathBuf> {
    let p = outdir.join("codon_comparison.csv");
    if !p.exists() {
        return None;
    }
    let mut rows: Vec<CodonCompRow> = read_csv(&p).ok()?;
    rows.sort_by(|a, b| {
        b.enrichment_pred_over_abs
            .partial_cmp(&a.enrichment_pred_over_abs)
            .unwrap()
    });

    if rows.is_empty() {
        return None;
    }

    let y_min = rows
        .iter()
        .map(|r| r.enrichment_pred_over_abs)
        .fold(f64::INFINITY, f64::min)
        .min(0.8)
        - 0.05;
    let y_max = rows
        .iter()
        .map(|r| r.enrichment_pred_over_abs)
        .fold(f64::NEG_INFINITY, f64::max)
        * 1.1;

    let n = rows.len();
    let out = save_path(outdir, "plot_codon_enrichment");
    let root = BitMapBackend::new(&out, (W_WIDE, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    // Use an integer x-axis (one slot per codon) so we can map each index
    // back to its codon name via x_label_formatter.
    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Codon enrichment (predicted / absolute)",
            ("sans-serif", 22),
        )
        .margin(10)
        .x_label_area_size(90)
        .y_label_area_size(60)
        .build_cartesian_2d(0f64..n as f64, y_min..y_max)
        .ok()?;

    // Bar geometry: each slot is 1.0 wide; bar fills 0.8, gap is 0.2.
    let bar_w = 0.8f64;
    let gap = (1.0 - bar_w) / 2.0; // 0.1 padding on each side

    // Suppress all x-axis tick labels from the mesh — we draw them manually
    // below so every codon is guaranteed to appear without skipping.
    chart
        .configure_mesh()
        .disable_x_mesh()
        .max_light_lines(0)
        .x_labels(0)
        .y_desc("enrichment")
        .x_desc("codons")
        .draw()
        .ok()?;

    // reference line at y = 1.0
    chart
        .draw_series(LineSeries::new(
            [(0f64, 1.0), (n as f64, 1.0)],
            BLACK.stroke_width(1),
        ))
        .ok()?;

    for (i, row) in rows.iter().enumerate() {
        let x0 = i as f64 + gap;
        let x1 = x0 + bar_w;

        // bar
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(x0, 0f64), (x1, row.enrichment_pred_over_abs)],
                FIREBRICK.mix(0.8).filled(),
            )))
            .ok()?;

        // codon label — drawn slightly below y_min so it clears the axis line
        let label_y = y_min - (y_max - y_min) * 0.03;
        let label_style = ("sans-serif", 12)
            .into_font()
            .transform(FontTransform::Rotate90)
            .color(&BLACK);
        chart
            .draw_series(std::iter::once(Text::new(
                row.codon.clone(),
                (i as f64 + 0.5, label_y),
                label_style,
            )))
            .ok()?;
    }

    root.present().ok()?;
    Some(out.clone())
}

// ══════════════════════════════════════════════════════════════════════════════
// U50 plots
// ══════════════════════════════════════════════════════════════════════════════

pub fn plot_u50_contig_lengths(outdir: &Path) -> Option<PathBuf> {
    let p = outdir.join("u50_contigs.csv");
    if !p.exists() {
        return None;
    }
    let rows: Vec<U50ContigRow> = read_csv(&p).ok()?;
    if rows.is_empty() {
        return None;
    }

    let orig: Vec<f64> = rows.iter().filter(|r| r.orig_length > 0.0).map(|r| r.orig_length).collect();
    let unique: Vec<f64> = rows.iter().filter(|r| r.unique_bp > 0.0).map(|r| r.unique_bp).collect();

    let (orig_edges, orig_counts) = log_histogram(&orig, 40);
    let (uniq_edges, uniq_counts) = log_histogram(&unique, 40);
    if orig_edges.len() < 2 && uniq_edges.len() < 2 {
        return None;
    }

    let x_min = orig_edges.first().copied().unwrap_or(1.0)
        .min(uniq_edges.first().copied().unwrap_or(1.0));
    let x_max = orig_edges.last().copied().unwrap_or(1e9)
        .max(uniq_edges.last().copied().unwrap_or(1e9));
    let y_max = orig_counts
        .iter()
        .chain(uniq_counts.iter())
        .max()
        .copied()
        .unwrap_or(1) as f64
        * 1.15;

    let out = save_path(outdir, "plot_u50_contig_lengths");
    let root = BitMapBackend::new(&out, (W_STD, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Contig lengths: original vs unique-bp (post-masking)",
            ("sans-serif", 20),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d((x_min..x_max).log_scale(), 0f64..y_max)
        .ok()?;

    chart
        .configure_mesh()
        .x_desc("length (bp, log10)")
        .y_desc("# contigs")
        .draw()
        .ok()?;

    for i in 0..orig_counts.len() {
        if i + 1 >= orig_edges.len() { break; }
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(orig_edges[i], 0f64), (orig_edges[i + 1], orig_counts[i] as f64)],
                STEEL_BLUE.mix(0.55).filled(),
            )))
            .ok()?;
    }
    for i in 0..uniq_counts.len() {
        if i + 1 >= uniq_edges.len() { break; }
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(uniq_edges[i], 0f64), (uniq_edges[i + 1], uniq_counts[i] as f64)],
                FIREBRICK.mix(0.55).filled(),
            )))
            .ok()?;
    }

    root.present().ok()?;
    Some(out.clone())
}

pub fn plot_u50_summary_bars(outdir: &Path) -> Option<PathBuf> {
    let summ = read_summary_map(&outdir.join("u50_summary.csv"))?;
    let metric_order = ["N50", "NG50", "U50", "UG50"];
    let lengths: Vec<f64> = metric_order
        .iter()
        .map(|m| summ.get(*m).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0))
        .collect();

    if lengths.iter().all(|&v| v == 0.0) {
        return None;
    }

    let ref_len: f64 = summ
        .get("ref_length")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.0);
    let ug50_pct: f64 = summ
        .get("UG50_percent")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.0);

    let y_max = lengths.iter().cloned().fold(ref_len, f64::max) * 1.15;
    let colours = [NAVY, GREEN, FIREBRICK, PURPLE];

    let out = save_path(outdir, "plot_u50_summary");
    let root = BitMapBackend::new(&out, (800, H_STD)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("Assembly metrics (UG50% = {ug50_pct:.2}%)"),
            ("sans-serif", 22),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(0f64..5f64, 0f64..y_max)
        .ok()?;

    chart
        .configure_mesh()
        .y_desc("length (bp)")
        .disable_x_mesh()
        .draw()
        .ok()?;

    for (i, (&label, &val)) in metric_order.iter().zip(lengths.iter()).enumerate() {
        let x = i as f64 + 0.1;
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(x, 0f64), (x + 0.8, val)],
                colours[i].mix(0.85).filled(),
            )))
            .ok()?;
        // value label
        chart
            .draw_series(std::iter::once(Text::new(
                format!("{}", val as u64),
                (x + 0.4, val),
                ("sans-serif", 13).into_font(),
            )))
            .ok()?;
        // x-axis label (metric name)
        chart
            .draw_series(std::iter::once(Text::new(
                label.to_string(),
                (x + 0.3, -y_max * 0.04),
                ("sans-serif", 14).into_font(),
            )))
            .ok()?;
    }

    if ref_len > 0.0 {
        chart
            .draw_series(LineSeries::new(
                [(0f64, ref_len), (5f64, ref_len)],
                Into::<ShapeStyle>::into(RGBColor(130, 130, 130))
                    .stroke_width(1),
            ))
            .ok()?;
    }

    root.present().ok()?;
    Some(out.clone())
}

pub fn plot_u50_coverage(outdir: &Path, bins: usize) -> Option<PathBuf> {
    let summ = read_summary_map(&outdir.join("u50_summary.csv"))?;
    let gi_path = outdir.join("u50_gap_intervals.csv");
    if !gi_path.exists() {
        return None;
    }

    let ref_len: usize = summ.get("ref_length")?.parse().ok().filter(|&v| v > 0)?;

    let mut mask = vec![1u8; ref_len];
    let gaps: Vec<GapIntervalRow> = read_csv(&gi_path).ok()?;
    for g in &gaps {
        let end = g.end.min(ref_len);
        if g.start < end {
            for b in &mut mask[g.start..end] {
                *b = 0;
            }
        }
    }

    let bin_size = (ref_len / bins).max(1);
    let n_bins = ref_len / bin_size;
    if n_bins == 0 {
        return None;
    }

    let cov: Vec<f64> = (0..n_bins)
        .map(|i| {
            let sl = &mask[i * bin_size..(i + 1) * bin_size];
            sl.iter().map(|&b| b as f64).sum::<f64>() / sl.len() as f64
        })
        .collect();

    let out = save_path(outdir, "plot_u50_coverage");
    let root = BitMapBackend::new(&out, (W_WIDE, H_SHORT)).into_drawing_area();
    root.fill(&WHITE).ok()?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Reference coverage by mapped contigs (gaps carved)",
            ("sans-serif", 22),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0f64..(n_bins as f64), 0f64..1.05)
        .ok()?;

    chart
        .configure_mesh()
        .x_desc(format!("reference position (bins of {bin_size} bp)"))
        .y_desc("fraction covered")
		.max_light_lines(0)
        .draw()
        .ok()?;

    chart
        .draw_series(AreaSeries::new(
            cov.iter().enumerate().map(|(i, &v)| (i as f64, v)),
            0f64,
            STEEL_BLUE.mix(0.8),
        ))
        .ok()?;

    root.present().ok()?;
    Some(out.clone())
}

// ══════════════════════════════════════════════════════════════════════════════
// HTML report
// ══════════════════════════════════════════════════════════════════════════════

const HTML_HEAD: &str = r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>rustyomestats report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }
  h1 { border-bottom: 2px solid #444; padding-bottom: .3em; }
  h2 { margin-top: 2em; color: #2a4b7c; }
  .meta { color: #666; font-size: .9em; }
  figure { margin: 1em 0 2em 0; text-align: center; }
  figure img { max-width: 100%; height: auto;
               border: 1px solid #ddd; border-radius: 4px; }
  figcaption { font-size: .9em; color: #555; margin-top: .4em; }
  table { border-collapse: collapse; margin: 1em 0; }
  th, td { border: 1px solid #ccc; padding: .4em .8em; text-align: left; }
  th { background: #f5f5f5; }
  code { background: #f5f5f5; padding: .1em .3em; border-radius: 3px; }
</style>
</head>
<body>
"#;

const HTML_TAIL: &str = "</body>\n</html>\n";

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn embed_png(png_path: &Path) -> Result<String> {
    let bytes = fs::read(png_path)
        .with_context(|| format!("reading {}", png_path.display()))?;
    let b64 = B64.encode(&bytes);
    let alt = escape_html(png_path.file_stem().unwrap_or_default().to_str().unwrap_or(""));
    Ok(format!(
        r#"<img src="data:image/png;base64,{b64}" alt="{alt}">"#
    ))
}

fn render_summary_table(csv_path: &Path, title: &str) -> String {
    #[derive(Deserialize)]
    struct Row {
        metric: String,
        value: String,
    }
    let rows: Vec<Row> = match read_csv(csv_path) {
        Ok(r) => r,
        Err(_) => return String::new(),
    };
    let mut html = format!("<h2>{}</h2>\n", escape_html(title));
    html.push_str("<table><thead><tr><th>metric</th><th>value</th></tr></thead><tbody>\n");
    for r in &rows {
        html.push_str(&format!(
            "<tr><td><code>{}</code></td><td>{}</td></tr>\n",
            escape_html(&r.metric),
            escape_html(&r.value),
        ));
    }
    html.push_str("</tbody></table>\n");
    html
}

fn write_html_report(
    outdir: &Path,
    figures: &[(&str, Option<PathBuf>, &str)],
    when: &str,
) -> Result<PathBuf> {
    let mut out = String::new();

    out.push_str(HTML_HEAD);
    out.push_str("<h1>rustyomestats report</h1>\n");
    out.push_str(&format!(
        "<p class=\"meta\">Generated {} from <code>{}</code>.</p>\n",
        escape_html(when),
        escape_html(&outdir.canonicalize().unwrap_or_else(|_| outdir.to_path_buf()).to_string_lossy()),
    ));

    out.push_str(&render_summary_table(
        &outdir.join("summary_stats.csv"),
        "Genome summary stats (summary_stats.csv)",
    ));
    out.push_str(&render_summary_table(
        &outdir.join("u50_summary.csv"),
        "U50 assembly metrics (u50_summary.csv)",
    ));

    for (title, maybe_png, caption) in figures {
        let png = match maybe_png.as_ref().filter(|p| p.exists()) {
            Some(p) => p,
            None => continue,
        };
        let img_tag = embed_png(png)?;
        out.push_str(&format!(
            "<h2>{}</h2>\n<figure>{img_tag}<figcaption>{}</figcaption></figure>\n",
            escape_html(title),
            escape_html(caption),
        ));
    }

    out.push_str(HTML_TAIL);

    let html_path = outdir.join("report.html");
    fs::write(&html_path, &out).with_context(|| format!("writing {}", html_path.display()))?;
    Ok(html_path)
}

// ══════════════════════════════════════════════════════════════════════════════
// CLI
// ══════════════════════════════════════════════════════════════════════════════

#[derive(Parser, Debug)]
#[command(
    name = "rustyomestats_plots",
    about = "Plot rustyomestats CSV outputs; emit PNGs and a self-contained HTML report"
)]
struct Args {
    /// rustyomestats output directory
    #[arg(short, long, default_value = "rustyomestats_out")]
    dir: PathBuf,

    /// Skip report.html generation (PNGs still written)
    #[arg(long, default_value_t = false)]
    no_html: bool,
}

// ══════════════════════════════════════════════════════════════════════════════
// main
// ══════════════════════════════════════════════════════════════════════════════

pub fn create_plots(
	dir : &std::path::Path,
	html_report : bool,
) -> Result<()> {
    //let args = Args::parse();
    //if !args.dir.is_dir() {
    //    return Err(anyhow!("error: {} is not a directory", args.dir.display()));
    //}

    let d = dir;

    let figures: Vec<(&str, Option<PathBuf>, &str)> = vec![
        (
            "Length histogram",
            plot_length_hist(d),
            "Distribution of sequence lengths (log10 x-axis).",
        ),
        (
            "GC distribution",
            plot_gc_distribution(d),
            "Per-sequence GC% (ambiguous bases excluded from the denominator).",
        ),
        (
            "GC vs length",
            plot_gc_vs_length(d),
            "Per-sequence GC% plotted against log-length.",
        ),
        (
            "Codon usage (6-frame vs FragGeneScan)",
            plot_codon_usage_bar(d),
            "64-codon density. Compare the 6-frame background (absolute) to predicted ORFs (FGS).",
        ),
        (
            "6-frame codon heatmap",
            plot_codon_heatmap_by_frame(d),
            "Codon density across all 6 reading frames.",
        ),
        (
            "Codon enrichment (predicted / absolute)",
            plot_codon_enrichment(d),
            "Codons overrepresented in predicted ORFs relative to the 6-frame background.",
        ),
        (
            "U50: contig length distributions",
            plot_u50_contig_lengths(d),
            "Original contig lengths vs unique-bp lengths after greedy masking.",
        ),
        (
            "U50: assembly metric summary",
            plot_u50_summary_bars(d),
            "N50 / NG50 / U50 / UG50. Dotted line = reference length.",
        ),
        (
            "U50: reference coverage",
            plot_u50_coverage(d, 500),
            "Binned fraction of the reference covered by at least one contig.",
        ),
    ];

    let produced: Vec<&PathBuf> = figures
        .iter()
        .filter_map(|(_, p, _)| p.as_ref())
        .collect();

    println!(
        "[rustyomestats] wrote {} PNG(s) to {}",
        produced.len(),
        d.display()
    );
    for p in &produced {
        println!("                {}", p.file_name().unwrap_or_default().to_string_lossy());
    }

    if html_report {
        let when = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
        let html_path = write_html_report(d, &figures, &when)?;
        println!(
            "[rustyomestats] wrote {}",
            html_path.file_name().unwrap_or_default().to_string_lossy()
        );
    }

    Ok(())
}
