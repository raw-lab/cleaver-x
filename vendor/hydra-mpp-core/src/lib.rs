//! rustyomestats — genome statistics, codon density, and U50-family
//! assembly metrics.
//!
//! Modules:
//! * [`io_utils`] — FASTA discovery / loading
//! * [`stats`]    — length / GC / N-L stats (polars output)
//! * [`codon`]    — 6-frame translation + absolute/predicted codon density
//! * [`fgs`]      — FragGeneScanRs subprocess wrapper for predicted ORFs
//! * [`u50`]      — Castro et al. (2016) N50/U50 assembly metrics from a
//!                  reference FASTA + mapped-contigs BED
//!
//! Vendored build note: the pure length/GC/N-L statistics in [`stats`]
//! (`compute_nl`, `BasicStats`) compile with no external dependencies. The
//! heavier functionality (everything below, plus the bio/polars/rayon paths in
//! [`stats`]) is gated behind the `full` feature.

#[cfg(feature = "full")]
pub mod codon;
#[cfg(feature = "full")]
pub mod fgs;
#[cfg(feature = "full")]
pub mod io_utils;
pub mod stats;
#[cfg(feature = "full")]
pub mod u50;
#[cfg(feature = "full")]
pub mod visualization;
