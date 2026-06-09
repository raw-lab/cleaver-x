//! rustyomestats — genome statistics, codon density, and U50-family
//! assembly metrics.
//!
//! The pure stat kernels ([`stats::compute_nl`], the codon counters in
//! [`codon`], and the N-stat helpers in [`u50`]) build with no external deps on
//! rustc 1.75. The file-I/O and CSV-output layers (polars/bio) and the
//! standalone CLI are behind the `full` feature (rustc >= 1.85).

pub mod codon;
pub mod stats;

#[cfg(feature = "full")]
pub mod fgs;
#[cfg(feature = "full")]
pub mod io_utils;
#[cfg(feature = "full")]
pub mod u50;
