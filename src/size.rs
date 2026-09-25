//! Human-readable size parsing.
//!
//! A faithful port of the original `human2bytes` helper, with two deliberate
//! improvements:
//!   * intermediate maths use `u128` so the high tiers never overflow/panic;
//!   * single-letter customary units are case-insensitive (so `m`, `g`, `t`
//!     work, and the historical `k`->`K` alias falls out for free).
//!
//! Everything else matches the original contract, including truncation toward
//! zero (`0.1 byte` -> `0`).

use anyhow::{bail, Result};

const CUSTOMARY: [&str; 9] = ["B", "K", "M", "G", "T", "P", "E", "Z", "Y"];
const CUSTOMARY_EXT: [&str; 9] =
    ["byte", "kilo", "mega", "giga", "tera", "peta", "exa", "zetta", "iotta"];
const IEC: [&str; 9] = ["Bi", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi"];
const IEC_EXT: [&str; 9] =
    ["byte", "kibi", "mebi", "gibi", "tebi", "pebi", "exbi", "zebi", "yobi"];

/// Parse a human-readable size string into a number of bytes.
///
/// ```
/// use cleaver::parse_size;
/// assert_eq!(parse_size("0 B").unwrap(), 0);
/// assert_eq!(parse_size("1 K").unwrap(), 1024);
/// assert_eq!(parse_size("1 M").unwrap(), 1_048_576);
/// assert_eq!(parse_size("1 Gi").unwrap(), 1_073_741_824);
/// assert_eq!(parse_size("1 tera").unwrap(), 1_099_511_627_776);
/// assert_eq!(parse_size("0.5kilo").unwrap(), 512);
/// assert_eq!(parse_size("0.1  byte").unwrap(), 0);
/// assert_eq!(parse_size("1 k").unwrap(), 1024); // k is an alias for K
/// assert_eq!(parse_size("100m").unwrap(), 104_857_600); // lowercase customary
/// assert!(parse_size("12 foo").is_err());
/// ```
pub fn parse_size(s: &str) -> Result<u64> {
    let init = s;
    let bytes = s.as_bytes();

    // Leading run of digits and '.' forms the number.
    let mut i = 0;
    while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
        i += 1;
    }
    let num: f64 = s[..i]
        .parse()
        .map_err(|_| anyhow::anyhow!("can't interpret {init:?}"))?;
    let letter = s[i..].trim();

    // Exact match against any symbol set (preserves IEC mixed-case + long forms).
    let sets: [&[&str]; 4] = [&CUSTOMARY, &CUSTOMARY_EXT, &IEC, &IEC_EXT];
    let mut idx: Option<usize> = None;
    for set in sets {
        if let Some(pos) = set.iter().position(|&x| x == letter) {
            idx = Some(pos);
            break;
        }
    }

    // Fallback: case-insensitive single-letter customary (k/m/g/t/...).
    let idx = match idx {
        Some(p) => p,
        None if letter.len() == 1 => CUSTOMARY
            .iter()
            .position(|&x| x == letter.to_ascii_uppercase())
            .ok_or_else(|| anyhow::anyhow!("can't interpret {init:?}"))?,
        None => bail!("can't interpret {init:?}"),
    };

    // factor = 1024^idx, in u128 so high tiers (Z, Y) cannot overflow.
    let factor = 1u128 << (idx * 10);
    let result = num * (factor as f64);
    if !result.is_finite() || result < 0.0 {
        bail!("can't interpret {init:?}");
    }
    if result > u64::MAX as f64 {
        bail!("size too large: {init:?}");
    }
    Ok(result as u64)
}

#[cfg(test)]
mod tests {
    use super::parse_size;

    #[test]
    fn customary_and_iec() {
        assert_eq!(parse_size("0 B").unwrap(), 0);
        assert_eq!(parse_size("1 K").unwrap(), 1024);
        assert_eq!(parse_size("1 M").unwrap(), 1_048_576);
        assert_eq!(parse_size("1 Gi").unwrap(), 1_073_741_824);
        assert_eq!(parse_size("1 tera").unwrap(), 1_099_511_627_776);
        assert_eq!(parse_size("1000M").unwrap(), 1_048_576_000);
    }

    #[test]
    fn fractions_and_aliases() {
        assert_eq!(parse_size("0.5kilo").unwrap(), 512);
        assert_eq!(parse_size("0.1  byte").unwrap(), 0);
        assert_eq!(parse_size("1 k").unwrap(), 1024);
        assert_eq!(parse_size("64m").unwrap(), 67_108_864);
        assert_eq!(parse_size("2g").unwrap(), 2_147_483_648);
    }

    #[test]
    fn rejects_garbage() {
        assert!(parse_size("12 foo").is_err());
        assert!(parse_size("M").is_err());
        assert!(parse_size("1.2.3M").is_err());
        assert!(parse_size("").is_err());
    }
}
