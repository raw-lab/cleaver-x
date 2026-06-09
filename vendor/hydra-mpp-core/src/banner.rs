//! The HydraMPP startup banner — a figlet-style wordmark in the RAW-Lab house
//! style (cf. SABER / RustyOmeStats). Printed once on `init()` unless quieted.

/// Figlet "Standard" wordmark for HydraMPP.
pub const BANNER: &str = r#"
  _   _           _           __  __ ____  ____
 | | | |_   _  __| |_ __ __ _|  \/  |  _ \|  _ \
 | |_| | | | |/ _` | '__/ _` | |\/| | |_) | |_) |
 |  _  | |_| | (_| | | | (_| | |  | |  __/|  __/
 |_| |_|\__, |\__,_|_|  \__,_|_|  |_|_|   |_|
        |___/    many-headed parallel processing
"#;

/// Render the banner together with the version and the active mode line.
pub fn render(version: &str, mode: &str, cpus: usize, gpus: usize) -> String {
    format!(
        "{BANNER}\n           HydraMPP v{version}  •  mode: {mode}  •  local CPUs: {cpus}  •  local GPUs: {gpus}\n           laptop  →  workstation  →  GPU cluster  •  one binary\n"
    )
}
