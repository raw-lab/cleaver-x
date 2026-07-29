#!/usr/bin/env bash
# Cleaver installer — builds the release binary and drops it on your PATH.
#
#   ./install.sh                 # default backend (rayon)
#   ./install.sh hydra           # HydraMPP backend (multi-core / cross-node)
#   ./install.sh hydra,gpu       # HydraMPP + CUDA stats kernel (needs CUDA)
#   DEST=/usr/local/bin ./install.sh    # choose an explicit install dir
#
# With no DEST, the binary goes to the first of ~/.cargo/bin or ~/.local/bin
# that is already on PATH; otherwise ~/.local/bin (and we offer to add it).
set -euo pipefail

FEATURES="${1:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

command -v cargo >/dev/null 2>&1 || {
    echo "error: cargo not found. Install Rust from https://rustup.rs (need rustc >= 1.75)." >&2
    exit 1
}

echo ">> building cleaver (release${FEATURES:+, features: $FEATURES}) ..."
if [ -n "$FEATURES" ]; then
    cargo build --release --locked -p cleaver-cli --features "$FEATURES"
else
    cargo build --release --locked -p cleaver-cli
fi
BIN="$HERE/target/release/cleaver"
[ -x "$BIN" ] || { echo "error: build did not produce $BIN" >&2; exit 1; }

on_path() { case ":$PATH:" in *":$1:"*) return 0 ;; *) return 1 ;; esac; }

if [ -n "${DEST:-}" ]; then
    INSTALL_DIR="$DEST"
elif on_path "$HOME/.cargo/bin"; then
    INSTALL_DIR="$HOME/.cargo/bin"
elif on_path "$HOME/.local/bin"; then
    INSTALL_DIR="$HOME/.local/bin"
else
    INSTALL_DIR="$HOME/.local/bin"
fi

mkdir -p "$INSTALL_DIR"
install -m 0755 "$BIN" "$INSTALL_DIR/cleaver"
echo ">> installed: $INSTALL_DIR/cleaver"

if ! on_path "$INSTALL_DIR"; then
    SHELL_RC="$HOME/.bashrc"
    [ -n "${ZSH_VERSION:-}" ] && SHELL_RC="$HOME/.zshrc"
    echo ">> $INSTALL_DIR is not on PATH; appending an export to $SHELL_RC"
    printf '\n# added by cleaver install.sh\nexport PATH="%s:$PATH"\n' "$INSTALL_DIR" >> "$SHELL_RC"
    echo ">> open a new shell or: export PATH=\"$INSTALL_DIR:\$PATH\""
fi

echo ">> verifying ..."
"$INSTALL_DIR/cleaver" version || true
echo ">> done. Try: cleaver doctor"
