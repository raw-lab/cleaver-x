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
    cargo build --release --locked --features "$FEATURES"
else
    cargo build --release --locked
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
    # Never edit the user's shell rc files behind their back: print the line and
    # let them add it (set CLEAVER_EDIT_RC=1 to opt in to the old behaviour).
    if [ "${CLEAVER_EDIT_RC:-0}" = "1" ]; then
        SHELL_RC="$HOME/.bashrc"
        [ -n "${ZSH_VERSION:-}" ] && SHELL_RC="$HOME/.zshrc"
        printf '\n# added by cleaver install.sh\nexport PATH="%s:$PATH"\n' "$INSTALL_DIR" >> "$SHELL_RC"
        echo ">> appended a PATH export to $SHELL_RC (CLEAVER_EDIT_RC=1)"
    else
        echo ">> NOTE: $INSTALL_DIR is not on your PATH. Add this line to your shell rc:"
        echo "       export PATH=\"$INSTALL_DIR:\$PATH\""
    fi
fi

# Verify for real: a failure here must fail the install, not be swallowed.
echo ">> verifying ..."
if ! "$INSTALL_DIR/cleaver" version; then
    echo "error: installed binary failed to run: $INSTALL_DIR/cleaver" >&2
    exit 1
fi
if ! "$INSTALL_DIR/cleaver" doctor; then
    echo "error: 'cleaver doctor' reported a failure (see above)" >&2
    exit 1
fi
echo ">> done. cleaver is installed and healthy."
