#!/usr/bin/env bash
# Build the WASM target via Emscripten.
#
# Requires emcc / emcmake in PATH. Outputs vikngs-core.{js,wasm} into
# web/public/ where the browser app loads them.
set -euo pipefail

REPO="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$REPO/build-wasm"

if ! command -v emcmake >/dev/null 2>&1; then
    echo "ERROR: emcmake not found in PATH." >&2
    echo "  Ubuntu: sudo apt-get install emscripten" >&2
    echo "  Upstream: https://emscripten.org/docs/getting_started/downloads.html" >&2
    exit 1
fi

mkdir -p "$REPO/web/public"
emcmake cmake -B "$BUILD_DIR" -S "$REPO" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" -j

echo
echo "Build complete. Output:"
ls -la "$REPO/web/public/"vikngs-core.* 2>/dev/null || true
