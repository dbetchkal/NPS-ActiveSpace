#!/usr/bin/env bash
# Stage a minimal NMSim runtime into vendor/nmsim-runtime/ for Docker+Wine runs.
#
# Usage:
#   docker/stage_nmsim_runtime.sh /path/to/NMSim-install
#   NMSIM_SOURCE=/Volumes/NPS_ADSB_Data/.../binaries/NMSim docker/stage_nmsim_runtime.sh
#
# Copies: Nord2000batch.exe, *.dll, RND/ (required). Does not copy example cases.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$(pwd)"
DEST="${NMSIM_RUNTIME:-$REPO/vendor/nmsim-runtime}"
SRC="${1:-${NMSIM_SOURCE:-}}"

if [[ -z "$SRC" ]]; then
  echo "Usage: $0 /path/to/NMSim-install" >&2
  echo "  or set NMSIM_SOURCE=/path/to/NMSim-install" >&2
  exit 1
fi
if [[ ! -e "$SRC/Nord2000batch.exe" ]]; then
  echo "ERROR: Nord2000batch.exe not found under $SRC" >&2
  exit 1
fi
if [[ ! -e "$SRC/RND/directories.ini" ]]; then
  echo "ERROR: RND/directories.ini not found under $SRC (required at runtime)" >&2
  exit 1
fi

mkdir -p "$DEST"
echo "[stage] $SRC -> $DEST"
cp -f "$SRC/Nord2000batch.exe" "$DEST/"
cp -f "$SRC"/*.dll "$DEST/" 2>/dev/null || true
rm -rf "$DEST/RND"
cp -R "$SRC/RND" "$DEST/RND"
echo "[stage] done ($(du -sh "$DEST" | cut -f1))"
