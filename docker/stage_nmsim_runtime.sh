#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Stage a minimal NMSim runtime into vendor/nmsim-runtime/ for Docker+Wine runs.

Usage:
  docker/stage_nmsim_runtime.sh [-h] /path/to/NMSim-install
  NMSIM_SOURCE=/path/to/NMSim docker/stage_nmsim_runtime.sh

Environment:
  NMSIM_SOURCE    source install dir (alternative to positional arg)
  NMSIM_RUNTIME   destination (default: ./vendor/nmsim-runtime)

Copies Nord2000batch.exe, *.dll, and RND/ (not example cases).
EOF
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
esac

cd "$(dirname "$0")/.."
REPO="$(pwd)"
DEST="${NMSIM_RUNTIME:-$REPO/vendor/nmsim-runtime}"
SRC="${1:-${NMSIM_SOURCE:-}}"

if [[ -z "$SRC" ]]; then
  usage >&2
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
