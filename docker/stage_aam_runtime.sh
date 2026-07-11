#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Stage AAM runtime into vendor/aam-runtime/ for Docker+Wine runs.

Usage:
  docker/stage_aam_runtime.sh [-h] /path/to/AAM-install
  AAM_SOURCE=/path/to/AAM docker/stage_aam_runtime.sh

Environment:
  AAM_SOURCE      source install dir (alternative to positional arg)
  AAM_RUNTIME     destination (default: ./vendor/aam-runtime)
  AAM_EXE         executable name to copy (default: AAM_3.0.0.exe)

Copies AAM_3.0.0.exe, NCfiles/, and noisecon.inp when present in the source dir
(needed for docker/validate_aam_smoke.py). Does not copy AAM.config.
EOF
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
esac

cd "$(dirname "$0")/.."
REPO="$(pwd)"
DEST="${AAM_RUNTIME:-$REPO/vendor/aam-runtime}"
SRC="${1:-${AAM_SOURCE:-}}"
AAM_EXE="${AAM_EXE:-AAM_3.0.0.exe}"

if [[ -z "$SRC" ]]; then
  usage >&2
  exit 1
fi
if [[ ! -e "$SRC/$AAM_EXE" ]]; then
  echo "ERROR: $AAM_EXE not found under $SRC" >&2
  exit 1
fi
if [[ ! -d "$SRC/NCfiles" ]]; then
  echo "ERROR: NCfiles/ not found under $SRC" >&2
  exit 1
fi

mkdir -p "$DEST"
echo "[stage] $SRC -> $DEST"
cp -f "$SRC/$AAM_EXE" "$DEST/"
rm -rf "$DEST/NCfiles"
cp -R "$SRC/NCfiles" "$DEST/NCfiles"

if [[ -f "$SRC/noisecon.inp" ]]; then
  cp -f "$SRC/noisecon.inp" "$DEST/"
  echo "[stage] copied noisecon.inp (AAM smoke test input)"
else
  echo "[stage] WARNING: noisecon.inp not in source — smoke test needs it in vendor/aam-runtime/" >&2
fi

echo "[stage] done ($(du -sh "$DEST" | cut -f1))"
