#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Confirm Wine can launch the staged acoustic-model binary.

Usage:
  docker/smoke.sh        # Nord2000batch.exe via /usr/local/bin/nord2000
  docker/smoke.sh aam    # AAM_3.0.0.exe via /usr/local/bin/aam (needs noisecon.inp)

Does not run the active-space pipeline. After this passes:
  docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py ...
EOF
}

cd "$(dirname "$0")/.."

case "${1:-nmsim}" in
  -h|--help) usage; exit 0 ;;
  nmsim)
    exec docker/run_activespace.sh docker/validate_nmsim_smoke.py
    ;;
  aam)
    exec docker/run_activespace.sh -m aam docker/validate_aam_smoke.py
    ;;
  *)
    echo "usage: docker/smoke.sh [nmsim|aam]" >&2
    exit 2
    ;;
esac
