#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Confirm the Docker+Wine image can run an acoustic model.

Usage:
  docker/smoke.sh        # NMSim on DENATRLA example data
  docker/smoke.sh aam    # AAM binary (needs staged vendor/aam-runtime)

Requires: docker/build.sh, staged runtime, container.config (see docker/README.md).
EOF
}

cd "$(dirname "$0")/.."

case "${1:-nmsim}" in
  -h|--help) usage; exit 0 ;;
  nmsim)
    exec docker/run_activespace.sh docker/validate_active_space.py
    ;;
  aam)
    exec docker/run_activespace.sh -m aam docker/validate_aam_smoke.py
    ;;
  *)
    echo "usage: docker/smoke.sh [nmsim|aam]" >&2
    exit 2
    ;;
esac
