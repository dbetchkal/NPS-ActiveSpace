#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run a script inside the nps-activespace:linux container (Python + Wine + NMSim).

Usage:
  docker/run_activespace.sh [-h] <script.py> [args...]
  docker/run_activespace.sh [-h] bash

Mounts:
  repo          -> /repo      (read-write)
  NMSim runtime -> /opt/nmsim (read-only)
  DATA_DRIVE    -> /data      (optional, read-only)

Environment:
  NMSIM_RUNTIME   path to Nord2000batch.exe + RND/ (default: ./vendor/nmsim-runtime)
  DATA_DRIVE      optional host data directory

Examples:
  docker/run_activespace.sh docker/validate_active_space.py -u DENA -s TRLA -y 2025
  docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py -e container ...
EOF
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
esac

cd "$(dirname "$0")/.."
REPO="$(pwd)"
NMSIM_RUNTIME="${NMSIM_RUNTIME:-$REPO/vendor/nmsim-runtime}"

if [[ ! -e "$NMSIM_RUNTIME/Nord2000batch.exe" ]]; then
  echo "[run] ERROR: NMSim runtime not found at $NMSIM_RUNTIME" >&2
  echo "[run] Populate it: docker/stage_nmsim_runtime.sh /path/to/NMSim" >&2
  echo "[run] Or set NMSIM_RUNTIME=/path/to/runtime" >&2
  exit 1
fi
if [[ ! -e "$NMSIM_RUNTIME/RND/directories.ini" ]]; then
  echo "[run] ERROR: $NMSIM_RUNTIME/RND/directories.ini missing (incomplete runtime)" >&2
  exit 1
fi

mounts=(-v "$REPO:/repo" -v "$NMSIM_RUNTIME:/opt/nmsim:ro")
if [[ -n "${DATA_DRIVE:-}" ]]; then
  if [[ -d "$DATA_DRIVE" ]]; then
    mounts+=(-v "$DATA_DRIVE:/data:ro")
  else
    echo "[run] WARNING: DATA_DRIVE=$DATA_DRIVE not found; continuing without /data" >&2
  fi
fi

if [[ $# -eq 0 ]]; then set -- bash; fi
if [[ "$1" == *.py ]]; then
  set -- /opt/venv/bin/python -u -W ignore "$@"
fi

docker_cmd=(
  docker run --rm --platform=linux/amd64
  -e PYTHONUNBUFFERED=1
  -e MPLBACKEND=Agg
  -e MPLCONFIGDIR=/tmp/matplotlib
  "${mounts[@]}"
  -w /repo nps-activespace:linux
  bash -c 'exec "$@"' _
  "$@"
)

{
  printf '[run] docker run --rm --platform=linux/amd64 -e PYTHONUNBUFFERED=1 -e MPLBACKEND=Agg -e MPLCONFIGDIR=/tmp/matplotlib'
  for m in "${mounts[@]}"; do printf ' %q' "$m"; done
  printf ' -w /repo nps-activespace:linux bash -c %q _' 'exec "$@"'
  for arg in "$@"; do printf ' %q' "$arg"; done
  echo
} >&2
echo '[run] launching container...' >&2
exec "${docker_cmd[@]}"
