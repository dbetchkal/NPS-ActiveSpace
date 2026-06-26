#!/usr/bin/env bash
# Run an NPS-ActiveSpace script inside the Linux/Wine container.
#
# Mounts:
#   repo            -> /work        (read-write)
#   NMSim runtime   -> /opt/nmsim   (read-only; see vendor/nmsim-runtime/README.md)
#   data drive      -> /data        (read-only, optional — only if DATA_DRIVE is set)
#
# Usage:
#   docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py -e container ...
#   docker/run_activespace.sh bash
#
# Environment:
#   NMSIM_RUNTIME   path to Nord2000batch.exe + RND/ (default: ./vendor/nmsim-runtime)
#   DATA_DRIVE      optional host path mounted at /data (no default)
set -euo pipefail
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

mounts=(-v "$REPO:/work" -v "$NMSIM_RUNTIME:/opt/nmsim:ro")
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

exec docker run --rm --platform=linux/amd64 "${mounts[@]}" -w /work nps-activespace:linux \
  bash -c 'exec "$@"' _ "$@"
