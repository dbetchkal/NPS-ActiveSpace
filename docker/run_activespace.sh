#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run a script inside the nps-activespace:linux container (Python + Wine + acoustic model).

Usage:
  docker/run_activespace.sh [-h] [-m nmsim|aam] <script.py> [args...]
  docker/run_activespace.sh [-h] [-m nmsim|aam] bash

Options:
  -m, --model     acoustic model runtime to mount: nmsim (default) or aam

Mounts (depends on model):
  repo              -> /repo       (read-write)
  NMSim runtime     -> /opt/nmsim  (read-only, when model is nmsim)
  AAM runtime       -> /opt/aam    (read-only, when model is aam)
  DATA_DRIVE        -> /data       (optional, read-only)

Environment:
  ACOUSTIC_MODEL    same as -m (env wins if set before parsing; -m is clearer on CLI)
  NMSIM_RUNTIME     path to Nord2000batch.exe + RND/ (default: ./vendor/nmsim-runtime)
  AAM_RUNTIME       path to AAM exe + NCfiles/ (default: ./vendor/aam-runtime)
  DATA_DRIVE        optional host data directory

Examples:
  docker/run_activespace.sh docker/validate_active_space.py -u DENA -s TRLA -y 2025
  docker/run_activespace.sh -m aam docker/validate_aam_smoke.py
  docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py -e container ...
EOF
}

_validate_nmsim_runtime() {
  local rt="$1"
  [[ -e "$rt/Nord2000batch.exe" && -e "$rt/RND/directories.ini" ]]
}

_validate_aam_runtime() {
  local rt="$1"
  local exe="${AAM_EXE:-AAM_3.0.0.exe}"
  [[ -e "$rt/$exe" && -d "$rt/NCfiles" ]]
}

ACOUSTIC_MODEL="${ACOUSTIC_MODEL:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    -m|--model)
      [[ $# -ge 2 ]] || { echo "[run] ERROR: -m requires nmsim or aam" >&2; exit 1; }
      ACOUSTIC_MODEL="$2"
      shift 2
      ;;
    --) shift; break ;;
    -*) echo "[run] ERROR: unknown option: $1" >&2; usage >&2; exit 1 ;;
    *) break ;;
  esac
done

cd "$(dirname "$0")/.."
REPO="$(pwd)"
ACOUSTIC_MODEL="${ACOUSTIC_MODEL:-nmsim}"

case "$ACOUSTIC_MODEL" in
  nmsim|aam) ;;
  *)
    echo "[run] ERROR: model must be nmsim or aam (got: $ACOUSTIC_MODEL)" >&2
    exit 1
    ;;
esac

mounts=(-v "$REPO:/repo")
docker_env=(
  -e PYTHONUNBUFFERED=1
  -e MPLBACKEND=Agg
  -e MPLCONFIGDIR=/tmp/matplotlib
  -e ACOUSTIC_MODEL="$ACOUSTIC_MODEL"
)

if [[ "$ACOUSTIC_MODEL" == "nmsim" ]]; then
  NMSIM_RUNTIME="${NMSIM_RUNTIME:-$REPO/vendor/nmsim-runtime}"
  if ! _validate_nmsim_runtime "$NMSIM_RUNTIME"; then
    echo "[run] ERROR: NMSim runtime not found at $NMSIM_RUNTIME" >&2
    echo "[run] Populate it: docker/stage_nmsim_runtime.sh /path/to/NMSim" >&2
    echo "[run] Or set NMSIM_RUNTIME=/path/to/runtime" >&2
    exit 1
  fi
  mounts+=(-v "$NMSIM_RUNTIME:/opt/nmsim:ro")
  docker_env+=(-e NMSIM_HOME=/opt/nmsim)
fi

if [[ "$ACOUSTIC_MODEL" == "aam" ]]; then
  AAM_RUNTIME="${AAM_RUNTIME:-$REPO/vendor/aam-runtime}"
  if ! _validate_aam_runtime "$AAM_RUNTIME"; then
    echo "[run] ERROR: AAM runtime not found at $AAM_RUNTIME" >&2
    echo "[run] Populate it: docker/stage_aam_runtime.sh /path/to/AAM" >&2
    echo "[run] Or set AAM_RUNTIME=/path/to/runtime" >&2
    exit 1
  fi
  mounts+=(-v "$AAM_RUNTIME:/opt/aam:ro")
  docker_env+=(-e AAM_HOME=/opt/aam)
fi

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
  "${docker_env[@]}"
  "${mounts[@]}"
  -w /repo nps-activespace:linux
  bash -c 'exec "$@"' _
  "$@"
)

{
  printf '[run] docker run --rm --platform=linux/amd64'
  for e in "${docker_env[@]}"; do printf ' %q' "$e"; done
  for m in "${mounts[@]}"; do printf ' %q' "$m"; done
  printf ' -w /repo nps-activespace:linux bash -c %q _' 'exec "$@"'
  for arg in "$@"; do printf ' %q' "$arg"; done
  echo
} >&2
echo "[run] model=$ACOUSTIC_MODEL — launching container..." >&2
exec "${docker_cmd[@]}"
