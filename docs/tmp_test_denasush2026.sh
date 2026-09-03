#!/usr/bin/env bash
# DENASUSH 2026 3D pipeline — edit CONFIG below.
# Usage: bash docs/tmp_test_denasush2026.sh [aam|nmsim|both]
set -euo pipefail
cd "$(dirname "$0")/.."

# --- CONFIG (edit these) ---
ENV=container
UNIT=DENA
SITE_CODE=SUSH
YEAR=2026
MIN_ALTITUDE=1200
MAX_ALTITUDE=2400
DENSITY=10
# Omni ladder step is 0.5 dB (0-2 = gains 0, 0.5, 1, 1.5, 2).
OMNI_MIN=0
OMNI_MAX=2
HEADINGS=0
ANNOT=DENASUSH2026_saved_annotations.geojson
# Set to 1 to record wall/cpu in docker/denasush2026_test_metrics.json
USE_METRICS=0
# --- end CONFIG ---

GEN3D=nps_active_space/scripts/generate_3d_active_space.py
SITE=(-e "$ENV" -u "$UNIT" -s "$SITE_CODE" -y "$YEAR" -a nvspl)
ALT=(--min-altitude "$MIN_ALTITUDE" --max-altitude "$MAX_ALTITUDE")
EXTRA=(
  --headings "$HEADINGS"
  --omni-min "$OMNI_MIN" --omni-max "$OMNI_MAX"
  --density "$DENSITY"
  --cleanup
  --annotation-file "$ANNOT"
)

OUTDIR=docker
OUT="${OUTDIR}/denasush2026_test_metrics.json"
JOB_FILE=docs/tmp_platform_test_denasush2026_job.json

RUN_AAM=1
RUN_NMSIM=1
case "${1:-both}" in
  aam) RUN_NMSIM=0 ;;
  nmsim) RUN_AAM=0 ;;
  both|"") ;;
  *)
    echo "Usage: $0 [aam|nmsim|both]" >&2
    exit 1
    ;;
esac

echo "DENASUSH2026 3D: layers ${MIN_ALTITUDE}-${MAX_ALTITUDE} m (300 m step), density=${DENSITY}, omni ${OMNI_MIN}-${OMNI_MAX}"
echo "Models: AAM=${RUN_AAM} NMSim=${RUN_NMSIM}"
echo

if [[ "$USE_METRICS" == 1 ]]; then
  rm -f "$OUT"
fi

run_3d() {
  local model="$1"
  local label="$2"
  if [[ "$USE_METRICS" == 1 ]]; then
    docker/run_activespace.sh -m "$model" docs/tmp_platform_run_metrics.py \
      --out "/repo/${OUT}" --label "$label" \
      --job-file "/repo/${JOB_FILE}" \
      --launcher docker+wine \
      -- "$GEN3D" --model "$model" "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
  else
    if [[ "$model" == "aam" ]]; then
      docker/run_activespace.sh -m aam "$GEN3D" --model aam "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
    else
      docker/run_activespace.sh "$GEN3D" --model nmsim "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
    fi
  fi
}

if [[ "$RUN_AAM" == 1 ]]; then
  echo "=== AAM 3D (generate + batch + fit) ==="
  run_3d aam aam_3d
fi

if [[ "$RUN_NMSIM" == 1 ]]; then
  echo "=== NMSim 3D (generate + batch + fit) ==="
  run_3d nmsim nmsim_3d
fi

echo "Done. fits.csv under example_data/site_projects/ (or your [project] dir)"
if [[ "$USE_METRICS" == 1 ]]; then
  echo "Metrics: ${OUT}"
fi
