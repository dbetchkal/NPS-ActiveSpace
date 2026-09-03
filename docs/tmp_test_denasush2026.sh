#!/usr/bin/env bash
# DENASUSH 2026 smoke — same flags as tmp_windows_aam_macparity.bat, site DENA/SUSH/2026.
set -euo pipefail
cd "$(dirname "$0")/.."

OUTDIR=docker
OUT="${OUTDIR}/denasush2026_test_metrics.json"
JOB_FILE=docs/tmp_platform_test_denasush2026_job.json
GEN=nps_active_space/scripts/generate_active_space.py
BASE=(-e container -u DENA -s SUSH -y 2026 -l 1500 --headings 0 --omni-min 0 --omni-max 0 --density 10 --cleanup --annotation-file DENASUSH2026_saved_annotations.geojson)

echo "DENASUSH2026 1500 m smoke (density 10, omni +000) — AAM then NMSim"
rm -f "$OUT" \
  "${OUTDIR}/denasush2026_test_aam.json" \
  "${OUTDIR}/denasush2026_test_nmsim.json"

docker/run_activespace.sh -m aam docs/tmp_platform_run_metrics.py \
  --out "/repo/${OUT}" --label aam \
  --job-file "/repo/${JOB_FILE}" \
  --launcher docker+wine \
  -- "$GEN" --model aam "${BASE[@]}" \
  --results-out "/repo/${OUTDIR}/denasush2026_test_aam.json"

docker/run_activespace.sh docs/tmp_platform_run_metrics.py \
  --out "/repo/${OUT}" --label nmsim \
  --job-file "/repo/${JOB_FILE}" \
  --launcher docker+wine \
  -- "$GEN" "${BASE[@]}" \
  --results-out "/repo/${OUTDIR}/denasush2026_test_nmsim.json"
