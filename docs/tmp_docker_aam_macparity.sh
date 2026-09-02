#!/usr/bin/env bash
# Same 1500 m smoke as docs/tmp_windows_aam_macparity.bat, via Docker+Wine.
set -euo pipefail
cd "$(dirname "$0")/.."

OUTDIR=docker
OUT="${OUTDIR}/denatrla_1500m_macparity_metrics.json"
JOB_FILE=docs/tmp_platform_aam_macparity_job.json
GEN=nps_active_space/scripts/generate_active_space.py
BASE=(-e container -u DENA -s TRLA -y 2025 -l 1500 --headings 0 --omni-min 0 --omni-max 0 --density 10 --cleanup --annotation-file DENATRLA2025_saved_annotations.geojson)

echo "Mac-parity 1500 m smoke (density 10, omni +000) — AAM then NMSim"
rm -f "$OUT" \
  "${OUTDIR}/denatrla_1500m_macparity_aam.json" \
  "${OUTDIR}/denatrla_1500m_macparity_nmsim.json"

docker/run_activespace.sh -m aam docs/tmp_platform_run_metrics.py \
  --out "/repo/${OUT}" --label aam \
  --job-file "/repo/${JOB_FILE}" \
  --launcher docker+wine \
  -- "$GEN" --model aam "${BASE[@]}" \
  --results-out "/repo/${OUTDIR}/denatrla_1500m_macparity_aam.json"

docker/run_activespace.sh docs/tmp_platform_run_metrics.py \
  --out "/repo/${OUT}" --label nmsim \
  --job-file "/repo/${JOB_FILE}" \
  --launcher docker+wine \
  -- "$GEN" "${BASE[@]}" \
  --results-out "/repo/${OUTDIR}/denatrla_1500m_macparity_nmsim.json"
