#!/usr/bin/env bash
# Same 3D job as docs/tmp_windows_aam_3d.bat (Docker+Wine on Mac).
set -euo pipefail
cd "$(dirname "$0")/.."

GEN3D=nps_active_space/scripts/generate_3d_active_space.py
SITE=(-e container -u DENA -s TRLA -y 2025 -a nvspl)
ALT=(--min-altitude 900 --max-altitude 2100)
DENSITY=48
EXTRA=(
  --headings 0
  --omni-min 0 --omni-max 2
  --density "$DENSITY"
  --cleanup
  --annotation-file DENATRLA2025_saved_annotations.geojson
)

run_aam=1
run_nmsim=1
case "${1:-both}" in
  aam) run_nmsim=0 ;;
  nmsim) run_aam=0 ;;
  both|"") ;;
  *) echo "Usage: $0 [aam|nmsim|both]" >&2; exit 2 ;;
esac

echo "DENATRLA 3D: layers 900-2100 m, density ${DENSITY}, omni 0-2"

if [[ "$run_aam" -eq 1 ]]; then
  echo "=== AAM 3D (generate + batch + fit) ==="
  docker/run_activespace.sh -m aam "$GEN3D" --model aam "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
fi

if [[ "$run_nmsim" -eq 1 ]]; then
  echo "=== NMSim 3D (generate + batch + fit) ==="
  docker/run_activespace.sh "$GEN3D" --model nmsim "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
fi

echo "Done. fits.csv -> example_data/site_projects/fits.csv"
