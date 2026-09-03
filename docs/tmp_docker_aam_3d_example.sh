#!/usr/bin/env bash
# Mac Docker twin of docs/tmp_windows_aam_3d_example.bat.
set -euo pipefail
cd "$(dirname "$0")/.."

GEN3D=nps_active_space/scripts/generate_3d_active_space.py
SITE=(-e container -u DENA -s TRLA -y 2025 -a nvspl)
ALT=(--min-altitude 1200 --max-altitude 1800)
ANNOT=DENATRLA2025_saved_annotations.geojson
EXTRA=(
  --headings 0
  --omni-min 0 --omni-max 0
  --density 10
  --cleanup
  --annotation-file "$ANNOT"
)

run_aam=1
run_nmsim=1
case "${1:-both}" in
  aam) run_nmsim=0 ;;
  nmsim) run_aam=0 ;;
  both|"") ;;
  *) echo "Usage: $0 [aam|nmsim|both]" >&2; exit 2 ;;
esac

annot_path="example_data/site_projects/DENATRLA/$ANNOT"
if [[ ! -f "$annot_path" ]]; then
  echo "ERROR: example annotations not found: $annot_path" >&2
  exit 1
fi

echo "DENATRLA 3D example smoke: layers 1200-1800 m, density 10, omni +000"
echo "Annotations: $annot_path"

if [[ "$run_aam" -eq 1 ]]; then
  echo "=== AAM 3D example (generate + batch + fit) ==="
  docker/run_activespace.sh -m aam "$GEN3D" --model aam "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
fi

if [[ "$run_nmsim" -eq 1 ]]; then
  echo "=== NMSim 3D example (generate + batch + fit) ==="
  docker/run_activespace.sh "$GEN3D" --model nmsim "${SITE[@]}" "${ALT[@]}" "${EXTRA[@]}"
fi

echo "Done. fits.csv -> example_data/site_projects/fits.csv"
