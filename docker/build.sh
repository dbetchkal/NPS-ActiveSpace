#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build the NPS-ActiveSpace Linux/Wine image (nps-activespace:linux).

Usage:
  docker/build.sh [-h]

Build context is the repo root (pyproject.toml, nps_active_space/, docker/nord2000).
First build is slow (~13 min on Apple Silicon under Rosetta).
EOF
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
esac

cd "$(dirname "$0")/.."
AAM_SRC="${AAM_TRANSLATOR:-$(dirname "$(pwd)")/aam_translator}"
STAGE="vendor/.docker-build/aam_translator"
if [[ ! -d "$AAM_SRC" ]]; then
  echo "ERROR: aam_translator not found at $AAM_SRC (set AAM_TRANSLATOR=...)" >&2
  exit 1
fi
rm -rf vendor/.docker-build
mkdir -p vendor/.docker-build
cp -R "$AAM_SRC" "$STAGE"
docker build --platform=linux/amd64 -f docker/Dockerfile -t nps-activespace:linux .
echo "[build] done -> nps-activespace:linux"
