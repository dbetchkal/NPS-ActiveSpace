#!/usr/bin/env bash
# Build the NPS-ActiveSpace Linux/Wine image. Build context = repo root so the
# Dockerfile can COPY requirements.linux.txt and docker/nord2000.
#   docker/build.sh
set -euo pipefail
cd "$(dirname "$0")/.."
docker build --platform=linux/amd64 -f docker/Dockerfile -t nps-activespace:linux .
echo "[build] done -> nps-activespace:linux"
