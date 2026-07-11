# Mac / Linux: Docker + Wine for acoustic models

NMSim (and optionally AAM) are Windows-only. On Mac/Linux, we run inside a container with
Python 3.12 + GDAL and execute the Windows binaries through Wine.

Note that given the additional layers of indirection to run the docker setup, there may be performance slowdowns vs. running natively on Windows.

**Prerequisites:** [Docker Desktop](https://docs.docker.com/get-started/get-docker/) for running the containerized Wine setup; on Apple Silicon enable Rosetta for amd64 emulation.

## One-time setup

```bash
# NMSim runtime (~10 MB; not in git — populate vendor/nmsim-runtime/)
docker/stage_nmsim_runtime.sh /path/to/NMSim-install

# Build the image (~13 min first time; installs deps from pyproject.toml)
docker/build.sh

# Config (if not already present)
cp nps_active_space/config/container_example.config nps_active_space/config/container.config
```

NMSim and AAM binaries are **not redistributable** (NPS internal). Do not commit them to a public repo.

## Run (NMSim — default)

```bash
# Smoke test (DENATRLA example data, no annotations required)
docker/run_activespace.sh docker/validate_active_space.py \
  -u DENA -s TRLA -y 2025 --gains 0 --altitude 1000 --density 10

# Full active-space script (needs annotations for gain selection)
docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000
```

Optional: mount a data drive at `/data` inside the container:

```bash
DATA_DRIVE=/Volumes/NPS_ADSB_Data docker/run_activespace.sh ...
```

Override runtime location: `NMSIM_RUNTIME=/path/to/runtime docker/run_activespace.sh ...`

Use `-e container` with absolute `/repo/...` paths in config. Native viz and ground-truthing on the host use `-e DENA_example` (or your own config) and a local venv — see [example_data/README.md](../example_data/README.md).

Windows setup is unchanged — see root [README.md](../README.md) Installation.

## AAM smoke test (not wired into pipeline)

AAM runtime lives in `vendor/aam-runtime/` (gitignored). Stage from a directory that
includes `AAM_3.0.0.exe`, `NCfiles/`, and `noisecon.inp` (e.g. experiments
`runs/aam_noisecon/` or the vendor install):

```bash
docker/stage_aam_runtime.sh /path/to/AAM_v3_dec2020
# or: docker/stage_aam_runtime.sh ~/dev/nmsim-aam-experiments/runs/aam_noisecon

docker/run_activespace.sh -m aam docker/validate_aam_smoke.py
```

Model selection is a **CLI flag** on `run_activespace.sh` (`-m nmsim|aam`), not a separate
config file — the pipeline still uses `container.config` with `project.nmsim` only.

| | NMSim (default) | AAM |
|-|-------|-----|
| Select | `docker/run_activespace.sh ...` | `docker/run_activespace.sh -m aam ...` |
| Staging | `docker/stage_nmsim_runtime.sh` | `docker/stage_aam_runtime.sh` |
| Local dir | `vendor/nmsim-runtime/` | `vendor/aam-runtime/` |
| Mount | `/opt/nmsim` | `/opt/aam` |
| Shim | `/usr/local/bin/nord2000` | `/usr/local/bin/aam` |
