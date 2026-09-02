# Mac / Linux: Docker + Wine for acoustic models

NMSim (and optionally AAM) are Windows-only. On Mac/Linux, we run inside a container with
Python 3.12 + GDAL and execute the Windows binaries through Wine.

Note that given the additional layers of indirection to run the docker setup, there may be performance slowdowns vs. running natively on Windows.

**Prerequisites:** [Docker Desktop](https://docs.docker.com/get-started/get-docker/) for running the containerized Wine setup; on Apple Silicon enable Rosetta for amd64 emulation.

## One-time setup

```bash
# NMSim runtime (~10 MB; not in git — populate vendor/nmsim-runtime/)
docker/stage_nmsim_runtime.sh /path/to/NMSim-install

# Build the image (~13 min first time; pip installs aam-translator from GitHub per pyproject.toml)
docker/build.sh

# Config (if not already present)
cp nps_active_space/config/container_example.config nps_active_space/config/container.config
```

NMSim and AAM binaries are **not redistributable** (NPS internal). Do not commit them to a public repo.

## Run (NMSim — default)

```bash
# Primary smoke test — one gain, coarse mesh (DENATRLA example data)
docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000 \
  --omni-min 0 --omni-max 0 --density 10 --headings 0

# Full active-space script (needs annotations for gain selection)
docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000
```

### AAM via `generate_active_space` (production path)

AAM runs need **both** layers of model selection:

- `docker/run_activespace.sh -m aam` — mounts the AAM Wine runtime at `/opt/aam` and caps Wine omni workers (`AAM_PARALLEL_N=2`)
- `--model aam` on the Python script — selects `AamPropagationModel` and `Output_Data/aam/` layout

```bash
# Primary AAM smoke (same coarse mesh as NMSim smoke above)
docker/run_activespace.sh -m aam nps_active_space/scripts/generate_active_space.py \
  -e container --model aam -u DENA -s TRLA -y 2025 -l 1000 \
  --omni-min 0 --omni-max 0 --density 10 --headings 0

# Full AAM generation (gain sweep; writes Output_Data/aam/; run fit_3d for project fits.csv)
docker/run_activespace.sh -m aam nps_active_space/scripts/generate_active_space.py \
  -e container --model aam -u DENA -s TRLA -y 2025 -l 1000
```

Optional: mount a data drive at `/data` inside the container:

```bash
DATA_DRIVE=/Volumes/NPS_ADSB_Data docker/run_activespace.sh ...
```

Override runtime location: `NMSIM_RUNTIME=/path/to/runtime docker/run_activespace.sh ...`

Use `-e container` with absolute `/repo/...` paths in config. Native viz and ground-truthing on the host use `-e DENA_example` (or your own config) and a local venv — see [example_data/README.md](../example_data/README.md).

Windows setup is unchanged — see root [README.md](../README.md) Installation.

## AAM smoke test

AAM runtime lives in `vendor/aam-runtime/` (gitignored). Stage from a directory that
includes `AAM_3.0.0.exe`, `NCfiles/`, and `noisecon.inp` (e.g. experiments
`runs/aam_noisecon/` or the vendor install):

```bash
docker/stage_aam_runtime.sh /path/to/AAM_v3_dec2020
# or: docker/stage_aam_runtime.sh ~/dev/nmsim-aam-experiments/runs/aam_noisecon

docker/run_activespace.sh -m aam docker/validate_aam_smoke.py

# AAM propagation adapter (two-point ridge reciprocal run)
docker/stage_aam_runtime.sh tests/active_space/fixtures/two_point_ridge
docker/run_activespace.sh -m aam docker/validate_aam_propagation_model.py

# DENATRLA ELV pre-filter vs AAM below-ground (density 10, 1500 m)
docker/run_activespace.sh -m aam docker/validate_aam_below_ground.py
```

Adapter implementation and pipeline notes: [`docs/aam_integration_notes.md`](../docs/aam_integration_notes.md).

AAM terrain (`Input_Data/AAM/terrain/{mic}/scenario.elv`) is cached on disk; reruns skip
`write_terrain` when ELV is newer than the parent DEM and `terrain_cache.json` matches.
Prediction cache: `Output_Data/aam/predictions/` (NMSim: `Output_Data/nmsim/predictions/`).
Run log (terrain, batches, summaries): `Output_Data/aam/active_space.log` — points at
`runs/{job}/scenario.inp` and `scenario.txt`, not full deck contents.

Model-scoped outputs live under `Output_Data/nmsim/` and `Output_Data/aam/` (active spaces,
prediction cache, precision-recall plots). NMSim readers still fall back to legacy
`Output_Data/TIG_TIS/` and flat `Output_Data/ACTIVESPACES/` when present.

## Optional: `validate_active_space.py` (debug harness)

`docker/validate_active_space.py` is an integration/debug harness — **not** the production
fit path. Use `generate_active_space.py` (and the 3D batch workflow) for real runs and
project-level `fits.csv` output.

```bash
# Quick NMSim integration check
docker/run_activespace.sh docker/validate_active_space.py \
  -u DENA -s TRLA -y 2025 --gains 0 --altitude 1000 --density 10 --heading 0

# AAM debug: annotation gain sweep + PR plot (both -m and --model required)
docker/run_activespace.sh -m aam docker/validate_active_space.py --model aam \
  -u DENA -s TRLA -y 2025 --fit --omni-min 0 --omni-max 2 \
  --altitude 1000 --density 10 --heading 0
```

The `--fit` debug path loads `*saved_annotations*.geojson`, sweeps gains, and can write
site-local diagnostics. For canonical fits, use `fit_3d_active_space.py` → project
`fits.csv` with a `Model` column (see [`docs/aam_integration_notes.md`](../docs/aam_integration_notes.md)).

Model selection is a **CLI flag** on `run_activespace.sh` (`-m nmsim|aam`), not a separate
config file — the pipeline still uses `container.config` with `project.nmsim` only.

| | NMSim (default) | AAM |
|-|-------|-----|
| Container | `docker/run_activespace.sh ...` | `docker/run_activespace.sh -m aam ...` |
| Pipeline script | (default) | add `--model aam` on `generate_active_space.py` etc. |
| Staging | `docker/stage_nmsim_runtime.sh` | `docker/stage_aam_runtime.sh` |
| Local dir | `vendor/nmsim-runtime/` | `vendor/aam-runtime/` |
| Mount | `/opt/nmsim` | `/opt/aam` |
| Shim | `/usr/local/bin/nord2000` | `/usr/local/bin/aam` |
