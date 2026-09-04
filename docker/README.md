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

The image already installs `[aam]` on build. For a host venv: `pip install -e ".[dev,aam]"` (GitHub pin in `pyproject.toml`).

Active space is one fixed microphone plus many source positions, so the adapter uses AAM `COMPUTEPOI` with 1 POI and N track points. General AAM I/O and runtime gotchas live in [aam-translator](https://github.com/elliott-ruebush/aam-translator) (`docs/reading_aam_output.md`).

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

### On-disk layout

AAM terrain is cached; reruns skip `write_terrain` when the ELV is newer than the parent DEM and `terrain_cache.json` matches. The run log points at on-disk artifacts (`scenario.inp`, `scenario.txt`), not full deck contents.

| Path | Role |
|------|------|
| `Input_Data/aam/terrain/{mic}/` | Cached `.ELV` / `.IMP` + `terrain_cache.json` (rebuild when DEM/AOI changes) |
| `Input_Data/aam/NCfiles/` | Generated per-omni NetCDF (`OMNI_000.nc`, …); runtime cache, not committed |
| `Output_Data/aam/predictions/` | Incremental spectral cache CSVs |
| `Output_Data/aam/runs/{job}/` | Per-batch scratch (`.inp`, `.POI`, `scenario.txt`). Job folders are hashed (`x{hash12}`) — full names in `active_space.log` |
| `Output_Data/aam/active_space.log` | Append-only run log |
| `Output_Data/aam/ACTIVESPACES/` | Final GeoJSON |

Legacy: `Input_Data/AAM/terrain/{mic}/` and `Input_Data/AAM/terrain_{mic}/` still read if present. Stale `Output_Data/TIG_TIS/*.csv` from early AAM runs should be deleted (they mixed with NMSim caches). New NMSim runs use `Output_Data/nmsim/predictions/`. NMSim readers still fall back to legacy `Output_Data/TIG_TIS/` and flat `Output_Data/ACTIVESPACES/` when present.

### Omni NetCDF staging

NMSim omni ladder files (`O_±XXX.src` + sibling `.avg`) are converted to AAM NetCDF at first use:

- **Template:** vendor `OMNI_200.nc` on the AAM install (read-only). Set `AAM_NC` if `Bin\NCfiles` is an empty stub.
- **Cache:** `Input_Data/aam/NCfiles/{token}.nc` where `O_+000` → `OMNI_000`, `O_-100` → `OMNIM100`.
- Each AAM run gets `work_dir/NCfiles/` with **only** that job's omni NetCDF on `ROTOR_NOISE` / `AAM_NC`. Pointing AAM at the whole site cache triggers `Same Profile data found in multiple spheres` when omni `--max` > 0.
- No post-hoc dB offset; levels come from the generated NetCDF.

Omni groups use the same `mp.Pool` as NMSim. Docker+Wine is capped at 2 workers (`AAM_PARALLEL_N` in `run_activespace.sh`). Do not set that in user config.

### Track batching

AAM aborts an entire `ONE TRACK` if any vertex or interpolated hop is below the ELV surface (often empty `.POI`). The adapter avoids those decks rather than retrying (`AamPropagationModel.predict`):

- Filter vertices against the ELV grid (`split_below_aam_terrain`)
- Snake the lattice so hops stay local (`_order_source_pts_for_track`)
- Split hops that still clip terrain (`split_safe_aam_track_runs`), then chunk at `AAM_CHUNK_SIZE` (default 400)
- Pad a leftover 1-vertex track ~1 m (`_pad_single_point_track`) — AAM 3.0.0 crashes on a single vertex
- On remaining failure, skip that chunk (points treated inaudible). Do **not** bisect below-ground aborts. Fortran FPA-bounds errors on a long high-altitude track **are** retried by halving.
- Run dirs are hashed because AAM's Fortran `FILENAME` buffer is 140 chars (`short_aam_work_dir_name`)

### Audibility / 12.5 kHz

AAM POI has no 12.5 kHz band; the adapter writes `NaN` (`poi_history_to_predictions_df`). `spectrum_is_audible` still compares through 12.5 kHz against `max(ambience, hearing threshold)`, so that band never tips audibility from AAM.

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

The `--fit` debug path loads `*saved_annotations*.geojson`, sweeps gains, and can write site-local diagnostics. Canonical 3D fits follow the production path in [`nps_active_space/scripts/README.md`](../nps_active_space/scripts/README.md): `plot_altitudes.py` → `generate_3d_active_space.py --model` → `generate_active_space_batch.py` → `fit_3d_active_space.py --model` → project `fits.csv`. That project `fits.csv` has a `Model` column (`nmsim` / `aam`), keyed by Designator + Model + Altitude_m. Per-site `{site}/fits.csv` from this debug harness is not authoritative.

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
