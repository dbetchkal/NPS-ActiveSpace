# Mac / Linux: Docker + Wine for acoustic models

Canonical setup for running **Windows-only** NMSim (and optionally AAM binaries) on Mac/Linux.
Python 3.12 + GDAL run inside the container; binaries execute through Wine. Expect slower
runs than native Windows.

**Prerequisites:** [Docker Desktop](https://docs.docker.com/get-started/get-docker/).
On Apple Silicon, enable Rosetta for amd64 emulation.

NMSim and AAM binaries are **not redistributable** (NPS internal). Do not commit them.
See [vendor/nmsim-runtime/README.md](../vendor/nmsim-runtime/README.md) for what gets staged locally.

## Setup (Wine + binaries)

From the repo root. **`docker/smoke.sh` does not use `container.config`** — only staged runtimes
and a built image.

```bash
docker/stage_nmsim_runtime.sh /path/to/NMSim-install   # copies into vendor/nmsim-runtime/
docker/build.sh                                        # first build ~13 min (linux/amd64)
docker/smoke.sh                                        # Nord2000batch.exe via Wine shim
```

Optional AAM binary check (same container; not wired into the Python pipeline in this PR):

```bash
docker/stage_aam_runtime.sh /path/to/AAM_v3_dec2020    # needs AAM_3.0.0.exe, NCfiles/, noisecon.inp
docker/smoke.sh aam
```

## Run (NMSim pipeline)

After smoke passes, copy a container config and run generation inside Docker.
Paths must be absolute **`/repo/...`** (repo is bind-mounted at `/repo`).

Set **`project.nmsim`** to the Wine shim **`/usr/local/bin/nord2000`**, not the runtime mount
at **`/opt/nmsim`** (the shim `cd`s there for `RND/`).

```bash
cp nps_active_space/config/container_example.config nps_active_space/config/container.config

docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000
```

**Host vs container:** ground-truthing, fit, and viz use a local venv and `-e DENA_example`
(or your config) — see [example_data/README.md](../example_data/README.md) and root
[README.md](../README.md). Windows NMSim install is unchanged.

### Overrides

| Variable | When | Purpose |
|----------|------|---------|
| `NMSIM_RUNTIME` | `stage_nmsim_runtime.sh` and `run_activespace.sh` | Staged copy location (default `vendor/nmsim-runtime/`). Use the **same** path for both. |
| `AAM_RUNTIME` | `stage_aam_runtime.sh` and `run_activespace.sh -m aam` | Same pattern for AAM (default `vendor/aam-runtime/`). |
| `DATA_DRIVE` | `run_activespace.sh` only | Optional host directory mounted read-only at `/data`. |

### Reference

| | NMSim (default) | AAM (optional) |
|-|-------|-----|
| Setup check | `docker/smoke.sh` | `docker/smoke.sh aam` |
| Active-space pipeline | `run_activespace.sh … generate_active_space.py` | follow-up PRs |
| Stage into | `docker/stage_nmsim_runtime.sh` | `docker/stage_aam_runtime.sh` |
| Default local dir | `vendor/nmsim-runtime/` | `vendor/aam-runtime/` |
| Container mount | `/opt/nmsim` | `/opt/aam` |
| Shim in config / PATH | `/usr/local/bin/nord2000` | `/usr/local/bin/aam` |
