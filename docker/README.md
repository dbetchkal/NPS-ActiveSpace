# Docker + Wine (Mac/Linux)

Windows-only NMSim runs in a container (Python 3.12 + GDAL) through Wine. Expect slower runs than native Windows.

**Prerequisites:** [Docker Desktop](https://docs.docker.com/get-started/get-docker/). On Apple Silicon, enable Rosetta for amd64 emulation.

Binaries are not in git — see [vendor/nmsim-runtime/README.md](../vendor/nmsim-runtime/README.md).

## Setup

From the repo root. `docker/smoke.sh` does not use `container.config`.

```bash
docker/stage_nmsim_runtime.sh /path/to/NMSim-install
docker/build.sh
docker/smoke.sh
```

Optional: `docker/stage_aam_runtime.sh` + `docker/smoke.sh aam` (Wine launch check only).

## Run

```bash
cp nps_active_space/config/container_example.config nps_active_space/config/container.config

docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000
```

Use absolute `/repo/...` paths in the config. Mounts and `project.nmsim` are described in [container_example.config](../nps_active_space/config/container_example.config).

### Overrides

| Variable | Purpose |
|----------|---------|
| `NMSIM_RUNTIME` | Staged NMSim tree (default `vendor/nmsim-runtime/`). Use the same path for `stage_nmsim_runtime.sh` and `run_activespace.sh`. |
| `DATA_DRIVE` | Optional host directory mounted read-only at `/data`. |
