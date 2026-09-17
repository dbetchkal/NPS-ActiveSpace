# Docker + Wine (Mac/Linux)

NMSim is Windows-only. On Mac/Linux it runs in a container (Python 3.12 + GDAL) through Wine—slower than a native Windows install.

Install [Docker Desktop](https://docs.docker.com/get-started/get-docker/). On Apple Silicon, turn on Rosetta for amd64 images.

Binaries are not in git. See [vendor/nmsim-runtime/README.md](../vendor/nmsim-runtime/README.md).

## Setup

Run from the repo root. `docker/smoke.sh` does not use `container.config`.

```bash
docker/stage_nmsim_runtime.sh /path/to/NMSim-install
docker/build.sh
docker/smoke.sh
```

Optional AAM smoke: `docker/stage_aam_runtime.sh` + `docker/smoke.sh aam`.

## Run

```bash
cp nps_active_space/config/container_example.config nps_active_space/config/container.config

docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000
```

Use `/repo/...` paths in the config. Mounts and shims: [container_example.config](../nps_active_space/config/container_example.config).

### Overrides

| Variable | Purpose |
|----------|---------|
| `NMSIM_RUNTIME` | Staged NMSim tree (default `vendor/nmsim-runtime/`). Same path for `stage_nmsim_runtime.sh` and `run_activespace.sh`. |
| `DATA_DRIVE` | Optional host directory mounted read-only at `/data`. |
