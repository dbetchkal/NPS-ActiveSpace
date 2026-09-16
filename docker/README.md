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

## AAM

`docker/build.sh` installs the `[aam]` extra. On the host: `pip install -e ".[dev,aam]"`.

Stage AAM, then use `-m aam` and `--model aam`:

```bash
docker/stage_aam_runtime.sh /path/to/AAM_v3_dec2020
docker/run_activespace.sh -m aam nps_active_space/scripts/generate_active_space.py \
  -e container --model aam -u DENA -s TRLA -y 2025 -l 1000
```

[propagation_model/aam/README.md](../nps_active_space/propagation_model/aam/README.md) · [scripts/README.md](../nps_active_space/scripts/README.md)

Under Wine, `run_activespace.sh` sets `AAM_PARALLEL_N=2` (do not set it in your config).
