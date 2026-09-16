# Mac / Linux: Docker + Wine for acoustic models

NMSim (and optionally AAM) are Windows-only. On Mac/Linux, run them in a container
(Python 3.12 + GDAL) through Wine. Extra layers can be slower than native Windows.

**Prerequisites:** [Docker Desktop](https://docs.docker.com/get-started/get-docker/).
On Apple Silicon, enable Rosetta for amd64 emulation.

NMSim and AAM binaries are **not redistributable** (NPS internal). Do not commit them.

## Setup

```bash
docker/stage_nmsim_runtime.sh /path/to/NMSim-install   # ~10 MB into vendor/nmsim-runtime/
docker/build.sh                                        # ~13 min first time
cp nps_active_space/config/container_example.config nps_active_space/config/container.config
docker/smoke.sh                                        # NMSim on DENATRLA example data
```

`container.config` must set `project.nmsim` to the Wine shim (`/usr/local/bin/nord2000`).
`/opt/nmsim` is the runtime bind-mount, not the executable path.

Optional AAM (binary-only until later PRs wire it into the pipeline):

```bash
docker/stage_aam_runtime.sh /path/to/AAM_v3_dec2020    # needs AAM_3.0.0.exe, NCfiles/, noisecon.inp
docker/smoke.sh aam
```

## Run

```bash
docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py \
  -e container -u DENA -s TRLA -y 2025 -l 1000
```

Ground-truthing, fit, and viz stay on the host (`-e DENA_example`) — see
[example_data/README.md](../example_data/README.md). Windows install is unchanged
(root [README.md](../README.md)).

Optional: `DATA_DRIVE=/Volumes/NPS_ADSB_Data docker/run_activespace.sh ...` mounts `/data`.
Override runtime with `NMSIM_RUNTIME=` / `AAM_RUNTIME=`. Use `-e container` with
absolute `/repo/...` paths.

| | NMSim (default) | AAM |
|-|-------|-----|
| Setup check | `docker/smoke.sh` | `docker/smoke.sh aam` |
| Run | `docker/run_activespace.sh …` | `docker/run_activespace.sh -m aam …` |
| Staging | `docker/stage_nmsim_runtime.sh` | `docker/stage_aam_runtime.sh` |
| Local dir | `vendor/nmsim-runtime/` | `vendor/aam-runtime/` |
| Mount | `/opt/nmsim` | `/opt/aam` |
| Shim | `/usr/local/bin/nord2000` | `/usr/local/bin/aam` |
