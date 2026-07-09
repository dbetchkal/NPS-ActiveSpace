# Mac / Linux: Docker + Wine for NMSim

NMSim is Windows-only. On Mac/Linux, we can run the full pipeline inside a container that
includes Python 3.12 + GDAL and executes the `Nord2000batch.exe` Windows binary through Wine.

Note that there are some performance implications here given the additional layers of indirection to run the docker setup.
# TODO: do some more testing to precisely quantify the performance difference if we can do that easily.

**Prerequisites:** Docker Desktop; on Apple Silicon enable Rosetta for amd64 emulation.
# TODO: should we have instructions on how to install this? Maybe at least a link to the docker install page.

## One-time setup

```bash
# 1) Stage the NMSim runtime locally (~10 MB; not in git — see vendor/nmsim-runtime/README.md)
docker/stage_nmsim_runtime.sh /path/to/NMSim-install

# 2) Build the image (~13 min first time; installs requirements.container.txt)
docker/build.sh

# 3) Config (if not already present)
cp nps_active_space/config/container_example.config nps_active_space/config/container.config
```

## Run

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

Use `-e container` with absolute `/work/...` paths in config. Native viz and ground-truthing on the host use `-e DENA_example` (or your own config) and a local venv — see [example_data/README.md](../example_data/README.md).

Windows setup is unchanged — see root [README.md](../README.md) Installation.
