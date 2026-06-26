# NMSim runtime (local, not in git)

The Windows NMSim binary and its runtime files are **not committed to this repository**.
Populate this directory before running the Docker+Wine pipeline on Mac/Linux.

## Size

Minimal runtime needed for omni-source active-space runs (~**10 MB**):

| Component | Size |
|-----------|------|
| `Nord2000batch.exe` | ~1.3 MB |
| `Nord2000-DLL.dll` | ~1.2 MB |
| `RNM_SRC.dll` | ~512 KB |
| `netcdf.dll` | ~204 KB |
| `RND/` (character files, `directories.ini`, …) | ~7 MB |

A full vendor install with example cases (`FortTiCase`, `Sources/`, …) is ~11 MB.

## Setup

From a machine that has the NMSim install (e.g. the NPS data drive or a Windows box):

```bash
docker/stage_nmsim_runtime.sh /path/to/NMSim
# or: NMSIM_SOURCE=/path/to/NMSim docker/stage_nmsim_runtime.sh
```

Then build and run:

```bash
docker/build.sh
docker/run_activespace.sh docker/validate_active_space.py -u DENA -s TRLA -y 2025 --gains 0
```

Override the runtime location with `NMSIM_RUNTIME=/other/path docker/run_activespace.sh …`.

## Redistribution / licensing

NMSim (`Nord2000batch.exe`) is **NPS internal / government software**, not open source.
Do **not** commit the binaries to a public GitHub repo unless NPS legal explicitly approves
redistribution. Size (~10 MB) is git-friendly; **licensing is the constraint**.

Omni tuning sources (`.src`/`.avg`) used by the pipeline live in-repo at
`nps_active_space/data/tuning/` and are separate from the NMSim executable runtime.
