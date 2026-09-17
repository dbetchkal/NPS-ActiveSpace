# NMSim runtime (local, not in git)

The Windows NMSim binary and its runtime files are **not committed to this repository**.
Populate **`vendor/nmsim-runtime/`** before Docker+Wine runs on Mac/Linux.

**Setup steps:** [docker/README.md](../../docker/README.md) (stage → build → smoke → pipeline).

## Size

Minimal runtime for omni-source active-space runs (~**10 MB**):

| Component | Size |
|-----------|------|
| `Nord2000batch.exe` | ~1.3 MB |
| `Nord2000-DLL.dll` | ~1.2 MB |
| `RNM_SRC.dll` | ~512 KB |
| `netcdf.dll` | ~204 KB |
| `RND/` (character files, `directories.ini`, …) | ~7 MB |

A full vendor install with example cases (`FortTiCase`, `Sources/`, …) is ~11 MB.

## Staging

From a machine that has the NMSim install (e.g. NPS data drive or Windows box):

```bash
docker/stage_nmsim_runtime.sh /path/to/NMSim
```

To keep the staged copy outside `vendor/nmsim-runtime/`, set **`NMSIM_RUNTIME`** to the
same directory when staging and when calling **`docker/run_activespace.sh`** (see
[docker/README.md](../../docker/README.md)).

## Redistribution / licensing

NMSim (`Nord2000batch.exe`) is **NPS internal / government software**, not open source.
Do **not** commit the binaries to a public GitHub repo unless NPS legal explicitly approves
redistribution. Size (~10 MB) is git-friendly; **licensing is the constraint**.

Omni tuning sources (`.src`/`.avg`) used by the pipeline live in-repo under
`nps_active_space/data/tuning/` and are separate from this executable runtime.
