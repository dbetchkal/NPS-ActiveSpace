# NMSim runtime (local, not in git)

The Windows NMSim binary is not in this repo. Stage a copy into `vendor/nmsim-runtime/` before Docker runs — [docker/README.md](../../docker/README.md).

## Size

Typical staged copy (~10 MB):

| Component | Size |
|-----------|------|
| `Nord2000batch.exe` | ~1.3 MB |
| `Nord2000-DLL.dll` | ~1.2 MB |
| `RNM_SRC.dll` | ~512 KB |
| `netcdf.dll` | ~204 KB |
| `RND/` (character files, `directories.ini`, …) | ~7 MB |

A full vendor install with example cases (`FortTiCase`, `Sources/`, …) is ~11 MB.

## Redistribution / licensing

NMSim (`Nord2000batch.exe`) is NPS internal software, not open source. Do not commit binaries to a public repo without NPS approval.

Omni `.src`/`.avg` tuning files live in `nps_active_space/propagation_model/nmsim/data/tuning/` (not part of the staged executable tree).
