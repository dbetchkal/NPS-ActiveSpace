# AAM runtime (local, not in git)

The Windows AAM binary is not in this repo. Stage a copy into `vendor/aam-runtime/` before Docker runs — [docker/README.md](../../docker/README.md).

Typical install includes `AAM_3.0.0.exe`, `NCfiles/` (template NetCDF), and `noisecon.inp` for smoke tests. `docker/stage_aam_runtime.sh` copies those from your vendor tree.

## Redistribution / licensing

AAM is NPS internal software, not open source. Do not commit binaries to a public repo without NPS approval.
