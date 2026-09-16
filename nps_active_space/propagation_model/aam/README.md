# AAM propagation adapter

Runs [AAM 3.0](https://github.com/elliott-ruebush/aam-translator) under the shared `PropagationModel` interface. Active space is one fixed microphone and many source positions; decks use AAM `COMPUTEPOI` with one POI and N track points. POI output semantics: [aam-translator `docs/reading_aam_output.md`](https://github.com/elliott-ruebush/aam-translator/blob/main/docs/reading_aam_output.md).

Docker staging and `-m aam` / `--model aam`: [docker/README.md](../../../docker/README.md).

## On-disk layout

Terrain is cached under `Input_Data/aam/terrain/{mic}/` (`.ELV`, `.IMP`, `terrain_cache.json`). Reruns skip `write_terrain` when the ELV is newer than the parent DEM and the cache metadata matches.

Generated omni NetCDF lives in `Input_Data/aam/NCfiles/` (not committed). Predictions, per-batch scratch under `Output_Data/aam/runs/{job}/`, `active_space.log`, and final GeoJSON under `Output_Data/aam/ACTIVESPACES/` follow the path helpers in `nps_active_space.utils.paths`.

Legacy `Input_Data/AAM/` paths and flat `Output_Data/ACTIVESPACES/` are still read when present. NMSim caches moved to `Output_Data/nmsim/predictions/`; delete stale mixed `Output_Data/TIG_TIS/` from early AAM experiments.

## Omni NetCDF

NMSim omni `.src`/`.avg` pairs convert to per-token NetCDF on first use (`OMNI_000`, `OMNIM100`, …). Template: vendor `OMNI_200.nc` on the AAM install (set `AAM_NC` if `Bin\NCfiles` is empty). Each run uses a job-local `NCfiles/` with **one** omni file on `ROTOR_NOISE` / `AAM_NC` — pointing at the full site cache can trigger `Same Profile data found in multiple spheres` when omni `--max` > 0. See `stage_run_ncfiles` in `source.py`.

## Track batching and audibility

Track splitting, terrain filtering, chunk size, and retry policy are documented on `AamPropagationModel.predict` in `model.py`. The generator may pass up to `DEFAULT_MAX_POINTS_PER_PREDICT` points per call; AAM chunks internally at `AAM_CHUNK_SIZE` (default 400).

AAM POI has no 12.5 kHz band; `poi_history_to_predictions_df` writes `NaN` there so audibility is not tipped by that band (`output.py`).
