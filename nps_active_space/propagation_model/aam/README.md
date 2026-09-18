# AAM propagation adapter

AAM 3.0 ([aam-translator](https://github.com/elliott-ruebush/aam-translator)) behind `PropagationModel`. One fixed mic, many source points; input decks use `COMPUTEPOI` with one POI. Reading POI files: [reading_aam_output.md](https://github.com/elliott-ruebush/aam-translator/blob/main/docs/reading_aam_output.md).

Docker: [docker/README.md](../../../docker/README.md).

## On-disk layout

Terrain cache: `Input_Data/aam/terrain/{mic}/` (`.ELV`, `.IMP`, `terrain_cache.json`). Skips `write_terrain` when ELV is newer than the DEM and metadata matches.

Omni NetCDF: `Input_Data/aam/NCfiles/` (generated, not committed). Other paths under `Output_Data/aam/` — see `nps_active_space.utils.paths`.

Old `Input_Data/AAM/` layouts and flat `Output_Data/ACTIVESPACES/` still work. NMSim predictions live under `Output_Data/nmsim/predictions/`. Remove stale `Output_Data/TIG_TIS/` if you mixed early AAM and NMSim runs.

## Omni NetCDF

`.src`/`.avg` omni files become per-token NetCDF on first use (`OMNI_000`, `OMNIM100`, …). Template: `OMNI_200.nc` on the AAM install (`AAM_NC` if `Bin\NCfiles` is empty). Each run copies one omni into job `NCfiles/` for `ROTOR_NOISE` / `AAM_NC`. Using the whole site cache can cause `Same Profile data found in multiple spheres` when `--omni-max` > 0. See `stage_run_ncfiles` in `source.py`.

## Track batching and audibility

See `AamPropagationModel.predict` in `model.py` and mesh ordering in `track.py`. Up to `DEFAULT_MAX_POINTS_PER_PREDICT` points per call; internal chunk size `AAM_CHUNK_SIZE` (default 400).

No 12.5 kHz in AAM POI output — `poi_history_to_predictions_df` uses `NaN` there (`output.py`).
