# AAM propagation adapter

AAM 3.0 ([aam-translator](https://github.com/elliott-ruebush/aam-translator)) behind `PropagationModel`. One fixed mic, many source points; input decks use `COMPUTEPOI` with one POI. Reading POI files: [reading_aam_output.md](https://github.com/elliott-ruebush/aam-translator/blob/main/docs/reading_aam_output.md).

Docker: [docker/README.md](../../../docker/README.md).

## Module map

| Module | Role |
|--------|------|
| `model.py` | `AamPropagationModel` — `prepare_site`, chunking, FPA split retry |
| `config.py` | `AAM_CHUNK_SIZE` / `resolve_aam_chunk_size`, subprocess timeout |
| `track.py` | Mesh snake order, single-point track padding |
| `terrain_cache.py` | DEM warp, ELV/IMP cache, `ensure_aam_terrain` |
| `terrain_sampling.py` | ELV bilinear sample, below-ground filter, clear-hop packing |
| `terrain.py` | Re-exports cache + sampling (stable import path) |
| `source.py` | Omni NetCDF cache, vendor `NCfiles/`, subprocess env |
| `run_batch.py` | Stage `scenario.inp`, Wine/native subprocess, read POI |
| `run_log.py` | Site log, stderr summaries, work-dir hashing |
| `output.py` | POI history → NMSim-shaped prediction DataFrame |

## On-disk layout

Terrain cache: `Input_Data/aam/terrain/{mic}/` (`.ELV`, `.IMP`, `terrain_cache.json`). Skips `write_terrain` when ELV is newer than the DEM and metadata matches.

Omni NetCDF: `Input_Data/aam/NCfiles/` (generated, not committed). Other paths under `Output_Data/aam/` — see `nps_active_space.utils.paths`.

NMSim legacy flat `Output_Data/ACTIVESPACES/` and `Output_Data/TIG_TIS/` are still read for older NMSim-only sites; new runs use `Output_Data/nmsim/`. Remove stale `TIG_TIS/` if an old tree was mixed with the new layout.

## Omni NetCDF

`.src`/`.avg` omni files become per-token NetCDF on first use (`OMNI_000`, `OMNIM100`, …). Template: `OMNI_200.nc` on the AAM install (`AAM_NC` if `Bin\NCfiles` is empty). Each run copies one omni into job `NCfiles/` for `ROTOR_NOISE` / `AAM_NC`. Using the whole site cache can cause `Same Profile data found in multiple spheres` when `--omni-max` > 0. See `stage_run_ncfiles` in `source.py`.

## Track batching and audibility

See `AamPropagationModel.predict` in `model.py`. Pipeline: `track.order_source_pts_for_track` → ELV vertex filter → clear-hop packs (`terrain_sampling`) → `config.resolve_aam_chunk_size()` (default 400) → `run_batch.execute_aam_batch`. Generator batch cap remains `DEFAULT_MAX_POINTS_PER_PREDICT` (4000) in `protocol.py`.

No 12.5 kHz in AAM POI output — `poi_history_to_predictions_df` uses `NaN` there (`output.py`).
