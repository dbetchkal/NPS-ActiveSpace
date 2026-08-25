# Temporary: Windows AAM vs NMSim at 1500 m (DENATRLA 2025)

Scratch notes for a native Windows compare. NMSim’s strongest per-layer F1 on this site was **1500 m (0.60)** in `DENATRLA2025_commands_output.csv`. Use that altitude so the compare is apples-to-apples.

Do **not** set `ACOUSTIC_MODEL`, `AAM_PARALLEL_N`, or `AAM_CHUNK_SIZE` in your shell. Those are Docker/script internals. User settings live in `nps_active_space/config/<env>.config`.

## One-time setup

1. Python 3.12 venv at the repo root:

   ```bat
   py -3.12 -m venv .venv
   .venv\Scripts\activate
   pip install -e ".[dev,aam]"
   ```

2. Copy config and fill paths (same pattern as NMSim):

   ```bat
   copy nps_active_space\config\template.config nps_active_space\config\windows.config
   ```

   Required `[project]` / `[data]` fields for this run:

   | Key | Windows example |
   |-----|-----------------|
   | `[project] dir` | `C:\...\example_data\site_projects` |
   | `[project] aam` | `C:\path\to\AAM_v3\AAM_3.0.0.exe` |
   | `[project] nmsim` | path to `Nord2000batch.exe` (only if you also run NMSim) |
   | `[data] dem` | `...\DENATRLA\Input_Data\01_ELEVATION\elevation_m_nad83_utm6.tif` |

   `NCfiles\` must sit **next to** `AAM_3.0.0.exe` (the code sets noise-database env vars on the subprocess only).

3. Ambience pickle already used on Mac:

   `example_data\site_projects\DENATRLA\Output_Data\AMBIENCE\DENATRLA2025_ambience.pkl`

## Smoke (one gain, coarse mesh) — minutes

From repo root, venv active:

```bat
python -u nps_active_space\scripts\generate_active_space.py ^
  -e windows --model aam -u DENA -s TRLA -y 2025 -l 1500 ^
  --omni-min 0 --omni-max 0 --density 10 --headings 0 --cleanup ^
  -a example_data\site_projects\DENATRLA\Output_Data\AMBIENCE\DENATRLA2025_ambience.pkl ^
  --results-out denatrla_1500m_aam_smoke.json
```

Same flags without `--model aam` for NMSim, if you want a matching smoke.

## Full single-layer (matches existing NMSim 1500 m batch flags)

Same as `DENATRLA2025_commands.txt` 1500 m line, plus `--model aam`:

```bat
python -u nps_active_space\scripts\generate_active_space.py ^
  -e windows --model aam -u DENA -s TRLA -y 2025 -l 1500 ^
  --headings 0 --omni-min 0 --omni-max 2 --cleanup ^
  -a example_data\site_projects\DENATRLA\Output_Data\AMBIENCE\DENATRLA2025_ambience.pkl ^
  --results-out denatrla_1500m_aam_full.json
```

Default mesh density is 48 (same as the NMSim batch). Native AAM should use `cpu_count - 1` omni workers automatically.

## What to copy back

- `denatrla_1500m_aam_smoke.json` and/or `denatrla_1500m_aam_full.json` (`F1` field)
- Last ~80 lines of console if the script exits non-zero
- Any `[aam-run] FAILED` line that is **not** followed by `retrying halves` or `skipping 1 point` (those are expected isolation retries)

Do **not** run `fit_3d_active_space` until several altitude layers exist. This note is one layer only.

## Expected log noise

`[aam-filter] filtered N/M below AAM terrain` and binary-split `FAILED` / `retrying halves` are expected. The job should still finish and write geojson under `Output_Data\aam\ACTIVESPACES\DENATRLA2025_1500m\`.
