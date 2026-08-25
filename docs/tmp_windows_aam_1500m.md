# Windows: AAM vs NMSim at 1500 m (DENATRLA 2025)

Native Windows compare of one active-space layer. Use **1500 m** (best NMSim F1 on the Mac example was 0.60) unless the Windows annotations change the altitude mix.

Do not set `ACOUSTIC_MODEL`, `AAM_PARALLEL_N`, or `AAM_CHUNK_SIZE`. Put paths in `nps_active_space/config/windows.config`.

## Setup (once)

```bat
py -3.12 -m venv .venv
.venv\Scripts\activate
pip install -e ".[dev,aam]"
copy nps_active_space\config\template.config nps_active_space\config\windows.config
```

In `windows.config` set:

- `[project] dir` — folder that contains `DENATRLA\`
- `[project] aam` — `AAM_3.0.0.exe` (`NCfiles\` must sit next to it)
- `[project] nmsim` — `Nord2000batch.exe`
- `[data] dem` — site DEM GeoTIFF

## Same inputs on both runs

Keep a single `DENATRLA2025_saved_annotations.geojson` in the site directory root. Do not mix Mac example tracks with the Windows set; regenerate NMSim here.

Use the same ambience pickle (`-a` below). If more than one `*saved_annotations*.geojson` exists, keep `--annotation-file` on both commands.

## Run

From the repo root, venv active. Only `--model` and `--results-out` differ.

**AAM**

```bat
python -u nps_active_space\scripts\generate_active_space.py ^
  -e windows --model aam -u DENA -s TRLA -y 2025 -l 1500 ^
  --headings 0 --omni-min 0 --omni-max 2 --cleanup ^
  --annotation-file DENATRLA2025_saved_annotations.geojson ^
  -a example_data\site_projects\DENATRLA\Output_Data\AMBIENCE\DENATRLA2025_ambience.pkl ^
  --results-out denatrla_1500m_aam.json
```

**NMSim**

```bat
python -u nps_active_space\scripts\generate_active_space.py ^
  -e windows -u DENA -s TRLA -y 2025 -l 1500 ^
  --headings 0 --omni-min 0 --omni-max 2 --cleanup ^
  --annotation-file DENATRLA2025_saved_annotations.geojson ^
  -a example_data\site_projects\DENATRLA\Output_Data\AMBIENCE\DENATRLA2025_ambience.pkl ^
  --results-out denatrla_1500m_nmsim.json
```

Smoke (minutes): add `--density 10 --omni-max 0` to both.

Outputs: `Output_Data\aam\ACTIVESPACES\DENATRLA2025_1500m\` and `Output_Data\nmsim\ACTIVESPACES\`. Do not run `fit_3d_active_space` for a single altitude.

## Compare

- **F1 / gain:** the two JSON files (`F1`, `1/3rd Octave Gain (F1)`).
- **Shape:** `python nps_active_space\scripts\viz.py DENATRLA2025 -e windows -s -a --compare -g <gain>`  
  NMSim orange, AAM cyan. Uncheck non-1500 m layers. Viz `--annotation-file` needs a full path.

## If it fails

Send both JSON files and the last ~80 console lines.

Expected, not a crash: `[aam-filter] filtered N/M below AAM terrain`, `[aam-predict] isolating …` / `skipped 1 point`. Fortran stacks: `Output_Data\aam\runs\<job>\aam_stderr.txt`.
