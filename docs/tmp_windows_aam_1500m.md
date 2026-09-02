# Windows: AAM vs NMSim at 1500 m (DENATRLA 2025)

Native Windows compare of one active-space layer. Use **1500 m** (best NMSim F1 on the Mac example was 0.60) unless your annotations change the altitude mix.

Do not set `ACOUSTIC_MODEL`, `AAM_PARALLEL_N`, or `AAM_CHUNK_SIZE`. Paths live in `nps_active_space/config/AK.config`.

## Setup (once)

```bat
py -3.12 -m venv .venv
.venv\Scripts\activate
pip install -e ".[dev,aam]"
copy nps_active_space\config\AK_example.config nps_active_space\config\AK.config
```

**Mac-parity / platform compare:** use `AK_example.config` as-is (points at `example_data\`). Set only `[project] nmsim` and `[project] aam`.

**Production site on `T:\`:** keep your existing `[project] dir` and data paths; do not mix with Mac `example_data`.

In `AK.config` set:

- `[project] dir` — `example_data\site_projects` for parity, or your folder that contains `DENATRLA\`
- `[project] aam` — `AAM_3.0.0.exe` (`NCfiles\` next to the exe, or under parent `AAM\`; override with env `AAM_NC`)
- `[project] nmsim` — `Nord2000batch.exe`
- `[data] nvspl_archive` — `example_data\nvspl_archive` for parity

Generation reads `Input_Data\01_ELEVATION\elevation_m_nad83_utm*.tif` (and `.flt`/`.hdr`) from `project_setup`. `[data] dem` is only the parent raster for that setup step.

## Annotations (committed)

For the example slice, annotations are in git:

`example_data\site_projects\DENATRLA\DENATRLA2025_saved_annotations.geojson`

The parity bats pass `--annotation-file DENATRLA2025_saved_annotations.geojson` so only that file is used. Keep a single `*saved_annotations*.geojson` in the site root if you omit the flag.

## Mac-parity smoke (platform compare)

The full `docs\tmp_windows_aam_1500m.bat` is density **48**, omni **0–2**. For the same job Mac Docker ran, use:

```bat
docs\tmp_windows_aam_macparity.bat
```

That is **example_data**, 1500 m, density **10**, one omni (`+000`), heading 0, AAM then NMSim, live NVSPL. Same generate flags as:

```bash
docs/tmp_docker_aam_macparity.sh
```

Both write `docker/denatrla_1500m_macparity_{aam,nmsim,metrics}.json`.

## Run

From the repo root, venv optional (the script activates `.venv` if present):

```bat
docs\tmp_windows_aam_1500m.bat
docs\tmp_windows_aam_1500m.bat smoke
```

Runs use `-e AK`. JSON results land under `docker\` for the bats above.

Outputs: `Output_Data\aam\ACTIVESPACES\DENATRLA2025_1500m\` and `Output_Data\nmsim\ACTIVESPACES\` under the site dir from `[project] dir`.

## Compare

- **F1 / gain:** the two JSON files (`F1`, `1/3rd Octave Gain (F1)`).
- **Shape:** `python nps_active_space\scripts\viz.py DENATRLA2025 -e AK -s -a --compare -g <gain>`  
  NMSim orange, AAM cyan. Uncheck non-1500 m layers. Viz `--annotation-file` needs a full path.

## If it fails

**`AAM NCfiles/ not found`:** On the Windows host, check where NetCDF lives:

```bat
dir "T:\ResMgmt\Sound\Applications\AAM\AAM_v3_dec2020\AAM-v3-Software\AAM\Bin\NCfiles"
dir "T:\ResMgmt\Sound\Applications\AAM\AAM_v3_dec2020\AAM-v3-Software\AAM\NCfiles"
```

If `NCfiles` is only under `AAM\` (not `Bin\`), either set `AAM_NC` before running:

```bat
set AAM_NC=T:\ResMgmt\Sound\Applications\AAM\AAM_v3_dec2020\AAM-v3-Software\AAM\NCfiles
```

or create a junction next to the exe:

```bat
mklink /J "T:\...\AAM\Bin\NCfiles" "T:\...\AAM\NCfiles"
```

**`Permission denied` on `Output_Data\aam\active_space.log`:** Usually a network-share write/lock issue when several omni workers log at once. For a one-off rerun, try `set AAM_PARALLEL_N=1`.

**All AAM batches empty POI:** Delete `Input_Data\aam\terrain\` under the site and rerun (stale ELV cache). Inspect `Output_Data\aam\runs\*_r000\scenario.txt`.

Send both JSON files and the last ~80 console lines.

Expected, not a crash: `[aam-filter] filtered N/M below AAM terrain`, `[aam-predict] skipped …`. Fortran stacks: `Output_Data\aam\runs\<job>\aam_stderr.txt`.
