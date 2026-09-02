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
- `[project] aam` — `AAM_3.0.0.exe` (`NCfiles\` next to the exe, or under the parent `AAM\` folder; override with env `AAM_NC`)
- `[project] nmsim` — `Nord2000batch.exe`
- `[data] nvspl_archive` — NVSPL archive (ambience is computed from NVSPL; no pickle yet)

Generation uses `Input_Data\01_ELEVATION\elevation_m_nad83_utm*.tif` (and `.flt`/`.hdr`) from `project_setup`. `[data] dem` is only the parent raster for that setup step. Do not point generate at full-state `AKR_DEM.TIF`.

## Same inputs on both runs

Keep a single `DENATRLA2025_saved_annotations.geojson` in the site directory root. Do not mix Mac example tracks with the Windows set; regenerate NMSim here.

## Mac-parity smoke (platform compare)

The full `docs\tmp_windows_aam_1500m.bat` is density **48**, omni **0–2**. That is not what finished on Mac Docker. For a same-flags platform compare, use:

```bat
docs\tmp_windows_aam_macparity.bat
```

That is 1500 m, density **10**, one omni (`+000`), heading 0, AAM then NMSim, live NVSPL. Same generate flags as:

```bash
docs/tmp_docker_aam_macparity.sh
```

Both write `docker/denatrla_1500m_macparity_{aam,nmsim,metrics}.json`. The metrics file has wall time, peak process-tree RSS, CPU seconds, and 5 s samples.

## Run

From the repo root, venv optional (the script activates `.venv` if present):

```bat
docs\tmp_windows_aam_1500m.bat
docs\tmp_windows_aam_1500m.bat smoke
```

That runs AAM then NMSim. Omits `-a`, so both use default NVSPL ambience. JSON results: `denatrla_1500m_aam.json` and `denatrla_1500m_nmsim.json` at the repo root.

Outputs: `Output_Data\aam\ACTIVESPACES\DENATRLA2025_1500m\` and `Output_Data\nmsim\ACTIVESPACES\`. Do not run `fit_3d_active_space` for a single altitude.

## Compare

- **F1 / gain:** the two JSON files (`F1`, `1/3rd Octave Gain (F1)`).
- **Shape:** `python nps_active_space\scripts\viz.py DENATRLA2025 -e windows -s -a --compare -g <gain>`  
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

**`Permission denied` on `Output_Data\aam\active_space.log`:** Usually a network-share write/lock issue when several omni workers log at once. Confirm you can create/edit files under `Output_Data\aam\`. For a one-off rerun, try `set AAM_PARALLEL_N=1`. Recent code drops log lines instead of crashing if the share rejects append.

Send both JSON files and the last ~80 console lines.

Expected, not a crash: `[aam-filter] filtered N/M below AAM terrain`, `[aam-predict] isolating …` / `skipped 1 point`. Fortran stacks: `Output_Data\aam\runs\<job>\aam_stderr.txt`.
