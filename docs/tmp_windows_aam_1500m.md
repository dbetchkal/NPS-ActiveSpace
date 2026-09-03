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

**Lightweight** cross-platform job — density **10**, omni **+000** only (not the standard d48 mesh):

```bat
docs\tmp_windows_aam_macparity.bat
```

Same flags as Mac Docker `docs/tmp_docker_aam_macparity.sh`. Writes `docker/denatrla_1500m_macparity_{aam,nmsim,metrics}.json`.

## Standard single-layer 1500 m

**Standard** mesh/gain — density **48** (default), omni **0–2**:

```bat
docs\tmp_windows_aam_1500m.bat
```

Quick smoke variant (d10, omni 0): `docs\tmp_windows_aam_1500m.bat smoke`.

## 3D pipeline (AAM + NMSim)

Full 3D product path (900–2100 m layers, density **48**, omni **0–2**, batch + `fit_3d_active_space.py`):

```bat
docs\tmp_windows_aam_3d.bat
docs\tmp_windows_aam_3d.bat aam
```

See [tmp_windows_aam_3d.md](tmp_windows_aam_3d.md). Mac Docker: `docs/tmp_docker_aam_3d.sh`.

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

**`Same Profile data found in multiple spheres`:** AAM was pointed at the whole site ``Input_Data\aam\NCfiles`` folder (multiple ``OMNI_*.nc`` from omni 0–2). Current code stages one omni per run under ``Output_Data\aam\runs\{job}\NCfiles\``. Pull the branch fix or rerun after upgrade; older builds need ``--omni-max 0`` as a workaround.

**All AAM batches empty POI:** Often `Bin\NCfiles` is an empty stub without `OMNI_200.nc` (the read-only template). ActiveSpace now generates per-gain sources under `Input_Data\aam\NCfiles\` from sibling `.avg` files; the vendor tree only supplies the template. Check:

```bat
dir "T:\...\AAM\Bin\NCfiles\OMNI_200.nc"
dir "T:\...\AAM\NCfiles\OMNI_200.nc"
set AAM_NC=T:\...\AAM\NCfiles
dir "%SITE%\Input_Data\aam\NCfiles\OMNI_000.nc"
```

If the template exists but POI is still empty, open `Output_Data\aam\runs\*_r000\scenario.txt` and search `NETCDF`, `OMNI`, `Below ground`. Delete `Input_Data\aam\terrain\` if ELV looks stale.

Send both JSON files and the last ~80 console lines.

Expected, not a crash: `[aam-filter] filtered N/M below AAM terrain`, `[aam-predict] skipped …`. Fortran stacks: `Output_Data\aam\runs\<job>\aam_stderr.txt`.
