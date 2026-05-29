# Configuration

Environment **TOML** files hold paths and credentials for an active space environment (typically one park or shared data layout). Every pipeline script takes **`-e <name>`** to load `{name}.toml`.

Bundled keys and defaults live in [`template.toml`](template.toml). Schema: [`models.py`](models.py). Load/init: [`io.py`](io.py), [`paths.py`](paths.py).

## Quick start

After `pip install -e .` from the repo root:

```bash
nps-init-config -e <environment_name>
```

That copies the template to `~/.nps_active_space/<environment_name>.toml` (or to `NPS_ACTIVE_SPACE_CONFIG_DIR` if set). Edit paths and credentials, then run scripts with `-e <environment_name>`. Repeat with other names for additional config environments.

**Blank strings** in TOML are treated as unset (`None` after load).

Example config folder:

```text
NPS_ActiveSpace_config/
├── DENA.toml
├── OLYM.toml
└── GRSM.toml
```

## Config directory

| Where configs live | How tools find them |
| ------------------ | ------------------- |
| **Custom folder** (recommended if you'd like to share config files on a network drive): local path or **network share** | Set `NPS_ACTIVE_SPACE_CONFIG_DIR` once per user or machine. Only that folder is searched. Run `nps-init-config -e <environment_name>` after setting the variable (or pass `--config-dir`). |
| **Default per user**: `~/.nps_active_space` | Do not set the env var. On Windows: `%USERPROFILE%\.nps_active_space`. `nps-init-config` creates this folder if needed. |

**Windows (cmd)**

```bat
set NPS_ACTIVE_SPACE_CONFIG_DIR=\\fileshare\NPS-ActiveSpace-config
nps-ground-truthing -e DENA -u DENA -s TRLA -y 2018
```

Persist across sessions: `setx NPS_ACTIVE_SPACE_CONFIG_DIR "..."` (then open a new terminal).

**macOS/Linux:**

```bash
export NPS_ACTIVE_SPACE_CONFIG_DIR=/path/to/configs
nps-ground-truthing -e DENA -u DENA -s TRLA -y 2018
```

Install overview: [repository README — Step 5](../../../README.md#step-5-configuration). Per-script flags: [scripts README](../scripts/README.md).

## TOML fields

Paths are filesystem locations unless noted. Most scripts require **`[project] dir`**.

### `[database]`

PostgreSQL **overflights** database for GPS flight tracks (not a filesystem path).

| Key | Required for |
| --- | ------------ |
| `name`, `username`, `password`, `port`, `host` | `nps-ground-truthing` (`-t GPS`); `nps-audible-transits` (`-t GPS`) |

### `[data]`

Paths usually live **outside** `[project] dir` (shared archives and rasters).

| Key | Description | Required for |
| --- | ----------- | ------------ |
| `nvspl_archive` | Root of the iyore NVSPL acoustic monitoring archive (hourly `.txt`) | `nps-ground-truthing`; `nps-generate-active-space` (`-a nvspl`); `nps-generate-3d-active-space` (`-a nvspl`); `nps-acoustic-metrics`; `nps-geographic-metrics` |
| `adsb` | Folder of ADS-B track exports (`.tsv`) | `nps-ground-truthing` (`-t ADSB`); `nps-audible-transits` (`-t ADSB`); `nps-acoustic-metrics` (`-t ADSB`); `nps-geographic-metrics` (`-t ADSB`) |
| `dem` | Optional GeoTIFF DEM | `nps-generate-active-space`; `nps-generate-active-space-mesh` (site elevation often comes from `Input_Data/01_ELEVATION` via `load_DEM` in other scripts) |
| `mennitt` | Optional statewide ambient-noise GeoTIFF (Mennitt L50) | `nps-generate-active-space` (`-a mennitt`); `nps-generate-active-space-mesh` |
| `site_metadata` | Reserved | Not read by current scripts |

### `[project]`

| Key | Description | Required for |
| --- | ----------- | ------------ |
| `dir` | Project workspace root; site folders `{unit}{site}` (e.g. `DENATRLA`) live under this root | `nps-ground-truthing`; `nps-plot-altitudes`; `nps-generate-active-space`; `nps-generate-active-space-batch`; `nps-generate-3d-active-space`; `nps-generate-active-space-mesh`; `nps-fit-3d-active-space`; `nps-audible-transits`; `nps-acoustic-metrics`; `nps-geographic-metrics`; `nps-check-study-duration`; `nps-viz` |
| `nmsim_binary` | Nord2000batch (NMSIM) executable | `nps-generate-active-space`; `nps-generate-active-space-mesh` |
| `FAA_Releasable_db` | FAA MASTER.txt releasable aircraft registry | `nps-ground-truthing` (`-t GPS` or `-t ADSB`); `nps-audible-transits` (`-t ADSB`) |
| `FAA_type_corrections` | Optional JSON: Mode S / ICAO hex → aircraft type | Optional for `nps-ground-truthing` (`-t GPS` or `-t ADSB`); `nps-audible-transits` (`-t ADSB`) |

SPLAT and other post-processing tools are separate; see [scripts README](../scripts/README.md) for acoustic metrics.

## Data and File Directory Setup

Config paths (`nvspl_archive`, `adsb`, `[database]`, etc.) and how they relate to on-disk
layout are documented in [TOML fields](#toml-fields) above. The following describes the
**NMSIM site folder tree** under `[project] dir`.

### NMSIM project directory

The config file's `[project] dir` entry points at a *project directory*, which may contain many listening locations (sites). Each listening deployment is a subfolder, conventionally named following the format`{unit}{site}` (e.g. `DENATRLA`). The toolkit expects that structure as follows:

```bash
project_directory/
├── UNITSITE_A/
│   ├── UNITSITE_A_study_area.shp
│   ├── Input_Data/
│   │   ├── 01_ELEVATION/
│   │   └── 02_IMPEDANCE/
│   │   └── 05_SITES/
│   └── Output_Data/
│       ├── ASCII/
│       ├── AUDIBLE_TRANSITS/
│       ├── IMAGES/
│       ├── SITE/
│       └── TIG_TIS/
├── UNITSITE_B/
│   ├── UNITSITE_B_study_area.shp
│   └── ...
└── ...
```

As an observer-based audibility model, each `nps_active_space` site directory (above `UNITSITE_A, UNITSITE_B, ...`) corresponds with a physical **listening location** on Earth. The files composing each site directory amount to a geographic model of a listener with no audible sound source present—only the quiescent surrounding land surface. Such a geographic model is a required input for every possible `nps_active_space` configuration scenario. 

#### Within each site folder (`UNITSITE*/`)

The overarching input of `nps_active_space` is a study area. It is a required input. It must live in the site folder ({project dir}/{unit}{site}/) and be named like UNITSITE_study_area.shp (ESRI shapefile), where UNITSITE is the site folder name (e.g. DENATRLA for unit DENA and site TRLA). It is recommended that study area geometries are saved using `NMSIM`'s native coordinate reference system (crs), NAD83 GCS North American (EPSG:4269).

After creation, ground-truthing annotation files (`{unit}{site}{year}*saved_annotations*.geojson`) and clock drift correction files (`{unit}{site}{year}_clock_drift_{GPS|ADSB|AIS}.csv` from [`utils/clock_drift.py`](../utils/clock_drift.py)) belong in the same site folder (`{project dir}/{unit}{site}/`). Most scripts discover them there by default; the ground-truthing GUI can save elsewhere if you pass `--annotation-file` explicitly.

#### Within `01_ELEVATION`:

Terrain within the geographic model is represented by a portion of [the National Elevation Dataset (Gesch, 2002)](https://apps.nationalmap.gov/downloader/). It is a required input. It must be formatted as ESRI grid-float (*.flt*) and stored in the `01_ELEVATION` subdirectory. 

#### Within `02_IMPEDANCE`:

Landcover-related variability in acoustic impedance along the propagation path is represented by a portion of [the National Landcover Dataset (NLCD; Yang , 2018)](https://apps.nationalmap.gov/downloader/). It is an optional input. Like elevation data, flow resistivity data should also be formatted as ESRI grid-float (*.flt*) for use with `NMSIM`. Each landcover category may be coarsely mapped to a corresponding value of flow resistivity ($\frac{Pa \ s}{m^2} = \frac{kg}{m^3 \ s}$) following the guidelines of Table I (<a href="https://doi.org/10.1121/10.0030300">adapted</a> from Ikelheimer and Plotkin, 2005b; Embleton, 1996b; Plovsing and Kragh, 2001).

<table>
  <caption>
    <strong>Table I.</strong> NLCD class mapped to flow resistivity value.
  </caption>
  <thead>
    <tr>
      <th>NLCD Class</th>
      <th>Description</th>
      <th>Flow Resistivity ($\frac{Pa \ s}{m^2}$)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>11</td>
      <td>Open water</td>
      <td>100&nbsp;000</td>
    </tr>
    <tr>
      <td>21</td>
      <td>Developed, open space</td>
      <td>200</td>
    </tr>
    <tr>
      <td>22</td>
      <td>Developed, low intensity</td>
      <td>300</td>
    </tr>
    <tr>
      <td>23</td>
      <td>Developed, medium intensity</td>
      <td>5&nbsp;000</td>
    </tr>
    <tr>
      <td>24</td>
      <td>Developed, high intensity</td>
      <td>10&nbsp;000</td>
    </tr>
    <tr>
      <td>31</td>
      <td>Barren (e.g., bedrock, talus)</td>
      <td>100&nbsp;000</td>
    </tr>
    <tr>
      <td>41</td>
      <td>Deciduous forest</td>
      <td>70</td>
    </tr>
    <tr>
      <td>42</td>
      <td>Evergreen forest</td>
      <td>70</td>
    </tr>
    <tr>
      <td>43</td>
      <td>Mixed forest</td>
      <td>70</td>
    </tr>
    <tr>
      <td>51</td>
      <td>Dwarf scrub</td>
      <td>200</td>
    </tr>
    <tr>
      <td>52</td>
      <td>Shrub/scrub</td>
      <td>200</td>
    </tr>
    <tr>
      <td>90</td>
      <td>Woody wetlands</td>
      <td>100&nbsp;000</td>
    </tr>
  </tbody>
</table>

#### Within `05_SITES`:

Any three-dimensional coordinate above Earth's surface may be used as a listener location. The listener location should be formatted as a `NMSIM` (*.sit*) file and stored in the `05_SITES` subdirectory.
