[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18111579.svg)](https://doi.org/10.5281/zenodo.18111579)

# NPS-ActiveSpace

An **_active space_** is a well-known sensory concept from bioacoustics ([Marten and Marler 1977](https://www.jstor.org/stable/pdf/4599136.pdf), [Gabriele et al. 2018](https://www.frontiersin.org/articles/10.3389/fmars.2018.00270/full)). It represents a geographic volume whose radii correspond to the limit of audibility for a specific signal in each direction. In other words, an active space provides an answer to the question, _"how far can you hear a certain sound source from a specific location on the Earth's surface?"_

This repository is designed to estimate active spaces for motorized noise sources transiting the U.S. National Park System. Aircraft are powerful noise sources audible over vast areas. Thus [considerable NPS management efforts have focused on protecting natural quietude from aviation noise intrusions](https://www.nps.gov/subjects/sound/overflights.htm). For coastal parks, vessels are similarly powerful noise sources of concern. For both transportation modalities `NPS-ActiveSpace` provides meaningful, quantitative spatial guides for noise mitigation and subsequent monitoring.

## Example

Consider an example active space, below. It was computed using data from a long term acoustic monitoring site in Denali National Park, DENAUWBT Upper West Branch Toklat ([Withers 2012](https://irma.nps.gov/DataStore/Reference/Profile/2184396)). The bold black polygon delineates an active space estimate for flights at 3000 meters altitude. Points interior to the polygon are predicted to be audible, those exterior, inaudible. <br>

Superposed over the polygon are colored flight track polylines. `NPS-ActiveSpace` includes an application that leverages the acoustic record to ground-truth audibility of co-variate vehicle tracks from GPS databases. Ground-truthing is used to "tune" an active space to the appropriate geographic extent via mathematical optimization.<br>

<br>
<img src="https://github.com/dbetchkal/NPS-ActiveSpace/blob/main/nps_active_space/img/NPS-ActiveSpace_example.png" alt="active space polygon example" width="200">

### Synthesis-oriented Design

At its heart, `nps_active_space` is organized around the idea of a **scientific geo-synthesis**:
a structured combination of heterogeneous spatial observations into a single, interpretable
geometric object.

Rather than treating sound propagation, audibility, and vehicle movement as
independent modeling problems, the toolkit treats them as *causally linked
components* of an observer-centered system:

- Vehicle trajectories represent potential causes
- Acoustic records represent observed effects
- Audibility is established through empirical association
- Geometry is used to reconcile these relationships in space

The resulting active space is not a purely predictive construct, nor a purely
descriptive one. It is a **mensurated estimate**: a spatial enclosure that is
consistent with observed audibility under specified environmental conditions.

## Installation

### Step 1: Clone the NPS-ActiveSpace repository.

```
git clone https://github.com/dbetchkal/NPS-ActiveSpace.git
cd NPS-ActiveSpace
```

### Step 2: Set Up a Virtual Environment.

Install Python, either via Anaconda/Miniconda, or directly. The repository has been tested with Python version 3.12.

You can use a Conda environment if you want, but all installation is managed by pip.

With Conda:

```
conda create --name active python=3.12.12
conda activate active
```

With venv in a Git Bash terminal:

```
python -m venv .venv
source .venv/bin/activate
```

With venv in a Windows Command Prompt terminal:

```
python -m venv .venv
source .venv\Scripts\activate.bat
```

### Step 3: Install Dependencies

Make sure you are inside your virtual environment, then:

```
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Historical note:
The GDAL dependency comes from a `.whl` file published [here](https://github.com/cgohlke/geospatial-wheels/releases). If the Python version is updated, the GDAL wheel URL in `requirements.txt` may need to be changed to reflect the updated version. For example, `gdal-3.11.1-cp312-cp312-win_amd64.whl` is GDAL version 3.11.1 for Python 3.12.

### Step 4: Install NPS-ActiveSpace

From the repository's root directory, inside the virtual environment:

```
pip install -e .
```

Try importing a python module to make sure this install worked, e.g. in a python file:

```
from nps_active_space.active_space import ActiveSpaceGenerator
```

### Step 5: Create Config File

All scripts require a configuration file saved in the config directory `nps_active_space/config`. Please copy the template config file, fill in the values required for the script(s) you will be running, and save it to the config directory as `<environment name>.config`. For example, a configuration file for Denali National Park and Preserve might be named `DENA.config` while a configuration file for Hawaii Volcanoes National Park might be named `HAVO.config` and have a different value for where the DEM file is stored than `DENA.config`

Currently, the template config file has the following data:

```text
[database:overflights] - Values required if pulling tracks from the database in run_ground_truthing.py or run_audible_transits.py
name = Database name.
username = Database credentials username.
password = Database credentials password.
port = Database port.
host = Database host.

[data]
site_metadata = Absolute path to the the file containing site metadata. Value required for all run_ground_truthing.py and generate_active_space.py
nvspl_archive = Absolute path to the directory where all NVSPL sound data is stored. Value required for all run_ground_truthing.py and generate_active_space.py
adsb = Absolute path to the directory where ADSB track data is stored.  Value required if pulling ADSB tracks in run_ground_truthing.py or run_audible_transits.py
dem = Absolute path to the source DEM GeoTIFF for project_setup and active space generation. Required for project_setup.py, generate_active_space.py, and generate_active_space_mesh.py
dem_elevation_units = Vertical units of the source DEM raster: ``feet`` or ``meters``. Defaults to ``feet``. Used by project_setup.py when clipping/reprojecting elevation to a site location.
mennitt = Absolute path to the mennitt ambience tif. Value required for generate_active_space.py and generate_active_space_mesh.py

[project]
dir = Absolute path to the directory where all NPS-ActiveSpace files are stored. Required for all scripts.
nmsim = Absolute path to the NMSIM Nord2000batch.exe file. Value required for generate_active_space.py and generate_active_space_mesh.py
FAA_Releasable_db = Absolute path to the FAA MASTER.txt database file downloaded from the [FAA website](https://www.faa.gov/licenses_certificates/aircraft_certification/aircraft_registry/releasable_aircraft_download). Required for run_audible_transits.py
FAA_type_corrections = Absolute path to a json file for correcting aircraft types in the FAA database. Keys are ICAO addresses, values are correct aircraft type. Required for run_ground_truthing.py and run_audible_transits.py
```

## Toolkit Architecture

```mermaid
graph LR

    %% Groups
    subgraph Scripts
        scripts[(scripts/)]
    end

    subgraph Synthesis
        active_space[(active_space/)]
        ground_truthing[(ground_truthing/)]
         validation[(validation/)]
    end

    subgraph Foundations
        config[(config/)]
        data[(data/)]
        utils[(utils/)]
    end

    %% Flow: scripts → core → analysis
    config --> scripts
    config --> ground_truthing
    data --> scripts
    utils --> scripts
    scripts --> active_space
    scripts --> ground_truthing
    scripts --> validation

    %% Foundations feed core
    config --> active_space
    data --> active_space
    utils --> active_space

    %% Foundations feed analysis modules
    data --> ground_truthing
    utils --> ground_truthing

    ground_truthing --> active_space
    active_space --> validation
    utils --> validation
```

### The Synthesis

At its heart, this toolkit structures a **scientific geo-synthesis** in two parts:

1) A GUI-based audibility measurement tool (`.ground_truthing`) that allows users to spatio-temporally $\text{JOIN}$ vehicle tracks (cause) and acoustic records (effect).
2) A geoprocess to $\text{ENCLOSE}$ the listening location and a user's audibility observations within an optimal 3-dimensional active space (`.active_space`).

This enables a causal geometric calculation (mensuration) of acoustic metrics. Tools to help validate a set of syntheses are provided (`.validation`).

### Scripted Control

The entire synthesis $\rightarrow$ validation $\rightarrow$ mensuration workflow is scripted in a Command Line Interface (`.scripts`). For more information on scripted control, see the detailed [control script documentation](https://github.com/dbetchkal/NPS-ActiveSpace/tree/main/nps_active_space/scripts#scripts).

### Foundations

To enable environmental scenario planning, global repository inputs are all structured together in a configuration file (`.config`) that spans the architecture. Utilities implement both data structure models and diverse computation tasks (`.utils`) that are used throughout the architecture. 

The toolkit also requires basic sound source and weather profile data to operate. The provided (`.data`) are a widely applicable example: a fixed-wing propeller aircraft source and a standard "acoustician's atmosphere" with dry adiabatic lapse conditions. The toolkit may be run using alternative, custom sound sources or weather profiles. 

### Data and File Directory Setup

#### `NMSIM` project directory setup:

The `.config` configuration file requires a *project directory*, which may contain many listening locations (sites), which in turn can be configured flexibly to produce test many scenarios. The toolkit expects the user's project directory to be structured as follows:

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

##### Within the root of the project directory:

The overarching input of `nps_active_space` is a study area. It is a required input. It must be contained within the root project directory and named similar to `UNITSITE_study_area.shp` (ESRI shapefile), where the variable geographic prefix `UNITSITE` matches the name of the project directory. It is recommended that study area geometries are saved using `NMSIM`'s native coordinate reference system (crs), NAD83 GCS North American (EPSG:4269).

After their creation, the annotation outputs of `.ground_truthing` (.geojson) also belong in the root of the project directory. Clock drift correction files produced via [the `.utils\clock_drift.py` utility](https://github.com/dbetchkal/NPS-ActiveSpace/blob/main/nps_active_space/utils/clock_drift.py) also belong in the root of the project directory.

##### Within `01_ELEVATION`:

Terrain within the geographic model is represented by a portion of [the National Elevation Dataset (Gesch, 2002)](https://apps.nationalmap.gov/downloader/). It is a required input. It must be formatted as ESRI grid-float (*.flt*) and stored in the `01_ELEVATION` subdirectory. 

##### Within `02_IMPEDANCE`:

Landcover-related variability in acoustic impedance along the propagation path is represented by a portion of [the National Landcover Dataset (NLCD; Yang , 2018)](https://apps.nationalmap.gov/downloader/). It is an optional input. Like elevation data, flow resistivity data should also be formatted as ESRI grid-float (*.flt*) for use with `NMSIM`. Each landcover category may be coasrely mapped to a corresponding value of flow resistivity ($\frac{Pa \ s}{m^2} = \frac{kg}{m^3 \ s}$) following the guidelines of Table I (<a href="https://doi.org/10.1121/10.0030300">adapted</a> from Ikelheimer and Plotkin, 2005b; Embleton, 1996b; Plovsing and Kragh, 2001).

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

##### Within `05_SITES`:

Any three-dimensional coordinate above Earth's surface may be used as a listener location. The listener location should be formatted as a `NMSIM` (*.sit*) file and stored in the `05_SITES` subdirectory.

## Intended Use & Status

`nps_active_space` is intended for **observer-based, post hoc audibility analysis** of transportation noise using simultaneously collected, long-term acoustic records and covariate vehicle tracking data. It is designed to support scientific research, environmental assessment, and resource management applications, particularly in large, heterogeneous landscapes such as U.S. National Parks.

The toolkit is especially suited to problems where:

- A fixed listening location (or small set of locations) is monitored over time
- Sound sources are mobile (e.g., aircraft, vessels, overland vehicles)
- Spatial questions are central (e.g., *where* is a source audible from, and *for how long?*)

`nps_active_space` is **not** a real-time noise monitoring system, nor a general-purpose
sound propagation model. It is a synthesis-oriented toolkit that combines
observations, geometry, and physical modeling to estimate spatial limits of
audibility under specified environmental conditions.

### Project Status

This software is:

- Research-grade and actively used in peer-reviewed studies
- Designed for transparency, reproducibility, and extensibility
- Under continued development as methods and use cases evolve

Interfaces and workflows may change as new synthesis pathways are added. Users are encouraged to treat outputs as **scientific estimates** whose validity depends on data quality, configuration choices, and underlying assumptions.

---

## License

### Public domain

This project is in the worldwide [public domain](LICENSE.md):

> This project is in the public domain within the United States,
> and copyright and related rights in the work worldwide are waived through the
> [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
>
> All contributions to this project will be released under the CC0 dedication.
> By submitting a pull request, you are agreeing to comply with this waiver of copyright interest.

## Publications

The original technical paper (supporting the v1.0.0 Release):
> Betchkal, D.H., J.A. Beeco, S.J. Anderson, B.A. Peterson, and D. Joyce. 2023. Using Aircraft Tracking Data to Estimate the Geographic Scope of Noise Impacts from Low-Level Overflights Above Parks and Protected Areas. Journal of Environmental Management 348(15): 119201 https://doi.org/10.1016/j.jenvman.2023.119201 <br><br>

Validation on jet aircraft source types (v2.1.0 Release):
> Gurung, B., Betchkal, D.H., Beeco, J.A., Peterson, B.A., Olstad, T.A., Anderson, S., Hutchinson, S., Jackson, S. and Joyce, D., 2026. Estimating Active Space Noise Extent from Two Aircraft Weight Classes over the Great Smoky Mountains National Park. Aerospace, 13(4), p.363. https://doi.org/10.3390/aerospace13040363 <br><br>

A testing presentation (supporting the v3.0.0 Release):
> Kosinski, A. and Betchkal, D. H. 2025. Reliability of the 3D Active Space Model. U.S. National Park Service. https://irma.nps.gov/DataStore/Reference/Profile/2316524 <br>


Management papers using `NPS-ActiveSpace`:
> Peterson, B.A., Coffey, L., Gurung, B., Hutchinson, J.S., Olstad, T.A., Betchkal, D.H., Anderson, S.J., Joyce, D. and Beeco, J.A., 2025. Understanding Overflight Travel Patterns Above Parks and Protected Areas: A Systematic Review. Journal of Park and Recreation Administration. https://doi.org/10.18666/JPRA-2025-13120 <br><br>
> Ferguson, L.A., Peterson, B.A., Crump, M., Taff, B.D., Newman, P., Betchkal, D.H., Hutchinson, J.S. and Beeco, J.A., 2025. Integrating aircraft tracking, acoustic data, and surveys to evaluate park aircraft noise. Tourism Geographies, pp.1-20. https://doi.org/10.1080/14616688.2025.2591832 <br><br>

