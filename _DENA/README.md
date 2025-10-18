# \_DENA

This directory contains scripts written specifically for Denali National Park and Preserve as it is where the
project was developed. None of the code contained in this directory was written with the intention for it to be imported
and used by others. However, these scripts can serve as examples of how the NPS-ActiveSpace modules can be used.

## Installation

1. Clone the NPS-ActiveSpace repository.

```bash
git clone https://github.com/dbetchkal/NPS-ActiveSpace.git
```

2. Install project dependencies.

Three dependencies, GDAL, Fiona, and rasterio need to be installed from `.whl` files. Please
[download](https://www.lfd.uci.edu/~gohlke/pythonlibs/) a GDAL, Fiona, and rastertio binary that matches the version
number specified in `requirements.txt` and that matches the python version you will be running the NPS-ActiveSpace code
with. For example `Fiona‑1.8.21‑cp310‑cp310‑win_amd64.whl` is Fiona version 1.8.21 for python 3.10. Save all three files in
the same location.

Then, run the following commands:

```bash
$ python -m pip install --upgrade pip
$ pip install --find-links </path/to/binaries> -r requirements.txt
```

3. Create config files.

All scripts in this directory require a configuration file. Please copy the template config file, fill
in the values required for the script(s) you will be running, and save it to the config directory as `<environment name>.config`.
For example, a DENA production configuration file might be named `dena_production.config` while a HAVO production
configuration file might be named `havo_production.config` and have a different value for where the DEM file is stored
than `dena_production.config`

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
dem = Absolute path to the DEM tif file to use for active space generation. Value required for generate_active_space.py and generate_active_space_mesh.py
mennitt = Absolute path to the mennitt ambience tif. Value required for generate_active_space.py and generate_active_space_mesh.py

[project]
dir = Absolute path to the directory where all NPS-ActiveSpace files are stored. Required for all scripts.
nmsim = Absolute path to the NMSIM Nord2000batch.exe file. Value required for generate_active_space.py and generate_active_space_mesh.py
FAA_Releasable_db = Absolute path to the FAA MASTER.txt database file downloaded from the [FAA website](https://www.faa.gov/licenses_certificates/aircraft_certification/aircraft_registry/releasable_aircraft_download). Required for run_audible_transits.py
FAA_type_corrections = Absolute path to a json file for correcting aircraft types in the FAA database. Keys are ICAO addresses, values are correct aircraft type. Required for run_ground_truthing.py and run_audible_transits.py
```

## Directories

`config/`: All `.config` should be placed here.

`resource/`: This directory contains helper functions that are used by multiple scripts.

`scripts/`: Home to scripts for the various sound management plan creation steps. See details of each script below.

## Scripts

### Ground Truthing

This script is used to launch the ground truthing application to annotate the audibility of sound source tracks.

| command-line arg        | description                                                                                                                                      |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `-e`, `--environment`   | **required.**<br/>The configuration environment to use. _Ex_: To use `production.config` pass `-e production`                                    |
| `-u`, `--unit`          | **required.**<br/>The 4 letter NPS unit code. _Ex_: Denali = DENA                                                                                |
| `-s`, `--site`          | **required.**<br/>The 4 letter site code. _Ex_: Cathedral = CATH                                                                                 |
| `-y`, `--year`          | **required.**<br/>The deployment year, YYYY. _Ex_: 2018                                                                                          |
| `-t`, `--database-type` | **_default GPS -> {GPS, ADSB, AIS}_**<br/>Which track source to use. Paths and login credentials for all source types are stored in config files |

Example executions:

```bash
$ python -u -W ignore _DENA/scripts/run_ground_truthing.py -e production -u DENA -s MOOS -y 2018
```

```bash
$ python -u -W ignore _DENA/scripts/run_ground_truthing.py -e production -u DENA -s TRLA -y 2018 -t ADSB
```

### Generate Active Space

This script is used to generate active spaces for a single site for a variety of omni sources to determine which
omni source produces the active space that most closely matches the ground truthed tracks.

Please note that the Precision-Recall plot that is shown at the end of a run is automatically saved.

| command-line arg      | description                                                                                                                                                                                                                                                             |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-e`, `--environment` | **required.**<br/>The configuration environment to use. _Ex_: To use `production.config` pass `-e production`                                                                                                                                                           |
| `-u`, `--unit`        | **required.**<br/>The 4 letter NPS unit code. _Ex_: Denali = DENA                                                                                                                                                                                                       |
| `-s`, `--site`        | **required.**<br/>The 4 letter site code. _Ex_: Cathedral = CATH                                                                                                                                                                                                        |
| `-y`, `--year`        | **required.**<br/>The deployment year, YYYY. _Ex_: 2018                                                                                                                                                                                                                 |
| `-a`, `--ambience`    | **_default nvspl -> {nvspl, mennitt, or .pkl file path}_**<br/>The ambience type to use when running NMSIM.                                                                                                                                                             |
| `--headings`          | **_default [0, 120, 240]_**<br/>A list of the active space headings that should be dissolved together to make the final active space. _Ex_: `--headings 0, 90, 180, 270`                                                                                                |
| `--omni-min`          | **_default -10.0_**<br/>The lowest gain to generate an active space for. Active spaces will be generated for all gains between `--omni-min` and `--omni-max`.                                                                                                           |
| `--omni-max`          | **_default 40.0_**<br/>The highest gain to generate an active space for. Active spaces will be generated for all gains between `--omni-min` and `--omni-max`.                                                                                                           |
| `-l`, `--altitude`    | Use this flag to generate the active spaces at a particular altitude (in meters). _Ex_: `-l 1524` generates active spaces at 1524 meters or 5000 feet.<br/>If not passed, the average altitude of the valid, audible ground-truthed tracks will be calculated and used. |
| `-b`, `--beta`        | **_default 1.0_**<br/>the beta value to use when calculating the f-beta for each active space.<br/>https://en.wikipedia.org/wiki/F-score#F%CE%B2_score)                                                                                                                 |
| `--cleanup`           | If this flag is added, all intermediary control and batch files will be deleted upon script completion.                                                                                                                                                                 |
| `--annotation-file`   | If provided, basename of GEOJSON annotations file to use instead of the default. File should be in the site directory.                                                                                                                                                  |

Example executions:

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space.py -e production -u DENA -s MOOS -y 2018 --cleanup
```

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space.py -e production -u DENA -s TRLA -y 2017  -a mennitt --headings 0 --omni-min -5 --omni-max 10.5 -l 3658 -b .5
```

If you would like the command output to be shown in the console and saved to a text file add the following after your command:

```bash
<command> | Tee-Object -FilePath "C:\Path\To\Output.txt"
```

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space.py -e production -u DENA -s MOOS -y 2018 --cleanup | Tee-Object -FilePath "C:\Path\To\active_space_output_DENAMOOS2018.txt"
```

### Generate Active Space Mesh

| command-line arg      | description                                                                                                                                                                   |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-e`, `--environment` | **required.**<br/>The configuration environment to use. _Ex_: To use `production.config` pass `-e production`                                                                 |
| `-n`, `--name`        | **required.**<br/>Name of the directory where intermediary and output files will be stored. _Ex_: `-n DENAFULL`                                                               |
| `-s`, `--study-area`  | **required.**<br/>Absolute path to the shapefile of the study area. _Ex_: `-s C:/Users/yourname/Desktop/DENA.shp`                                                             |
| `--headings`          | **_default [0, 120, 240]_**<br/>A list of the active space headings that should be dissolved together to make the final active space. _Ex_: `--headings 0, 90, 180, 270`      |
| `--omni-source`       | **_default 0_**<br/>Gain to generate the mesh with.                                                                                                                           |
| `--mesh-spacing`      | **_default 1_**<br/>How far apart, in km, mesh square centroids should be.                                                                                                    |
| `--mesh-size`         | **_default 25_**<br/>How large, in km, each mesh square should be. Mesh squares will be mesh-size x mesh-size.                                                                |
| `-l`, `--altitude`    | **_default 3658_**<br/>Use this flag to generate the active spaces at a particular altitude (in meters). _Ex_: `-l 1524` generates active spaces at 1524 meters or 5000 feet. |
| `--cleanup`           | If this flag is added, all intermediary control and batch files will be deleted upon script completion.                                                                       |

Example executions:

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space_mesh.py -e production -n DENAFULL -s C:/Users/yourname/Desktop/DENA.shp --cleanup
```

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space_mesh.py -e production -n DENAFULL -s C:/Users/yourname/Desktop/DENA.shp --headings 0 180 --omni-source -12.5 --mesh-spacing 10 --mesh-size 20 -l 1524
```

### Audible Transits

This script is used to estimate the geographic intersection of a set of tracks with an active space. It simplifies and enriches the geographic data with information necessary to produce an audiblity time series.

| command-line arg           | description                                                                                                                                      |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `-e`, `--environment`      | **required.**<br/>The configuration environment to use. _Ex_: To use `production.config` pass `-e production`                                    |
| `-u`, `--unit`             | **required.**<br/>The 4 letter NPS unit code. _Ex_: Denali = DENA                                                                                |
| `-s`, `--site`             | **required.**<br/>The 4 letter site code. _Ex_: Cathedral = CATH                                                                                 |
| `-y`, `--year`             | **required.**<br/>Which year's active space to use, YYYY. _Ex_: 2018                                                                             |
| `-g`, `--gain`             | **required.**<br/>The signed gain of the optimal active space fit, float. _Ex._: -20.5                                                           |
| `-t0`, `--begintracks`     | **required.**<br/>Date to begin parsing the position record, YYYY-MM-DD. _Ex._: 2018-01-01                                                       |
| `-tf`, `--endtracks`       | **required.**<br/>Date to stop parsing the position record, YYYY-MM-DD. _Ex._: 2018-06-05                                                        |
| `-t`, `--track-source`     | **_default GPS -> {GPS, ADSB, AIS}_**<br/>Which track source to use. Paths and login credentials for all source types are stored in config files |
| `-o`, `--output`           | **default ""**<br/>Output directory. Directory to store output files. Defaults to [project directory]/[unit][site]/Output_Data                   |
| `-garb`, `--exportgarbage` | **default 0 -> {0, 1}**<br/>Whether to export garbage tracks (1) or not (0). _Ex._: 1                                                            |
| `-v`, `--verbose`          | <br/>If provided, prints additional details to the console during processing steps.                                                              |

Example executions:

```bash
$ python _DENA/scripts/run_audible_transits.py -e production -u DENA -s FANG -y 2018 -g -1.5 -t0 2018-01-01 -tf 2024-08-20
```

```bash
$ python _DENA/scripts/run_audible_transits.py -e production -u DENA -s FANG -y 2018 -g -1.5 -t0 2018-01-01 -tf 2024-08-20 -t ADSB
```
