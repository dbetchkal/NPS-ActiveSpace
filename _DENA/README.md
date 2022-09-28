# _DENA

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
in the values, and save it to the config directory as `<environment name>.config`. For example, a production
configuration file might be named `production.config`.

## Directories

`config/`: All `.config` should be placed here.

`resource/`: This directory contains helper functions that are used by multiple scripts.

`scripts/`: Home to scripts for the various sound management plan creation steps. See details of each script below.

## Scripts

### Ground Truthing

This script is used to launch the ground truthing application to annotate the audibility of sound source tracks.

| command-line arg       | description                                                                                                                                                |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-e`, `--environment`  | **required.**<br/>The configuration environment to use. *Ex*: To use `production.config` pass `-e production`                                              |
| `-u`, `--unit`         | **required.**<br/>The 4 letter NPS unit code. *Ex*: Denali = DENA                                                                                          |
| `-s`, `--site`         | **required.**<br/>The 4 letter site code. *Ex*: Cathedral = CATH                                                                                           |
| `-y`, `--year`         | **required.**<br/>The deployment year, YYYY. *Ex*: 2018                                                                                                    |
| `-t`, `--track-source` | ***default Database -> {Database, ADSB, AIS}***<br/>Which track source to use. Paths and login credentials for all source types are stored in config files |

Example executions:

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space_mesh.py -e production -u DENA -s MOOS -y 2018
```

```bash
$ python -u -W ignore _DENA/scripts/generate_active_space_mesh.py -e production -u DENA -s TRLA -y 2018 -t ADSB
```

### Generate Active Space

| command-line arg      | description             |
|-----------------------|-------------------------|
| `-e`, `--environment` | **required.**           |
| `-u`, `--unit`        | **required.**           |
| `-s`, `--site`        | **required.**           |
| `-y`, `--year`        | **required.**           |
| `-a`, `--ambience`    | *default nvspl*         |
| `--headings`          | *default [0, 120, 240]* |
| `--omni-min`          | *default -20*           |
| `--omni-max`          | *default 30*            |
| `-l`, `--altitude`    |                         |
| `-b`, `--beta`        | *default 1.0*           |
| `--cleanup`           |        |

### Generate Active Space Mesh

| command-line arg      | description             |
|-----------------------|-------------------------|
| `-e`, `--environment` | **required.**           |
| `-n`, `--name`        | **required.**           |
| `-s`, `--study-area`  | **required.**           |
| `--headings`          | *default [0, 120, 240]* |
| `--omni-source`       | *default 0*             |
| `--mesh-spacing`      | *default 1*             |
| `--mesh-size`         | *default 25*            |
| `-l`, `--altitude`    |                         |
| `--cleanup`           |        |


```bash
$ python -u -W ignore .\_DENA/scripts/generate_active_space_mesh.py -e test -n TRLATEST -s C:/Users/azucker/Desktop/DENATRLA2021_study_area.shp --omni-min 0 --omni-max 0 --mesh-spacing 30
```
