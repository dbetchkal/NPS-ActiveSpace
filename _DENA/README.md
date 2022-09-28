# _DENA

This directory contains scripts written specifically for Denali National Park and Preserve as it is where the 
project was developed. However, these scripts can serve as examples of how the NPS-ActiveSpace modules can be used. 

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

All of the scripts in this directory require a configuration file. Please copy the template config file and fill
in the values.


## Scripts

### Ground Truthing

| command-line arg       | description        |
|------------------------|--------------------|
| `-e`, `--environment`  | **required.**      |
| `-u`, `--unit`         | **required.**      |
| `-s`, `--site`         | **required.**      |
| `-y`, `--year`         | **required.**      |
| `-t`, `--track-source` | *default Database* |


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
