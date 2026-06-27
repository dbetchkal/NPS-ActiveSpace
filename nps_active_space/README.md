# Installation

### Clone the NPS-ActiveSpace repository.

```
git clone https://github.com/dbetchkal/NPS-ActiveSpace.git
cd NPS-ActiveSpace
```

### Set Up a Virtual Environment.

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

### Install Dependencies

Make sure you are inside your virtual environment, then:

```
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Historical note:
The GDAL dependency comes from a `.whl` file published [here](https://github.com/cgohlke/geospatial-wheels/releases). If the Python version is updated, the GDAL wheel URL in `requirements.txt` may need to be changed to reflect the updated version. For example, `gdal-3.11.1-cp312-cp312-win_amd64.whl` is GDAL version 3.11.1 for Python 3.12.

### Install NPS-ActiveSpace

From the repository's root directory, inside the virtual environment:

```
pip install -e .
```

Try importing a python module to make sure this install worked, e.g. in a python file:

```
from nps_active_space.active_space import ActiveSpaceGenerator
```

### Create Config File

All scripts require a configuration file saved in the config directory `nps_active_space/config`. Please copy the template config file, fill in the values required for the script(s) you will be running, and save it to the config directory as `<environment name>.config`. For example, a DENA configuration file might be named `dena.config` while a HAVO configuration file might be named `havo.config` and have a different value for where the DEM file is stored than `dena.config`

Currently, the template config file has the following data:

TODO - check this / update this for new scripts. Explaining which scripts require what is a bit messy, especially with all the new scripts we added. Consider a better way to document this.

```text
[database:overflights] - Values required if pulling GPS tracks from the database in run_ground_truthing.py or run_audible_transits.py
name = Database name.
username = Database credentials username.
password = Database credentials password.
port = Database port.
host = Database host.

[data]
site_metadata = Absolute path to the the file containing site metadata. Value required for all run_ground_truthing.py and generate_active_space.py
nvspl_archive = Absolute path to the directory where all NVSPL sound data is stored. Value required for all run_ground_truthing.py and generate_active_space.py
adsb = Absolute path to the directory where ADSB track data is stored.  Value required if pulling ADSB tracks in run_ground_truthing.py or run_audible_transits.py
ais = Absolute path to the MXAK AIS archive directory (e.g. MXAK-AIS-GLBA). Value required if pulling AIS tracks in run_ground_truthing.py with -t AIS
dem = Absolute path to the DEM tif file to use for active space generation. Value required for generate_active_space.py and generate_active_space_mesh.py
mennitt = Absolute path to the mennitt ambience tif. Value required for generate_active_space.py and generate_active_space_mesh.py

[project]
dir = Absolute path to the directory where all NPS-ActiveSpace files are stored. Required for all scripts.
nmsim = Absolute path to the NMSIM Nord2000batch.exe file. Value required for generate_active_space.py and generate_active_space_mesh.py
FAA_Releasable_db = Absolute path to the FAA MASTER.txt database file downloaded from the [FAA website](https://www.faa.gov/licenses_certificates/aircraft_certification/aircraft_registry/releasable_aircraft_download). Required for run_audible_transits.py
FAA_type_corrections = Absolute path to a json file for correcting aircraft types in the FAA database. Keys are ICAO addresses, values are correct aircraft type. Required for run_ground_truthing.py and run_audible_transits.py
```

# Data and File Directory Setup

TODO

Define "project directory" as the directory containing the site folders. This is consistent with how the code defines it.

Please document:

- acoustic data (NVSPL, SPLAT)
- ADSB data
- GPS data (database setup and usage via config file)
- project directory
- site directories, including minimal manual setup steps (what directories the user needs to make vs. what's made automatically by the code)
  - mention that annotation files and clock drift correction files go in here
- NMSIM
