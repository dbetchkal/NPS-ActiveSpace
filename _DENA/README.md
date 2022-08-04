# _DENA

## About


## Installation

1. Clone the NPS-ActiveSpace repository.
```bash
git clone https://github.com/dbetchkal/NPS-ActiveSpace.git
```

2. Install project dependencies.

**For Windows users:**

https://www.lfd.uci.edu/~gohlke/pythonlibs/

Two dependencies, GDAL and Fiona, need to be installed from `.whl` files. Please download a GDAL
and Fiona binary that matches the version number specified in `pyproject.toml` and that matches
the python version you will be running the NPS-ActiveSpace code with. For example `Fiona‑1.8.21‑cp310‑cp310‑win32.whl` is Fiona version 1.8.21 for python 3.10.
Save both files in the same location.

Then, run the following commands:

```bash
$ python -m pip install --upgrade pip
$ pip install --find-links==</path/to/binaries> -r requirements.txt
```

**For all other users:**

```bash
$ python -m pip install --upgrade pip
$ pip install -r requirements.txt
```
