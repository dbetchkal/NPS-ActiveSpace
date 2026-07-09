# TODO: update with pyproject.toml-based install (can see the tests.yml 
# github action + readme for details)
# NPS-ActiveSpace on Linux/Mac: Debian + Wine (for the Windows NMSim binary) +
# Python geospatial stack (GDAL/rasterio/geopandas/...).
#
# Built for linux/amd64; on Apple Silicon it runs under Docker Desktop's Rosetta
# emulation. Linux Wine has real 32-bit (i386) support, which the 32-bit
# Nord2000batch.exe (NMSim) needs.
#
# The NMSim runtime (Nord2000batch.exe + DLLs + RND/ + Sources/) is NOT baked in;
# it is bind-mounted at /opt/nmsim at run time (see docker/run_activespace.sh).
#
# The current version of NPS ActiveSpace targets Python 3.12. 
# This Dockerfile should be updated if we move to a newer version of Python
FROM --platform=linux/amd64 python:3.12-slim-bookworm

# ---- System packages: Wine (i386 for 32-bit NMSim) + GDAL/PROJ/GEOS + build tools ----
RUN dpkg --add-architecture i386 \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
      wine wine32 wine64 winbind \
      cabextract ca-certificates xvfb xauth \
      procps findutils git \
      build-essential \
      gdal-bin libgdal-dev \
      proj-bin libproj-dev \
      libgeos-dev \
      libpq-dev \
 && rm -rf /var/lib/apt/lists/*

# ---- Python venv with the NPS-ActiveSpace dependencies ----
# numpy is installed first so the GDAL bindings build against the pinned numpy ABI.
# GDAL (osgeo) is pinned to the system libgdal version reported by gdal-config.
ENV VENV=/opt/venv
RUN python3 -m venv "$VENV"
ENV PATH="$VENV/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

COPY requirements.container.txt /tmp/requirements.container.txt
RUN pip install --no-cache-dir "numpy==2.3.1" \
 && pip install --no-cache-dir --no-build-isolation "GDAL==$(gdal-config --version)" \
 && pip install --no-cache-dir -r /tmp/requirements.container.txt

# ---- Wine shim so `project.nmsim` can be a plain command that calls the .exe ----
COPY docker/nord2000 /usr/local/bin/nord2000
RUN chmod +x /usr/local/bin/nord2000

# ---- Non-root user (Wine dislikes running as root) ----
RUN useradd -m -s /bin/bash wineuser
USER wineuser
ENV HOME=/home/wineuser \
    WINEPREFIX=/home/wineuser/.wine \
    WINEARCH=win64 \
    WINEDEBUG=-all \
    WINEDLLOVERRIDES="mscoree,mshtml=" \
    NMSIM_HOME=/opt/nmsim \
    PYTHONPATH=/work

# TODO: understand the below/what is Mono/Gecko
# Initialize the Wine prefix at build time so first run is fast (no Mono/Gecko prompts).
RUN wineboot --init && wineserver -w || true

WORKDIR /work
CMD ["/bin/bash"]
