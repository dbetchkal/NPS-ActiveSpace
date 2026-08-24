# NPS-ActiveSpace on Linux/Mac: Debian + Wine (for Windows acoustic-model binaries) +
# Python geospatial stack (GDAL/rasterio/geopandas/...).
#
# Built for linux/amd64; on Apple Silicon it runs under Docker Desktop's Rosetta
# emulation. Linux Wine has real 32-bit (i386) support, which the 32-bit
# Nord2000batch.exe (NMSim) needs.
#
# Acoustic-model runtimes (NMSim, AAM) are NOT baked in; they are bind-mounted at
# run time (see docker/run_activespace.sh and ACOUSTIC_MODEL).
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
ENV VENV=/opt/venv
RUN python3 -m venv "$VENV"
ENV PATH="$VENV/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

WORKDIR /repo
COPY pyproject.toml README.md LICENSE.md ./
COPY nps_active_space/ ./nps_active_space/
COPY vendor/.docker-build/aam_translator /tmp/aam_translator
RUN pip install --no-cache-dir "GDAL==$(gdal-config --version)" \
 && pip install --no-cache-dir /tmp/aam_translator \
 && pip install --no-cache-dir -e .

# ---- Wine shims (binaries are bind-mounted at run time) ----
COPY docker/nord2000 /usr/local/bin/nord2000
COPY docker/aam /usr/local/bin/aam
RUN chmod +x /usr/local/bin/nord2000 /usr/local/bin/aam

# ---- Non-root user (Wine dislikes running as root) ----
# Mono is part of Wine used for running .NET applications.
# Gecko is Wine's Internet Explorer engine for rendering web pages in Windows apps.
# NMSim is Fortran binaries only — we disable Mono/Gecko via WINEDLLOVERRIDES below.
RUN useradd -m -s /bin/bash wineuser
USER wineuser
ENV HOME=/home/wineuser \
    WINEPREFIX=/home/wineuser/.wine \
    WINEARCH=win64 \
    WINEDEBUG=-all \
    WINEDLLOVERRIDES="mscoree,mshtml=" \
    NMSIM_HOME=/opt/nmsim \
    AAM_HOME=/opt/aam \
    PYTHONPATH=/repo

# Initialize the Wine prefix at build time so first run is fast (no install prompts).
RUN wineboot --init && wineserver -w || true

CMD ["/bin/bash"]
