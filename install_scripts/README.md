# Install scripts

Dev-environment bootstrap at the repo root. These are **not** the scientific workflow CLIs under [`nps_active_space/scripts/`](../nps_active_space/scripts/README.md).

| Script | Purpose |
|--------|---------|
| [`install_gdal_python_binding.sh`](install_gdal_python_binding.sh) | macOS/Linux: build PyPI `GDAL` against system `libgdal`, then `pip install -e ".[dev,aam]"` (recommended) |
