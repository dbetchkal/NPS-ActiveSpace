"""Shared constants and helpers for setup tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from affine import Affine

from nps_active_space.setup import NMSIM_DST_CRS
from nps_active_space.setup.elevation import GRIDFLOAT_NODATA, NODATA_INT16, write_gridfloat
from nps_active_space.setup.site_writer import create_site_dir

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_PROJECT_DIR = REPO_ROOT / "example_data" / "site_projects"
STUDY_BOUNDS_4269 = (-136.2, 58.4, -135.8, 58.8)
# DENABULL-style bbox: mic in UTM zone 6, study-area western edge in zone 5.
DENABULL_BOUNDS_4269 = (-150.2, 63.1, -149.1, 63.6)
DENABULL_MIC_COORD = (-149.62901, 63.39192)
DENABULL_STUDYAREA_SW = (-150.02907, 63.21250)
DENABULL_STUDYAREA_NE = (-149.22895, 63.57134)


def write_source_dem(
    path: Path,
    bounds_4269: tuple[float, float, float, float],
    elevation: float,
    cellsize_deg: float = 0.02,
) -> None:
    """Write a minimal single-band DEM in NAD83 geographic (EPSG:4269)."""
    minx, miny, maxx, maxy = bounds_4269
    width = max(2, int(np.ceil((maxx - minx) / cellsize_deg)))
    height = max(2, int(np.ceil((maxy - miny) / cellsize_deg)))
    transform = Affine(cellsize_deg, 0.0, minx, 0.0, -cellsize_deg, maxy)
    data = np.full((height, width), elevation, dtype=np.float32)

    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs=NMSIM_DST_CRS,
        transform=transform,
        nodata=-9999.0,
    ) as ds:
        ds.write(data, 1)


def write_project_setup_elevation_artifacts(
    site_dir: Path,
    *,
    utm_suffix: str = "6",
    bounds_4269: tuple[float, float, float, float] = STUDY_BOUNDS_4269,
    elevation_m: int = 1524,
) -> tuple[Path, Path]:
    """
    Write minimal ``elevation_m_nad83_utm*.tif`` / ``.flt`` / ``.hdr`` under a site directory.

    Faster than ``setup_site`` when tests only need the artifact trio that
    ``get_project_setup_elevation`` discovers.
    """
    create_site_dir(site_dir)
    elev_dir = site_dir / "Input_Data" / "01_ELEVATION"
    base = elev_dir / f"elevation_m_nad83_utm{utm_suffix}"
    tif_path = base.with_suffix(".tif")
    minx, _, _, maxy = bounds_4269
    cellsize_deg = 0.02
    width = 2
    height = 2
    transform = Affine(cellsize_deg, 0.0, minx, 0.0, -cellsize_deg, maxy)
    data = np.full((height, width), elevation_m, dtype=np.int16)

    with rasterio.open(
        tif_path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="int16",
        crs=NMSIM_DST_CRS,
        transform=transform,
        nodata=int(NODATA_INT16),
    ) as ds:
        ds.write(data, 1)

    write_gridfloat(
        base,
        data,
        transform,
        width,
        height,
        NODATA_INT16,
        GRIDFLOAT_NODATA,
    )
    return tif_path, base.with_suffix(".flt")
