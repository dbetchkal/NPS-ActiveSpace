"""
Complete DENATRLA (TRLA) example elevation artifacts and optional TRLA-only source DEM.

Run from the repository root:

    python example_data/scripts/regenerate_denatrla_elevation.py
    python example_data/scripts/regenerate_denatrla_elevation.py --export-source-dem

The default action writes GridFloat (.flt/.hdr) from the checked-in
``elevation_m_nad83_utm6.tif`` and refreshes ``DENATRLA2025.sit``. This preserves
real topography without re-clipping from a global DEM.

``--export-source-dem`` additionally writes ``example_data/source_dems/DENATRLA_trla_m.tif``
(float32 **meters**, NAD83 geographic) for optional ``project_setup`` round-trips on TRLA only.
Use with ``dem_elevation_units = meters`` in ``DENA_example.config``. Do not set ``feet`` on
meter data — that would apply ``× 0.3048`` again and crush elevations.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import rasterio

from nps_active_space.setup.elevation import GRIDFLOAT_NODATA, NODATA_INT16, write_gridfloat
from nps_active_space.setup.site_writer import write_listener_site_file
from nps_active_space.utils.computation import NMSIM_bbox_utm
from nps_active_space.utils.helpers import get_deployment, load_studyarea

REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_DIR = REPO_ROOT / "example_data/site_projects"
SITE_DIR = PROJECT_DIR / "DENATRLA"
ELEV_TIF = SITE_DIR / "Input_Data/01_ELEVATION/elevation_m_nad83_utm6.tif"
SOURCE_DEM_DIR = REPO_ROOT / "example_data/source_dems"
SOURCE_DEM_PATH = SOURCE_DEM_DIR / "DENATRLA_trla_m.tif"


def complete_gridfloat_and_sit() -> None:
    if not ELEV_TIF.is_file():
        raise FileNotFoundError(f"Expected checked-in elevation GeoTIFF at {ELEV_TIF}")

    output_base = ELEV_TIF.with_suffix("")
    with rasterio.open(ELEV_TIF) as ds:
        data = ds.read(1)
        nodata = np.int16(ds.nodata)
        write_gridfloat(
            output_base,
            data,
            ds.transform,
            ds.width,
            ds.height,
            nodata,
            GRIDFLOAT_NODATA,
        )

    study_area = load_studyarea(str(PROJECT_DIR), "DENA", "TRLA", 2025)
    mic = get_deployment(str(PROJECT_DIR), "DENA", "TRLA", 2025, elevation=False)
    projected_mic = mic.to_crs(NMSIM_bbox_utm(study_area))
    flt_rel = output_base.with_suffix(".flt").relative_to(SITE_DIR)

    write_listener_site_file(
        SITE_DIR,
        mic.name,
        projected_mic.x,
        projected_mic.y,
        projected_mic.z,
        flt_rel,
    )

    print(f"Wrote {output_base.with_suffix('.flt')}")
    print(f"Wrote {output_base.with_suffix('.hdr')}")
    print(f"Updated {SITE_DIR / 'Input_Data/05_SITES/DENATRLA2025.sit'}")


def export_trla_source_dem() -> None:
    """Export a TRLA-only meter DEM from the checked-in GeoTIFF for ``project_setup``."""
    SOURCE_DEM_DIR.mkdir(parents=True, exist_ok=True)

    with rasterio.open(ELEV_TIF) as ds:
        data_m = ds.read(1)
        nodata_m = np.int16(ds.nodata)
        meters = np.full(data_m.shape, -9999.0, dtype=np.float32)
        valid = data_m != nodata_m
        meters[valid] = data_m[valid]

        profile = ds.profile.copy()
        profile.update(dtype="float32", nodata=-9999.0, compress="lzw")

        with rasterio.open(SOURCE_DEM_PATH, "w", **profile) as out:
            out.write(meters, 1)

    print(f"Wrote TRLA source DEM (meters) to {SOURCE_DEM_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--export-source-dem",
        action="store_true",
        help="Also export example_data/source_dems/DENATRLA_trla_m.tif (meters) for project_setup.",
    )
    args = parser.parse_args()

    complete_gridfloat_and_sit()
    if args.export_source_dem:
        export_trla_source_dem()


if __name__ == "__main__":
    main()
