"""NMSIM project setup: study area, elevation, site artifacts, and workflow orchestration."""

from nps_active_space.setup.elevation import (
    FEET_TO_METERS,
    GRIDFLOAT_NODATA,
    NMSIM_DST_CRS,
    NODATA_INT16,
    create_elevation_tif,
    parse_source_elevation_units,
    write_gridfloat,
)
from nps_active_space.setup.site import (
    NMSIM_SITE_SUBFOLDERS,
    NMSIM_SITES_DIR,
    create_site_dir,
    create_site_file,
    create_study_area,
    deployment_sit_name,
    parse_sit_coords_line,
    sit_file_path,
    write_listener_site_file,
    write_site_file,
)
from nps_active_space.setup.workflow import SetupSiteResult, setup_site
from nps_active_space.utils.enums import SourceElevationUnits

# Backward-compatible aliases for project_setup.py function names (PR #75).
create_NMSIM_site_dir = create_site_dir
create_NMSIM_site_file = create_site_file
create_NMSIM_elevation_tif = create_elevation_tif
create_gridfloat = write_gridfloat

__all__ = [
    "FEET_TO_METERS",
    "GRIDFLOAT_NODATA",
    "NMSIM_DST_CRS",
    "NMSIM_SITE_SUBFOLDERS",
    "NMSIM_SITES_DIR",
    "NODATA_INT16",
    "SourceElevationUnits",
    "SetupSiteResult",
    "create_NMSIM_elevation_tif",
    "create_NMSIM_site_dir",
    "create_NMSIM_site_file",
    "create_elevation_tif",
    "create_gridfloat",
    "create_site_dir",
    "create_site_file",
    "create_study_area",
    "deployment_sit_name",
    "parse_sit_coords_line",
    "parse_source_elevation_units",
    "setup_site",
    "sit_file_path",
    "write_gridfloat",
    "write_listener_site_file",
    "write_site_file",
]
