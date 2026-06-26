"""Cross-platform path helpers for the NPS-ActiveSpace project directory layout.

All filesystem paths in the pipeline should be built through these helpers (or
``os.path.join`` directly) rather than hard-coded backslashes or ``f"{a}/{b}"``
strings. Forward slashes often work on Windows too, but ``os.path.join`` is
explicit and keeps globs working on Linux/Mac.
"""
import glob
import os
from typing import List


def join(*parts: str) -> str:
    return os.path.join(*parts)


def deployment_id(unit: str, site: str, year) -> str:
    return f"{unit}{site}{year}"


def site_dir(project_dir: str, unit: str, site: str) -> str:
    return join(project_dir, f"{unit}{site}")


def site_path(project_dir: str, unit: str, site: str, *parts: str) -> str:
    return join(project_dir, f"{unit}{site}", *parts)


def input_data_dir(project_dir: str, unit: str, site: str) -> str:
    return site_path(project_dir, unit, site, "Input_Data")


def output_data_dir(project_dir: str, unit: str, site: str) -> str:
    return site_path(project_dir, unit, site, "Output_Data")


def activespaces_dir(project_dir: str, unit: str, site: str) -> str:
    return site_path(project_dir, unit, site, "Output_Data", "ACTIVESPACES")


def activespace_layer_dirs(project_dir: str, unit: str, site: str, year) -> List[str]:
    """Glob paths like ``.../ACTIVESPACES/DENATRLA2025_1000m``."""
    pattern = join(activespaces_dir(project_dir, unit, site), f"{deployment_id(unit, site, year)}_*m")
    return glob.glob(pattern)


def activespace_geojson(project_dir: str, unit: str, site: str, year, altitude_m: int,
                        gain_sign: str, gain_string: str) -> str:
    usy = deployment_id(unit, site, year)
    return join(activespaces_dir(project_dir, unit, site), f"{usy}_{altitude_m}m",
                f"{usy}_O_{gain_sign}{gain_string}.geojson")


def annotation_files(project_dir: str, unit: str, site: str, year) -> List[str]:
    """Annotation geojson files saved in the site directory root (ground-truthing output)."""
    usy = deployment_id(unit, site, year)
    pattern = join(site_dir(project_dir, unit, site), f"{usy}*saved_annotations*.geojson")
    return glob.glob(pattern)


def study_area_shapefile(project_dir: str, unit: str, site: str) -> str:
    site = site_dir(project_dir, unit, site)
    matches = glob.glob(join(site, "*study*.shp"))
    if not matches:
        matches = glob.glob(join(site, f"{unit}{site}*study*area*.shp"))
    if not matches:
        raise FileNotFoundError(f"No study area shapefile under {site}")
    return matches[0]


def dem_raster(project_dir: str, unit: str, site: str) -> str:
    elevation_dir = join(input_data_dir(project_dir, unit, site), "01_ELEVATION")
    matches = glob.glob(join(elevation_dir, "elevation_m_nad83_utm*.tif"))
    if not matches:
        raise FileNotFoundError(f"No elevation_m_nad83_utm*.tif under {elevation_dir}")
    return matches[0]


def site_file(project_dir: str, unit: str, site: str, year) -> str:
    return join(input_data_dir(project_dir, unit, site), "05_SITES",
                f"{unit}{site}{year}.sit")


def fits_csv(project_dir: str) -> str:
    return join(project_dir, "fits.csv")


def precision_recall_3d_plot(project_dir: str, unit: str, site: str, year) -> str:
    usy = deployment_id(unit, site, year)
    return join(site_dir(project_dir, unit, site), f"precision_recall_3d_f1_{usy}.png")


def altitude_histogram_plot(project_dir: str, unit: str, site: str, year) -> str:
    usy = deployment_id(unit, site, year)
    return join(site_dir(project_dir, unit, site), f"Altitude_Histogram_{usy}.png")


def all_annotation_files(project_dir: str) -> List[str]:
    """Find annotation geojson files under each site directory in ``project_dir``."""
    return glob.glob(join(project_dir, "*", "*saved_annotations*.geojson"))
