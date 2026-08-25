"""Cross-platform path helpers for the NPS-ActiveSpace project directory layout.

All filesystem paths in the pipeline should be built through these helpers (or
``os.path.join`` directly) rather than hard-coded backslashes or ``f"{a}/{b}"``
strings. Forward slashes often work on Windows too, but ``os.path.join`` is
explicit and keeps globs working on Linux/Mac.
"""
from __future__ import annotations

import glob
import os

from nps_active_space.active_space.propagation_model import (
    AAM_OUTPUT_SUBDIR,
    NMSIM_OUTPUT_SUBDIR,
)
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.utils.legacy_nmsim_paths import (
    NMSIM_ACTIVESPACES_SUBDIR,
    is_standard_altitude_layer_dir,
    resolve_activespace_geojson as _resolve_activespace_geojson,
    resolve_activespace_layer_dirs as _resolve_activespace_layer_dirs,
    resolve_nmsim_activespaces_dir,
)


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


def model_output_dir(site_dir_path: str, model: AcousticModel) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.NMSIM:
            subdir = NMSIM_OUTPUT_SUBDIR
        case AcousticModel.AAM:
            subdir = AAM_OUTPUT_SUBDIR
    return join(site_dir_path, subdir)


def model_activespaces_dir(site_dir_path: str, model: AcousticModel) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            return join(site_dir_path, AAM_OUTPUT_SUBDIR, "ACTIVESPACES")
        case AcousticModel.NMSIM:
            return join(site_dir_path, NMSIM_ACTIVESPACES_SUBDIR)


def activespaces_dir(project_dir: str, unit: str, site: str) -> str:
    """Legacy-compatible ACTIVESPACES root (prefers new NMSim layout when populated)."""
    return resolve_nmsim_activespaces_dir(site_dir(project_dir, unit, site))


def activespace_layer_dir(
    site_dir_path: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    model: AcousticModel,
) -> str:
    usy = deployment_id(unit, site, year)
    return join(model_activespaces_dir(site_dir_path, model), f"{usy}_{altitude_m}m")


def activespace_layer_dirs(
    project_dir: str,
    unit: str,
    site: str,
    year,
    model: AcousticModel = AcousticModel.NMSIM,
) -> list[str]:
    """Glob paths like ``.../ACTIVESPACES/DENATRLA2025_1000m``."""
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            usy = deployment_id(unit, site, year)
            pattern = join(
                model_activespaces_dir(
                    site_dir(project_dir, unit, site), AcousticModel.AAM,
                ),
                f"{usy}_*m",
            )
            matches = glob.glob(pattern)
            return sorted(
                d for d in matches if is_standard_altitude_layer_dir(d, usy)
            )
        case AcousticModel.NMSIM:
            return _resolve_activespace_layer_dirs(project_dir, unit, site, year)


def activespace_geojson(
    project_dir: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    gain_sign: str,
    gain_string: str,
    model: AcousticModel = AcousticModel.NMSIM,
) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            usy = deployment_id(unit, site, year)
            return join(
                model_activespaces_dir(
                    site_dir(project_dir, unit, site), AcousticModel.AAM,
                ),
                f"{usy}_{altitude_m}m",
                f"{usy}_O_{gain_sign}{gain_string}.geojson",
            )
        case AcousticModel.NMSIM:
            return _resolve_activespace_geojson(
                project_dir, unit, site, year, altitude_m, gain_sign, gain_string,
            )


def precision_recall_dir(site_dir_path: str, model: AcousticModel) -> str:
    return join(model_output_dir(site_dir_path, model), "PRECISION_RECALL")


def site_fits_csv(site_dir_path: str) -> str:
    return join(site_dir_path, "fits.csv")


def tested_points_dir(
    site_dir_path: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    model: AcousticModel,
) -> str:
    usy = deployment_id(unit, site, year)
    return join(
        model_output_dir(site_dir_path, model),
        "TESTED_POINTS",
        f"{usy}_{altitude_m}m",
    )


def annotation_files(project_dir: str, unit: str, site: str, year) -> list[str]:
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


def precision_recall_3d_plot(
    project_dir: str,
    unit: str,
    site: str,
    year,
    model: AcousticModel = AcousticModel.NMSIM,
    beta: float = 1.0,
) -> str:
    usy = deployment_id(unit, site, year)
    beta_str = str(beta).replace(".", "p")
    site_dir_path = site_dir(project_dir, unit, site)
    return join(
        model_output_dir(site_dir_path, model),
        "PRECISION_RECALL",
        f"PrecisionRecallPlot_3d_{usy}_f{beta_str}.png",
    )


def altitude_histogram_plot(project_dir: str, unit: str, site: str, year) -> str:
    usy = deployment_id(unit, site, year)
    return join(site_dir(project_dir, unit, site), f"Altitude_Histogram_{usy}.png")


def all_annotation_files(project_dir: str) -> list[str]:
    """Find annotation geojson files under each site directory in ``project_dir``."""
    return glob.glob(join(project_dir, "*", "*saved_annotations*.geojson"))
