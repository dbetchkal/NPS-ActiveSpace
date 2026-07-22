"""NMSIM site directory, study area, and .sit file helpers."""

from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

NMSIM_SITES_DIR = "Input_Data/05_SITES"

NMSIM_SITE_SUBFOLDERS = [
    "Input_Data",
    "Input_Data/01_ELEVATION",
    "Input_Data/02_IMPEDANCE",
    "Input_Data/03_TRAJECTORY",
    "Input_Data/04_LAYERS",
    NMSIM_SITES_DIR,
    "Input_Data/06_AMBIENCE",
    "Input_Data/07_WEATHER",
    "Input_Data/08_TREES",
    "Output_Data",
    "Output_Data/ASCII",
    "Output_Data/IMAGES",
    "Output_Data/SITE",
    "Output_Data/TIG_TIS",
]


def deployment_sit_name(unit: str, site: str, year: int) -> str:
    """Return the canonical NMSIM ``.sit`` basename (without extension) for a deployment."""
    return f"{unit}{site}{year}"


def sit_file_path(site_dir: str | Path, sit_name: str) -> Path:
    """Return ``{site_dir}/Input_Data/05_SITES/{sit_name}.sit``."""
    return Path(site_dir) / NMSIM_SITES_DIR / f"{sit_name}.sit"


def create_site_dir(site_dir: str | Path) -> None:
    """
    Create a canonical NMSIM site directory structure expected by downstream NPS-ActiveSpace modules.

    Inputs
    ------
    site_dir : str | Path
        Base directory where the canonical site structure will be created.

    Returns
    -------
    None
    """
    site_path = Path(site_dir)
    for folder in NMSIM_SITE_SUBFOLDERS:
        (site_path / folder).mkdir(parents=True, exist_ok=True)


def create_study_area(
    site_dir: str | Path,
    unit: str,
    site: str,
    year: int,
    study_area_sw_corner: tuple[float, float],
    study_area_ne_corner: tuple[float, float],
) -> gpd.GeoDataFrame:
    """
    Create the most important conceptual input of ``NPS-ActiveSpace``: a user-defined study area.

    Builds a rectangular polygon from the provided SW/NE corners and saves it as an ESRI shapefile (.shp).

    Inputs
    ------
    site_dir : str | Path
        A canonical NMSIM project directory.
    unit : str
        4-character NPS Alpha Code, e.g. ``BITH``, ``YUCH``.
    site : str
        Alpha-numeric acoustic monitoring site code, e.g., ``002``, ``TRLA``.
    year : int
        Four digit deployment year, e.g. 2018.
    study_area_sw_corner : tuple
        Study area SW corner x y (WGS84). Example: (-136.088360, 58.569310).
    study_area_ne_corner : tuple
        Study area NE corner x y (WGS84). Example: (-135.818994, 58.706095).

    Returns
    -------
    study_area_wgs84 : gpd.GeoDataFrame
        The rectangular study area geometry (WGS84), with columns ``["Unit", "Site", "Year"]``.
    """
    # Input args are WGS84 lon/lat; shapely box expects x1,y1,x2,y2 => lon/lat
    study_area_wgs84 = gpd.GeoDataFrame(
        [[unit, site, year]],
        geometry=[box(*study_area_sw_corner, *study_area_ne_corner)],
        crs="EPSG:4326",
        columns=["Unit", "Site", "Year"],
    )
    study_area_path = Path(site_dir) / (f"{unit}{site}_study_area.shp")
    study_area_proj = study_area_wgs84.to_crs("EPSG:4269")
    study_area_proj.to_file(study_area_path)
    return study_area_wgs84


def write_site_file(
    sit_path: str | Path,
    easting_m: float,
    northing_m: float,
    height_agl_m: float,
    label: str,
    dem_flt_path: str | Path,
) -> None:
    """
    Write the body of an NMSIM ``.sit`` listener-location file.

    Low-level helper used by :func:`create_site_file` and ``ActiveSpaceGenerator``.
    """
    sit_path = Path(sit_path)
    with open(sit_path, "w") as site_file:
        site_file.write("    0\n")
        site_file.write("    1\n")
        site_file.write(
            "{0:19.0f}.{1:9.0f}.{2:10.5f} {3:20}\n".format(
                easting_m, northing_m, height_agl_m, label
            )
        )
        site_file.write(str(dem_flt_path) + "\n")


def write_listener_site_file(
    site_dir: str | Path,
    sit_name: str,
    easting_m: float,
    northing_m: float,
    height_agl_m: float,
    dem_flt_path: str | Path,
    *,
    label: str | None = None,
) -> Path:
    """Write ``Input_Data/05_SITES/{sit_name}.sit`` for an NMSIM listener location."""
    sit_path = sit_file_path(site_dir, sit_name)
    write_site_file(sit_path, easting_m, northing_m, height_agl_m, label or sit_name, dem_flt_path)
    return sit_path


def create_site_file(
    site_dir: str | Path,
    unit: str,
    site: str,
    year: int,
    easting_m: float,
    northing_m: float,
    height: float = 1.5,
    dem_flt_path: str | Path | None = None,
) -> None:
    """
    Create an NMSIM site (.sit) file for a given NPS microphone deployment.
    The .sit file represents a "listener location" (microphone location), one of two conceptual inputs to NPS-ActiveSpace.

    Inputs
    ------
    site_dir : str | Path
        A path location of a canonical NMSIM site directory.
    unit : str
        4-character NPS Alpha Code, e.g. ``BITH``, ``YUCH``.
    site : str
        Alpha-numeric acoustic monitoring site code, e.g., ``002``, ``TRLA``.
    year : int
        Four digit deployment year, e.g. 2018.
    easting_m : float
        UTM easting in meters (microphone's UTM zone).
    northing_m : float
        UTM northing in meters (microphone's UTM zone).
    height : float
        Microphone height in meters; defaults to ANSI standard, 1.5 meters.
    dem_flt_path : str | Path, optional
        Path to the GridFloat ``.flt`` referenced by the ``.sit`` file. When omitted,
        the lexicographically first ``*.flt`` in ``Input_Data/01_ELEVATION`` is used.

    Returns
    -------
    None
    """
    # deployment location can subtly vary by year
    if dem_flt_path is None:
        elev_dir = Path(site_dir) / "Input_Data" / "01_ELEVATION"
        matches = sorted(elev_dir.glob("*.flt"))
        if not matches:
            raise FileNotFoundError(f"No .flt files found in {elev_dir}")
        dem_flt_path = matches[0]

    write_listener_site_file(
        site_dir,
        deployment_sit_name(unit, site, year),
        easting_m,
        northing_m,
        height,
        dem_flt_path,
        label=unit + site,
    )


def parse_sit_coords_line(line: str) -> tuple[float, float, float]:
    """Parse easting, northing, and height (m) from a ``.sit`` coordinates line."""
    coords_str = re.split(r"\s+", line.strip())[0:3]
    return float(coords_str[0]), float(coords_str[1]), float(coords_str[2])
