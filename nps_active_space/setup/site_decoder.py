"""Decode ``.sit`` UTM coordinates to geographic locations and diagnose legacy zone mismatches."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import geopandas as gpd
from pyproj import Transformer
from shapely.geometry import Point

from nps_active_space.utils.computation import NMSIM_bbox_utm, coords_to_utm

logger = logging.getLogger(__name__)

SitCoordsStatus = Literal["ok", "legacy", "failed"]


@dataclass(frozen=True)
class SitFileContents:
    """Parsed body of an NMSIM ``.sit`` listener-location file."""

    header_line_1: str
    header_line_2: str
    easting_m: float
    northing_m: float
    height_agl_m: float
    label: str
    dem_flt_path: str


def parse_sit_coords_line(line: str) -> tuple[float, float, float]:
    """Parse easting, northing, and height (m) from a ``.sit`` coordinates line."""
    coords_str = re.split(r"\s+", line.strip())[0:3]
    return float(coords_str[0]), float(coords_str[1]), float(coords_str[2])


def parse_sit_label(coords_line: str) -> str:
    """Parse the trailing listener label from a ``.sit`` coordinates line."""
    tokens = coords_line.strip().split()
    if len(tokens) < 4:
        return ""
    return tokens[3]


def read_sit_file(sit_path: str | Path) -> SitFileContents:
    """Read and parse an NMSIM ``.sit`` listener-location file."""
    sit_path = Path(sit_path)
    lines = sit_path.read_text().splitlines()
    if len(lines) < 4:
        raise ValueError(
            f"Microphone site file has invalid format (expected coordinates on line 3):\n"
            f"  {sit_path}"
        )
    coords_line = lines[2]
    easting_m, northing_m, height_agl_m = parse_sit_coords_line(coords_line)
    return SitFileContents(
        header_line_1=lines[0],
        header_line_2=lines[1],
        easting_m=easting_m,
        northing_m=northing_m,
        height_agl_m=height_agl_m,
        label=parse_sit_label(coords_line),
        dem_flt_path=lines[3],
    )


def study_area_corner_coords_wgs84(
    study_area: gpd.GeoDataFrame,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return study-area SW and NE corners as ``(lon, lat)`` in WGS84."""
    minx, miny, maxx, maxy = study_area.to_crs("EPSG:4326").total_bounds
    return (minx, miny), (maxx, maxy)


def format_project_setup_command(
    unit: str,
    site: str,
    year: int,
    mic_coord: tuple[float, float],
    studyarea_sw: tuple[float, float],
    studyarea_ne: tuple[float, float],
) -> str:
    """Return a copy-paste ``project_setup`` invocation for the given deployment."""
    mic_lon, mic_lat = mic_coord
    sw_lon, sw_lat = studyarea_sw
    ne_lon, ne_lat = studyarea_ne
    return (
        f"python -m nps_active_space.scripts.project_setup -e <environment> "
        f"-u {unit} -s {site} -y {year} "
        f"--mic-coord {mic_lon:.8f} {mic_lat:.8f} "
        f"--studyarea-sw {sw_lon:.8f} {sw_lat:.8f} "
        f"--studyarea-ne {ne_lon:.8f} {ne_lat:.8f}"
    )


def _utm_epsg_candidates_for_study_area(study_area: gpd.GeoDataFrame) -> list[str]:
    """UTM EPSG codes that may span the study-area bbox, with project zone first."""
    study_nad83 = study_area if study_area.crs.to_epsg() == 4269 else study_area.to_crs(4269)
    minx, miny, maxx, maxy = study_nad83.total_bounds
    project_utm = NMSIM_bbox_utm(study_area)
    candidates: list[str] = []

    for lon, lat in ((minx, miny), (minx, maxy), (maxx, miny), (maxx, maxy)):
        epsg, _ = coords_to_utm(lat=lat, lon=lon)
        if epsg not in candidates:
            candidates.append(epsg)

    if project_utm in candidates:
        candidates.remove(project_utm)
    return [project_utm, *candidates]


def _lonlat_in_study_area(lon: float, lat: float, study_area: gpd.GeoDataFrame) -> bool:
    study_polygon = study_area.to_crs("EPSG:4326").union_all()
    return study_polygon.contains(Point(lon, lat))


def _sit_utm_to_lonlat(
    easting_m: float,
    northing_m: float,
    utm_epsg: str,
) -> tuple[float, float]:
    """Transform ``.sit`` easting/northing in ``utm_epsg`` to WGS84 lon/lat."""
    return Transformer.from_crs(utm_epsg, "EPSG:4326", always_xy=True).transform(
        easting_m, northing_m
    )


def _decode_sit_lonlat(
    easting_m: float,
    northing_m: float,
    study_area: gpd.GeoDataFrame,
) -> tuple[float, float, str] | None:
    """Return ``(lon, lat, utm_epsg)`` for the first zone placing the point in the study area."""
    for utm_epsg in _utm_epsg_candidates_for_study_area(study_area):
        lon, lat = _sit_utm_to_lonlat(easting_m, northing_m, utm_epsg)
        if _lonlat_in_study_area(lon, lat, study_area):
            return lon, lat, utm_epsg
    return None


def diagnose_sit_coords(
    easting_m: float,
    northing_m: float,
    study_area: gpd.GeoDataFrame,
) -> tuple[SitCoordsStatus, float, float, str, str | None]:
    """
    Classify ``.sit`` easting/northing without raising.

    Returns ``(status, lon, lat, project_utm, decoded_utm)``.
    """
    project_utm = NMSIM_bbox_utm(study_area)
    decoded = _decode_sit_lonlat(easting_m, northing_m, study_area)
    if decoded is not None:
        lon, lat, decoded_utm = decoded
        if decoded_utm == project_utm:
            return "ok", lon, lat, project_utm, decoded_utm
        return "legacy", lon, lat, project_utm, decoded_utm

    lon, lat = _sit_utm_to_lonlat(easting_m, northing_m, project_utm)
    return "failed", lon, lat, project_utm, None


def project_zone_utm_coords(
    lon: float,
    lat: float,
    study_area: gpd.GeoDataFrame,
) -> tuple[float, float]:
    """Project WGS84 lon/lat to the NMSIM project UTM zone for ``study_area``."""
    project_utm = NMSIM_bbox_utm(study_area)
    easting_m, northing_m = Transformer.from_crs(
        "EPSG:4326", project_utm, always_xy=True
    ).transform(lon, lat)
    return easting_m, northing_m


def _raise_sit_coords_error(
    *,
    sit_path: Path,
    headline: str,
    detail: str,
    unit: str,
    site: str,
    year: int,
    mic_coord: tuple[float, float],
    study_area: gpd.GeoDataFrame,
) -> None:
    """Log a warning and raise with a copy-paste ``project_setup`` command."""
    studyarea_sw, studyarea_ne = study_area_corner_coords_wgs84(study_area)
    setup_cmd = format_project_setup_command(
        unit,
        site,
        year,
        mic_coord=mic_coord,
        studyarea_sw=studyarea_sw,
        studyarea_ne=studyarea_ne,
    )
    message = (
        f"{headline}\n"
        f"  {sit_path}\n"
        f"{detail}"
        f"Re-run project setup to rewrite the .sit file:\n"
        f"  {setup_cmd}"
    )
    logger.warning(message)
    raise ValueError(message)


def decode_sit_geographic_coords(
    easting_m: float,
    northing_m: float,
    study_area: gpd.GeoDataFrame,
    *,
    sit_path: str | Path,
    unit: str,
    site: str,
    year: int,
) -> tuple[float, float]:
    """
    Decode ``.sit`` easting/northing to WGS84 lon/lat using the NMSIM project UTM zone.

    Raises
    ------
    ValueError
        If coords appear to use a legacy mic-local UTM zone, with a ``project_setup`` command to rerun.
    """
    sit_path = Path(sit_path)
    status, lon, lat, project_utm, decoded_utm = diagnose_sit_coords(
        easting_m, northing_m, study_area
    )
    if status == "ok":
        return lon, lat
    if status == "legacy":
        _raise_sit_coords_error(
            sit_path=sit_path,
            headline="Microphone site file uses legacy mic-local UTM coordinates:",
            detail=(
                f"Decoded as {decoded_utm}, but NMSIM expects the project zone {project_utm} "
                f"(western edge of the study area).\n"
            ),
            unit=unit,
            site=site,
            year=year,
            mic_coord=(lon, lat),
            study_area=study_area,
        )

    _raise_sit_coords_error(
        sit_path=sit_path,
        headline=(
            "Microphone site file coordinates do not fall inside the study area when decoded "
            f"with the NMSIM project UTM zone {project_utm}:"
        ),
        detail=(
            f"Decoded location: lon={lon:.6f}, lat={lat:.6f}\n"
            "This deployment may have been created before the project-zone .sit fix.\n"
        ),
        unit=unit,
        site=site,
        year=year,
        mic_coord=(lon, lat),
        study_area=study_area,
    )
