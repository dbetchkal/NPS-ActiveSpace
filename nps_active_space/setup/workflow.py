"""
Orchestrate the full NMSIM project setup workflow.

This requires two conceptual inputs: (1) a study area polygon, and (2) a listening location point.
We use the study area to prepare a portion of a Digital Elevation Model for use in computation and visualization.
The elevation data is also converted to an antiquated ESRI GridFloat format required by ``NMSIM``.
We represent the listening location coordinates as a canonical NMSIM site (.sit) file.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Point

from nps_active_space.setup.elevation import (
    GRIDFLOAT_NODATA,
    NMSIM_DST_CRS,
    NODATA_INT16,
    create_elevation_tif,
    write_gridfloat,
)
from nps_active_space.setup.site import (
    create_site_dir,
    create_site_file,
    create_study_area,
    deployment_sit_name,
    sit_file_path,
)
from nps_active_space.utils.computation import coords_to_utm
from nps_active_space.utils.enums import SourceElevationUnits


@dataclass(frozen=True)
class SetupSiteResult:
    """Paths and in-memory artifacts produced by :func:`setup_site`."""

    study_area_wgs84: gpd.GeoDataFrame
    study_area_path: Path
    elevation_tif_path: Path
    sit_path: Path


def setup_site(
    site_dir: str | Path,
    unit: str,
    site: str,
    year: int,
    mic_coord: tuple[float, float],
    studyarea_sw: tuple[float, float],
    studyarea_ne: tuple[float, float],
    source_dem: str | Path,
    dst_crs: str = NMSIM_DST_CRS,
    mic_height_agl_m: float = 1.5,
    source_elevation_units: SourceElevationUnits = SourceElevationUnits.FEET,
) -> SetupSiteResult:
    """
    Run the full NMSIM site setup workflow for a new deployment directory.

    Creates the canonical site directory tree, study area shapefile, clipped/reprojected
    elevation GeoTIFF, GridFloat (.flt/.hdr), and microphone .sit file.

    Parameters
    ----------
    site_dir : str | Path
        Canonical NMSIM site directory (e.g. ``.../DENATRLA``).
    unit : str
        Four-character NPS Alpha Code, e.g. ``DENA``.
    site : str
        Alpha-numeric acoustic monitoring site code, e.g. ``TRLA``.
    year : int
        Four-digit deployment year, e.g. 2018.
    mic_coord : tuple[float, float]
        Microphone location as ``(lon, lat)`` in WGS84 decimal degrees.
    studyarea_sw : tuple[float, float]
        Study area southwest corner as ``(lon, lat)`` in WGS84.
    studyarea_ne : tuple[float, float]
        Study area northeast corner as ``(lon, lat)`` in WGS84.
    source_dem : str | Path
        Path to the source DEM raster. Elevation values must be in ``source_elevation_units``.
    dst_crs : str, default ``EPSG:4269``
        Coordinate reference system for the output elevation GeoTIFF.
    mic_height_agl_m : float, default 1.5
        Microphone height above ground in meters (ANSI standard).
    source_elevation_units : SourceElevationUnits, default ``FEET``
        Vertical units of ``source_dem``. ``project_setup.py`` reads ``dem_elevation_units``
        from the environment config; library callers may override.

    Returns
    -------
    SetupSiteResult
        Study area geometry plus paths to written study area, elevation, and ``.sit`` files.
    """
    site_dir = Path(site_dir)
    create_site_dir(site_dir)

    study_area_wgs84 = create_study_area(
        site_dir=site_dir,
        unit=unit,
        site=site,
        year=year,
        study_area_sw_corner=studyarea_sw,
        study_area_ne_corner=studyarea_ne,
    )

    # to create a .sit file, we will eventually need the UTM coordinates of the microphone...
    utm_epsg, _ = coords_to_utm(lat=mic_coord[1], lon=mic_coord[0])
    mic_utm = gpd.GeoSeries([Point(mic_coord)], crs="EPSG:4326").to_crs(utm_epsg)

    # quirk: `NMSIM` expects the elevation input's spatial reference to be the UTM Zone of the westernmost exent
    #         this is occasionally different than the microphone's UTM Zone in the .sit file
    _, utm_zone_str = coords_to_utm(lat=studyarea_sw[1], lon=studyarea_sw[0])
    output_base = site_dir / "Input_Data" / "01_ELEVATION" / f"elevation_m_nad83_utm{utm_zone_str}"
    output_tif = output_base.with_suffix(".tif")

    dst_int16, dst_transform, dst_width, dst_height = create_elevation_tif(
        source_dem=source_dem,
        study_area_wgs84=study_area_wgs84,
        dst_crs=dst_crs,
        output_tif=output_tif,
        nodata_int16=NODATA_INT16,
        source_elevation_units=source_elevation_units,
    )

    # TODO when a National Landcover Database (NLCD) mapping has been construed, an equivalent
    # `create_NMSIM_landcover_tif` function would belong here...

    write_gridfloat(
        output_base=output_base,
        dst_int16=dst_int16,
        dst_transform=dst_transform,
        dst_width=dst_width,
        dst_height=dst_height,
        nodata_int16=NODATA_INT16,
        gridfloat_nodata=GRIDFLOAT_NODATA,
    )

    create_site_file(
        site_dir=site_dir,
        unit=unit,
        site=site,
        year=year,
        easting_m=mic_utm.x.item(),
        northing_m=mic_utm.y.item(),
        height=mic_height_agl_m,
    )

    return SetupSiteResult(
        study_area_wgs84=study_area_wgs84,
        study_area_path=site_dir / f"{unit}{site}_study_area.shp",
        elevation_tif_path=output_tif,
        sit_path=sit_file_path(site_dir, deployment_sit_name(unit, site, year)),
    )
