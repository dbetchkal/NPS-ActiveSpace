"""DEM clip/reproject and ESRI GridFloat output for NMSIM."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from affine import Affine
from rasterio.mask import mask as rio_mask
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import mapping

from nps_active_space.utils.enums import FEET_TO_METERS, SourceElevationUnits

NODATA_INT16 = np.int16(32767)
GRIDFLOAT_NODATA = np.float32(-32768.0)
NMSIM_DST_CRS = "EPSG:4269"  # NAD 83 - https://epsg.io/4269


def parse_source_elevation_units(units: str) -> SourceElevationUnits:
    """Validate and normalize a ``dem_elevation_units`` config value."""
    return SourceElevationUnits.parse_config_value(units)


def _to_meters_elevation(values: np.ndarray, source_units: SourceElevationUnits) -> np.ndarray:
    """Convert raster elevation values from ``source_units`` to integer meters."""
    return np.rint(values * source_units.to_meters_scale())


def create_elevation_tif(
    source_dem: str | Path,
    study_area_wgs84: gpd.GeoDataFrame,
    dst_crs: str,
    output_tif: str | Path,
    nodata_int16: np.int16,
    *,
    source_elevation_units: SourceElevationUnits = SourceElevationUnits.FEET,
) -> tuple[np.ndarray, Affine, int, int]:
    """
    Clip a feet- or meters-elevation source DEM, reproject, and write int16 meters for ``NMSIM``.

    Reads a source raster, clips it to ``study_area_wgs84``, reprojects to ``dst_crs``,
    converts elevation to integer meters, and saves a GeoTIFF used downstream for
    computation and visualization.

    Inputs
    ------
    source_dem : str | Path
        Path to the input DEM raster. Elevation values must be in ``source_elevation_units``.
    study_area_wgs84 : gpd.GeoDataFrame
        Study area polygon (WGS84).
    dst_crs : str
        Coordinate reference system for the destination raster.
    output_tif : str | Path
        Path for output GeoTIFF file (elevation values in meters).
    nodata_int16 : np.int16
        The int16 value used as NODATA in the output GeoTIFF.
    source_elevation_units : SourceElevationUnits, default ``FEET``
        Vertical units of ``source_dem`` before conversion to meters. Scripts read
        ``dem_elevation_units`` from the environment config; library callers may override.

    Returns
    -------
    dst_int16 : np.ndarray
        A 2D array (dtype int16) of clipped/reprojected elevations in meters.
    dst_transform : affine.Affine
        The affine transformation matrix representing reprojection from WGS84 to NAD83 GCS North American.
    dst_width : int
        Number of columns in the output raster.
    dst_height : int
        Number of rows in the output raster.

    Notes
    -----
    ``NMSIM`` requires a 16-bit signed integer GeoTIFF with elevations in meters.
    """
    # Clip using source DEM CRS
    with rasterio.open(source_dem) as src:
        study_in_src = study_area_wgs84.to_crs(src.crs)
        clipped, clipped_transform = rio_mask(
            src,
            [mapping(g) for g in study_in_src.geometry],
            crop=True,
            filled=False,
        )
        clipped = clipped[0]
        src_profile = src.profile
        src_crs = src.crs

    # Compute output grid geometry
    left = clipped_transform.c
    top = clipped_transform.f
    right = left + clipped.shape[1] * clipped_transform.a
    bottom = top + clipped.shape[0] * clipped_transform.e

    dst_transform, dst_width, dst_height = calculate_default_transform(
        src_crs,
        dst_crs,
        clipped.shape[1],
        clipped.shape[0],
        left,
        bottom,
        right,
        top,
    )

    dst = np.full((dst_height, dst_width), nodata_int16, dtype=np.float32)

    reproject(
        source=clipped.filled(src_profile["nodata"]),
        destination=dst,
        src_transform=clipped_transform,
        src_crs=src_crs,
        src_nodata=src_profile["nodata"],
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        dst_nodata=np.float32(nodata_int16),
        resampling=Resampling.bilinear,
    )

    nodata_mask = dst == nodata_int16
    dst[~nodata_mask] = _to_meters_elevation(dst[~nodata_mask], source_elevation_units)
    dst_int16 = dst.astype(np.int16)

    profile = dict(
        driver="GTiff",
        dtype=rasterio.int16,
        count=1,
        compress="lzw",
        width=dst_width,
        height=dst_height,
        transform=dst_transform,
        crs=dst_crs,
        nodata=nodata_int16,
    )

    Path(output_tif).parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(output_tif, "w", **profile) as ds:
        ds.write(dst_int16, 1)

    return dst_int16, dst_transform, dst_width, dst_height


def write_gridfloat(
    output_base: str | Path,
    dst_int16: np.ndarray,
    dst_transform: Affine,
    dst_width: int,
    dst_height: int,
    nodata_int16: np.int16,
    gridfloat_nodata: np.float32,
) -> None:
    """
    Write elevation raster data into legacy ESRI GridFloat format (.flt + .hdr).

    ``NMSIM`` uses a legacy ESRI GridFloat for elevation/impedance inputs, which has two associated files:
        - .flt: binary grid values (float32)
        - .hdr ASCII header containing grid geometry and NODATA metadata

    Inputs
    ------
    output_base : str | Path
        Output path without extension; writes ``{output_base}.flt``, ``{output_base}.hdr``.
    dst_int16 : np.ndarray
        2D array of elevations in int16.
    dst_transform : affine.Affine
        The affine transformation matrix representing reprojection from WGS84 to NAD83 GCS North American.
    dst_width : int
        Number of columns in the raster.
    dst_height : int
        Number of rows in the raster.
    nodata_int16 : np.int16
        The int16 value used as NODATA in ``dst_int16``.
    gridfloat_nodata : np.float32
        The float32 NODATA value to store in the .flt grid.

    Returns
    -------
    None

    Notes
    -----
    [1] For more information on ESRI GridFloat format, see https://www.loc.gov/preservation/digital/formats/fdd/fdd000422.shtml
    [2] Header fields use transform geometry: xllcorner, yllcorner, and cellsize.
    """
    # first we'll write the grid, .flt
    grid = dst_int16.astype(np.float32)

    grid[dst_int16 == nodata_int16] = gridfloat_nodata

    grid.tofile(Path(output_base).with_suffix(".flt"))

    # then we'll write the header, .hdr
    transform = dst_transform

    xllcorner = transform.c

    yllcorner = transform.f + dst_height * transform.e

    cellsize = transform.a

    with open(Path(output_base).with_suffix(".hdr"), "w") as hdr:
        hdr.write(f"ncols         {dst_width}\n")
        hdr.write(f"nrows         {dst_height}\n")
        hdr.write(f"xllcorner     {xllcorner:.15f}\n")
        hdr.write(f"yllcorner     {yllcorner:.15f}\n")
        hdr.write(f"cellsize      {cellsize:.15f}\n")
        hdr.write(f"NODATA_value  {gridfloat_nodata:.0f}\n")
        hdr.write("byteorder     LSBFIRST\n")
