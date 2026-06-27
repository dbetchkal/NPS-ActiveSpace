import argparse
import geopandas as gpd
import glob
import numpy as np
import os
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import box, mapping
import nps_active_space.utils.config as cfg
from mps_active_space.utils.helpers import make, make_NMSIM_site_dir

"""
This script creates the overarching NMSIM project that is used by every module of `NPS-ActiveSpace`. 
In doing so, it clips a section of a Digital Elevation Model for use in computation and visualization.
It also saves a canonical NMSIM site (.sit) file, which represents the listening location.


"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    argparse.add_argument('-e', '--environment', required=True,
                          help="The configuration environment to run the script in.")
    argparse.add_argument('-u', '--unit', required=True,
                          help="Four letter unit code. E.g. DENA")
    argparse.add_argument('-s', '--site', required=True,
                          help="Four letter site code. E.g. TRLA")
    argparse.add_argument('-y', '--year', type=int, required=True,
                          help="Four digit year. E.g. 2018")
    argparse.add_argument('--mic-coord', type=str, required=True,
                          help="Float coordinates x,y of mic (WGS84). E.g. '-136.088360, 58.569310'")
    argparse.add_argument('--studyarea-sw', type=str, required=True,
                          help="Float coordinates x,y of study area's southwest/lower-left corner (WGS84). E.g. '-136.088360, 58.569310'")
    argparse.add_argument('--studyarea-ne', type=str, required=True,
                          help="Float coordinates x,y of study area's northeast/upper-right corner (WGS84). E.g. '-135.818994, 58.706095'")  

    args = argparse.parse_args()

    cfg.initialize(environment=args.environment)

    project_dir = f"{cfg.read('project', 'dir')}
    site_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    make_NMSIM_site_dir(project_dir) # generates a new template direcory structure


    study_area = gpd.GeoDataFrame([[u,s,y]],
                                  geometry=[box(*lower_left, *upper_right)],
                                  crs="EPSG:4326", columns=["Unit","Site","Year"])

    save_default_study_area(project_directory, study_area, unit=u, site=s)

    # ==============================================================================
    # USER SETTINGS
    # ==============================================================================

    source_dem = r"S:\Sound\NimSIM_DEMs\16_Bit\AKR_DEM.TIF"

    output_base = os.path.join(project_directory, rf"Input_Data\01_ELEVATION\elevation_m_nad83_utm{utm_zone_str}")
    output_tif = output_base + ".tif"

    feet_to_meters = 0.3048

    nodata_int16 = np.int16(32767)

    gridfloat_nodata = np.float32(-32768.0)

    dst_crs = "EPSG:4269"


    # ==============================================================================
    # READ + CLIP
    # ==============================================================================

    with rasterio.open(source_dem) as src:

        study = study_area.to_crs(src.crs)

        clipped, clipped_transform = mask(
            src,
            [mapping(g) for g in study.geometry],
            crop=True,
            filled=False,
        )

        clipped = clipped[0]

        src_profile = src.profile

        src_crs = src.crs


    # ==============================================================================
    # REPROJECT
    # ==============================================================================

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

    dst = np.full(
        (dst_height, dst_width),
        nodata_int16,
        dtype=np.float32,
    )

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

    mask = dst == nodata_int16

    dst[~mask] = np.rint(dst[~mask] * feet_to_meters)

    dst = dst.astype(np.int16)


    # ==============================================================================
    # WRITE GEOTIFF
    # ==============================================================================

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

    with rasterio.open(output_tif, "w", **profile) as ds:
        ds.write(dst, 1)


    # ==============================================================================
    # WRITE GRIDFLOAT
    # ==============================================================================

    grid = dst.astype(np.float32)

    grid[dst == nodata_int16] = gridfloat_nodata

    grid.tofile(output_base + ".flt")


    # ==============================================================================
    # WRITE HEADER
    # ==============================================================================

    transform = dst_transform

    xllcorner = transform.c

    yllcorner = transform.f + dst_height * transform.e

    cellsize = transform.a

    with open(output_base + ".hdr", "w") as hdr:

        hdr.write(f"ncols         {dst_width}\n")
        hdr.write(f"nrows         {dst_height}\n")
        hdr.write(f"xllcorner     {xllcorner:.15f}\n")
        hdr.write(f"yllcorner     {yllcorner:.15f}\n")
        hdr.write(f"cellsize      {cellsize:.15f}\n")
        hdr.write(f"NODATA_value  {gridfloat_nodata:.0f}\n")
        hdr.write("byteorder     LSBFIRST\n")


    # to create a .sit file, we need the UTM coordinates of the microphone
    utm_epsg,utm_zone_str = coords_to_utm(lat=mic_coord[1], lon=mic_coord[0])
    mic_utm = gpd.GeoSeries([Point(mic_coord)], crs="EPSG:4326").to_crs(utm_epsg)
    _,utm_zone_str = coords_to_utm(lat=studyarea[1], lon=lower_left[0])

    create_NMSIM_site_file(project_directory, unit=unit, site=site, year=year, long_utm=mic_utm.x[0], lat_utm=mic_utm.y[0], height=1.5)

    print("Finished.")    