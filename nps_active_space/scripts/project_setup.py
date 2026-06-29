import argparse
import geopandas as gpd
import glob
import numpy as np
import os
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import box, mapping, Point
import nps_active_space.utils.config as cfg
from nps_active_space.utils.computation import coords_to_utm
from nps_active_space.utils.helpers import make, make_NMSIM_site_dir, create_NMSIM_site_file

"""
This script creates the overarching NMSIM project that is used by every module of `NPS-ActiveSpace`. 
In doing so, it clips a section of a Digital Elevation Model for use in computation and visualization.
It also saves a canonical NMSIM site (.sit) file as a representation of the listening location.
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--environment', required=True,
                          help="The configuration environment to run the script in.")
    parser.add_argument('-u', '--unit', required=True,
                          help="Four letter unit code. E.g. DENA")
    parser.add_argument('-s', '--site', required=True,
                          help="Four letter site code. E.g. TRLA")
    parser.add_argument('-y', '--year', type=int, required=True,
                          help="Four digit year. E.g. 2018")
    parser.add_argument('--mic-coord', nargs=2, type=float, required=True,
                        help="Mic coordinates x y (WGS84). Example: -136.088360 58.569310")
    parser.add_argument('--studyarea-sw', nargs=2, type=float, required=True,
                        help="Study area SW corner x y (WGS84). Example: -136.088360 58.569310")
    parser.add_argument('--studyarea-ne', nargs=2, type=float, required=True,
                        help="Study area NE corner x y (WGS84). Example: -135.818994 58.706095")


    args = parser.parse_args()

    feet_to_meters = 0.3048

    nodata_int16 = np.int16(32767)

    gridfloat_nodata = np.float32(-32768.0)

    dst_crs = "EPSG:4269"

    cfg.initialize(environment=args.environment)

    project_dir = f"{cfg.read('project', 'dir')}"
    site_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    make_NMSIM_site_dir(site_dir) # generates a standardized template directory structure for `NMSIM` projects

    print(f"Created a new NMSIM site directory {args.unit}{args.site} in {project_dir}.")

    # the most important conceptual input of `NPS-ActiveSpace` is a user-defined study area.
    # TODO allow a user to pass a polygon, instead...
    study_area = gpd.GeoDataFrame([[args.unit,args.site,args.year]],
                                  geometry=[box(*args.studyarea_sw, *args.studyarea_ne)],
                                  crs="EPSG:4326", columns=["Unit","Site","Year"])
    study_area_proj = study_area.to_crs("EPSG:4269") # `NMSIM` uses NAD83 GSC North America, from U.S. Defense Mapping Agency TR8350.2 revision of August 1993
    study_area_proj.to_file(os.path.join(site_dir, args.unit+args.site+"_study_area.shp"))

    # to create a .sit file, we will eventually need the UTM coordinates of the microphone...
    utm_epsg,utm_zone_str = coords_to_utm(lat=args.mic_coord[1], lon=args.mic_coord[0])
    mic_utm = gpd.GeoSeries([Point(args.mic_coord)], crs="EPSG:4326").to_crs(utm_epsg)
    _,utm_zone_str = coords_to_utm(lat=args.studyarea_sw[1], lon=args.studyarea_sw[0])

    source_dem = cfg.read('data', 'dem')

    output_base = os.path.join(site_dir, rf"Input_Data\01_ELEVATION\elevation_m_nad83_utm{utm_zone_str}")
    output_tif = output_base + ".tif"

    # clip the study area from the source DEM
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

    # reproject raster into NAD83 GCS North American (EPSG:4269)
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

    dst[~mask] = np.rint(dst[~mask] * feet_to_meters) # we assume that the input raster is IN FEET

    dst = dst.astype(np.int16)

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

    print(f"Elevation data have been written to {output_tif}.")

    # `NMSIM` uses an antiquated "ESRI GridFloat" format for elevation and impedance data
    # https://www.loc.gov/preservation/digital/formats/fdd/fdd000422.shtml
    # it has two parts: (1) .flt, accompanied by (2) .hdr

    # first we'll write the grid, .flt
    grid = dst.astype(np.float32)

    grid[dst == nodata_int16] = gridfloat_nodata

    grid.tofile(output_base + ".flt")

    # then we'll write the header, .hdr
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

    print("Grid float elevation file components (.flt, .hdr) have been written.")    

    create_NMSIM_site_file(site_dir, 
                        unit=args.unit, 
                        site=args.site, 
                        year=args.year, 
                        long_utm=mic_utm.x[0], 
                        lat_utm=mic_utm.y[0], 
                        height=1.5)

    print("NMSIM site file has been written. Finished!")
    