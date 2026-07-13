from affine import Affine
import argparse
import geopandas as gpd
import glob
import numpy as np
from pathlib import Path
import rasterio
from rasterio.mask import mask as rio_mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import box, mapping, Point
import nps_active_space.utils.config as cfg
from nps_active_space.utils.computation import coords_to_utm

"""
This script creates the overarching NMSIM project that is used by every module of `NPS-ActiveSpace`.
This requires two conceptual inputs: (1) a study area polygon, and (2) a listening location point.
We use the study area to prepare a portion of a Digital Elevation Model for use in computation and visualization.
The elevation data is also converted to an antiquated ESRI GridFloat format required by `NMSIM`.
We represent the listening location coordinates as a canonical NMSIM site (.sit) file.
"""

__all__ = [
    'create_NMSIM_site_dir',
    'create_NMSIM_site_file',
    'create_study_area',
    'create_NMSIM_elevation_tif',
    'create_gridfloat'
]

def create_NMSIM_site_dir(siteDir):
    """
    Create a canonical NMSIM site directory structure expected by downstream NPS-ActiveSpace modules.

    Inputs
    ------
    siteDir (str | pathlib.Path): base directory where the canonical site structure will be created

    Returns
    -------
    None
    """
    subfolders = [
            "Input_Data",
            "Input_Data/01_ELEVATION",
            "Input_Data/02_IMPEDANCE",
            "Input_Data/03_TRAJECTORY",
            "Input_Data/04_LAYERS",
            "Input_Data/05_SITES",
            "Input_Data/06_AMBIENCE",
            "Input_Data/07_WEATHER",
            "Input_Data/08_TREES",
            "Output_Data",
            "Output_Data/ASCII",
            "Output_Data/IMAGES",
            "Output_Data/SITE",
            "Output_Data/TIG_TIS",
    ]

    siteDir = Path(siteDir)
    for folderExt in subfolders:
        Path(siteDir / folderExt).mkdir(parents=True, exist_ok=True)


def create_NMSIM_site_file(site_dir, unit, site, year, long_utm, lat_utm, height=1.5):

    """
    Create an NMSIM site (.sit) file for a given NPS microphone deployment.
    The .sit file represents a "listener location" (microphone location), one of two conceptual inputs to NPS-ActiveSpace.
    
    Inputs
    ------
    site_dir (str): a path location of a canonical NMSIM site directory
    unit (str): 4-character NPS Alpha Code, e.g. "BITH", "YUCH"
    site (str): alpha-numeric acoustic monitoring site code, e.g., "002", "TRLA"
    year (int): four digit deployment year, e.g. 2018
    long_utm (float): listening location longitude in meters (microphone's UTM zone)
    lat_utm (float): listening location latitude in meters (microphone's UTM zone)
    height (float): microphone height in meters; defaults to ANSI standard, 1.5 meters
    
    Returns
    -------
    None
    """
    
    out_path = Path(site_dir) / "Input_Data" / "05_SITES" / (unit + site + str(year) + ".sit") # deployment location can subtly vary by year
    elev_dir = Path(site_dir) / "Input_Data" / "01_ELEVATION"
    matches = list(elev_dir.glob("*.flt"))
    if not matches:
        raise FileNotFoundError(f"No .flt files found in {elev_dir}")

    # open a file and write to it
    with open(out_path, 'w') as site_file:

        site_file.write("    0\n")
        site_file.write("    1\n")
        site_file.write("{0:19.0f}.{1:9.0f}.{2:10.5f} {3:20}\n".format(long_utm, lat_utm, height, unit+site))
        site_file.write(str(next(elev_dir.glob("*.flt")))+"\n")


def create_study_area(site_dir, unit, site, year, study_area_sw_corner, study_area_ne_corner) -> gpd.GeoDataFrame:
    """
    Create the most important conceptual input of `NPS-ActiveSpace`: a user-defined study area.
    
    Builds a rectangular polygon from the provided SW/NE corners and saves it as an ESRI shapefile (.shp).

    Inputs
    ------
    site_dir (str|Path): a canonical NMSIM project directory
    unit (str): 4-character NPS Alpha Code, e.g. "BITH", "YUCH"
    site (str): alpha-numeric acoustic monitoring site code, e.g., "002", "TRLA"
    year (int): four digit deployment year, e.g. 2018
    study_area_sw_corner (tuple): study area SW corner x y (WGS84). Example: (-136.088360, 58.569310)
    study_area_ne_corner (tuple): study area NE corner x y (WGS84). Example: (-135.818994, 58.706095)

    Returns
    -------
    study_area_wgs84 (gpd.GeoDataFrame): the rectangular study area geometry (WGS84),
                                          with columns ["Unit", "Site", "Year"]
    """
    # Input args are WGS84 lon/lat; shapely box expects x1,y1,x2,y2 => lon/lat
    study_area_wgs84 = gpd.GeoDataFrame([[unit, site, year]],
                                        geometry=[box(*study_area_sw_corner, *study_area_ne_corner)],
                                        crs="EPSG:4326",
                                        columns=["Unit", "Site", "Year"])
    study_area_path = Path(site_dir) / (f"{unit}{site}_study_area.shp")
    study_area_proj = study_area_wgs84.to_crs("EPSG:4269")
    study_area_proj.to_file(study_area_path)
    return study_area_wgs84


def create_NMSIM_elevation_tif(source_dem, study_area_wgs84, dst_crs, output_tif, feet_to_meters, nodata_int16) -> tuple[np.ndarray, Affine, int, int]:
    """
    Clip a source DEM to the study area, reproject, and convert to int16 elevation suitable for `NMSIM`.
    Saves a .tif that is used downstream for computation and visualization tasks.

    Inputs
    ------
    source_dem (str | pathlib.Path): path to the input DEM raster
    study_area_wgs84 (gpd.GeoDataFrame): study area polygon (WGS84)
    dst_crs (str): coordinate reference system for the destination raster
    output_tif (str | pathlib.Path): path for output GeoTIFF file
    feet_to_meters (float): scalar constant converting feet to meters
    nodata_int16 (np.int16): the int16 value used as NODATA in the output GeoTIFF

    Returns
    -------
    dst_int16 (np.ndarray): a 2D array (dtype int16) of clipped/reprojected elevations
    dst_transform (affine.Affine): the affine transformation matrix representing reprojection from WGS84 to NAD83 GCS North American
    dst_width (int): number of columns in the output raster
    dst_height (int): number of rows in the output raster

    Notes
    -----
    [1] We assume the source DEM is in feet.
    [2] `NMSIM` requires a 16-bit signed integer GeoTIFF.
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

    # Reproject + unit conversion (assumes input raster in FEET)
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
    dst[~nodata_mask] = np.rint(dst[~nodata_mask] * feet_to_meters)
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


def create_gridfloat(output_base, dst_int16, dst_transform, dst_width,
                     dst_height, nodata_int16, gridfloat_nodata):
    """
    Write elevation raster data into legacy ESRI GridFloat format (.flt + .hdr).

    `NMSIM` uses a legacy ESRI GridFloat for elevation/impedance inputs, which has two associated files:
        - .flt: binary grid values (float32)
        - .hdr ASCII header containing grid geometry and NODATA metadata

    
    Inputs
    ------
    output_base (str | pathlib.Path): output path without extension; writes {output_base}.flt, {output_base}.hdr
    dst_int16 (np.ndarray): 2D array of elevations in int16
    dst_transform (affine.Affine): the affine transformation matrix representing reprojection from WGS84 to NAD83 GCS North American
    dst_width (int): number of columns in the raster
    dst_height (int): number of rows in the raster
    nodata_int16 (np.int16): the int16 value used as NODATA in `dst_int16`
    gridfloat_nodata (np.float32): the float32 NODATA value to store in the .flt grid
    
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

    grid.tofile(output_base.with_suffix(".flt"))

    # then we'll write the header, .hdr
    transform = dst_transform

    xllcorner = transform.c

    yllcorner = transform.f + dst_height * transform.e

    cellsize = transform.a

    with open(output_base.with_suffix(".hdr"), "w") as hdr:

        hdr.write(f"ncols         {dst_width}\n")
        hdr.write(f"nrows         {dst_height}\n")
        hdr.write(f"xllcorner     {xllcorner:.15f}\n")
        hdr.write(f"yllcorner     {yllcorner:.15f}\n")
        hdr.write(f"cellsize      {cellsize:.15f}\n")
        hdr.write(f"NODATA_value  {gridfloat_nodata:.0f}\n")
        hdr.write("byteorder     LSBFIRST\n")


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
    project_dir = Path(cfg.read('project', 'dir'))
    site_dir = Path(cfg.read('project', 'dir')) / f"{args.unit}{args.site}"

    create_NMSIM_site_dir(site_dir) # a standardized template directory structure for `NMSIM` projects
    print(f"Created a new NMSIM site directory {args.unit}{args.site} in {project_dir}.")

    study_area_wgs84 = create_study_area(site_dir=site_dir,
                                         unit=args.unit,
                                         site=args.site,
                                         year=args.year,
                                         study_area_sw_corner=args.studyarea_sw,
                                         study_area_ne_corner=args.studyarea_ne)
    study_area_path = site_dir / f"{args.unit}{args.site}_study_area.shp"
    print(f"Saved study area to {study_area_path}")

    # to create a .sit file, we will eventually need the UTM coordinates of the microphone...
    utm_epsg, _ = coords_to_utm(lat=args.mic_coord[1], lon=args.mic_coord[0]) 
    mic_utm = gpd.GeoSeries([Point(args.mic_coord)], crs="EPSG:4326").to_crs(utm_epsg)

    source_dem = cfg.read('data', 'dem')
    # quirk: `NMSIM` expects the elevation input's spatial reference to be the UTM Zone of the westernmost exent
    #         this is occasionally different than the microphone's UTM Zone in the .sit file
    _, utm_zone_str = coords_to_utm(lat=args.studyarea_sw[1], lon=args.studyarea_sw[0]) 

    output_base = site_dir / "Input_Data" / "01_ELEVATION" / f"elevation_m_nad83_utm{utm_zone_str}"
    output_tif = output_base.with_suffix(".tif")

    dst_int16, dst_transform, dst_width, dst_height = create_NMSIM_elevation_tif(source_dem=source_dem,
                                                                                 study_area_wgs84=study_area_wgs84,
                                                                                 dst_crs=dst_crs,
                                                                                 output_tif=output_tif,
                                                                                 feet_to_meters=feet_to_meters,
                                                                                 nodata_int16=nodata_int16)
    print(f"Elevation data have been written to {output_tif}.")
    # TODO when a National Landcover Database (NLCD) mapping has been construed, an equivalent `create_NMSIM_landcover_tif` function would belong here...
    
    create_gridfloat(output_base=output_base,
                     dst_int16=dst_int16,
                     dst_transform=dst_transform,
                     dst_width=dst_width,
                     dst_height=dst_height,
                     nodata_int16=nodata_int16,
                     gridfloat_nodata=gridfloat_nodata)
    print("Grid float elevation file components (.flt, .hdr) have been written.")

    create_NMSIM_site_file(site_dir,
                           unit=args.unit,
                           site=args.site,
                           year=args.year,
                           long_utm=float(mic_utm.x.iloc[0]),
                           lat_utm=float(mic_utm.y.iloc[0]),
                           height=1.5)
    print("NMSIM site file has been written. Finished!")