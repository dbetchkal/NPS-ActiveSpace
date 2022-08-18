import os
from typing import TYPE_CHECKING, Union

import geopandas as gpd
import rasterio
import rasterio.mask
from osgeo import gdal

from nps_active_space.utils import NMSIM_bbox_utm

if TYPE_CHECKING:
    from nps_active_space.utils import Microphone, Nvspl


class ActiveSpace:

    def computer_f1(self):
        pass

    def plot(self):
        pass


class ActiveSpaceGenerator:
    """
    # TODO

    Parameters
    ----------
    study_area : gpd.GeoDataFrame
        A gpd.GeoDataFrame of polygon(s) that make up the study area.
    dem_src : str
        Path to a DEM raster file to be cropped.
    ambience_src : Nvspl or str
        an NVSPL with third-octave band ambience or the absolute file path to a raster of ambience.
    broadband : bool, default False
    quantile : int, default 50
    """
    def __init__(self, study_area: gpd.GeoDataFrame, root_dir: str, dem_src: str,
                 ambience_src: Union['Nvspl', str],
                 broadband: bool = False, quantile: int = 50,):

        self.root_dir = root_dir
        self.study_area = study_area.to_crs('epsg:4269')  # Convert study area to NAD83 (NMSIM preference)
        self.broadband = broadband
        self.quantile = quantile

        self._make_dir_tree()
        # self.ambience = self._set_ambience(ambience_src) # TODO
        self.dem_file = self._write_dem(dem_src)

        # TODO: raster to polygon...?

    def _set_ambience(self, ambience_src: Union['Nvspl', str]):
        """
        Load ambience

        Parameters
        ----------
        ambience_src : Nvspl or str

        """
        # TODO; not sure how this gets used yet...

        if type(ambience_src) == Nvspl:
            if self.broadband:
                Lx = ambience_src.loc[:, "12.5":"20000"].quantile(1 - (self.quantile / 100))
            else:
                Lx = ambience_src.loc[:, 'dbA'].quantile(1 - (self.quantile / 100))

        # else: # tODO
        #     # If the ambience source is an ambience raster.
        #     with rasterio.open(ambience_src) as raster:
        #         if raster.crs != self.crs:
        #             mic = self.mic
        #             mic.to_crs(raster.crs)
        #         band1 = raster.read(1)  # Look up the broadband level at that location
        #         Lx = band1[raster.index(mic.x, mic.y)]

        return Lx

    def _make_dir_tree(self):
        """Create a canonical NMSIM project directory. Copied from NMSIM_Create_Base_Layers.py"""
        directories = [
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
            "Output_Data/TIG_TIS"
        ]

        for directory in directories:
            if not os.path.exists(f"{self.root_dir}/{directory}"):
                os.makedirs(f"{self.root_dir}/{directory}")

    def _write_dem_files(self, dem_src: str) -> str:
        """
        Project and mask a DEM raster to be used as input into NMSIM.
        Follows this tutorial: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html

        Parameters
        ----------
        dem_src : str
            Absolute path to a DEM raster to be projected and clipped for use by NMSIM.

        Returns
        -------
        dem_out : str
            Absolute path to masked DEM raster.
        """
        # Project and mask the DEM file to the study area.
        dem_projected_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation_NAD83.tif"
        gdal.Warp(dem_projected_filename, dem_src, dstSRS='EPSG:4269')

        study_area_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/study_area_NAD83.shp"
        self.study_area.to_file(study_area_filename)

        dem_masked_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation_NAD83_mask.tif"
        gdal.Warp(
            dem_masked_filename,
            dem_src,
            cutlineDSName=study_area_filename,
            cropToCutline=True
        )

        # Create a parallel .flt DEM file.
        flt_filename = dem_masked_filename.replace('.tif', '.flt')

        # gdal_translate.exe command line args; `ot` output type, `of` output format `a_nodata` no-data filler
        base_cmd = "-ot Float32 -of ehdr -a_nodata -9999"
        f"{} {base_cmd} {dem_masked_filename} {flt_filename}"

        def quote(item):
            # this avoids any issues with file paths including white space as command line arguments
            return "\"" + item + "\""

        fullCmd = ' '.join([gdal_translate_path, cmd, quote(src), quote(dst)])
        subprocess.call(fullCmd)  # call GDAL

        # the header file doesn't write correctly... manually overwrite this:
        hdr = os.path.splitext(src)[0] + '.hdr'  # removes the .tif file extension and appends '.hdr'
        old_hdr = pd.read_csv(hdr, header=None, delim_whitespace=True, index_col=0).T

        # compute new lower left corner y-val
        yllcorner = float(old_hdr.ULYMAP) - float(old_hdr.NROWS) * float(old_hdr.XDIM)

        # write a new header file exactly as output by ESRI
        with open(hdr, 'w') as header:
            header.write("{:14}{:}\n".format("ncols", old_hdr.NCOLS.values[0]))
            header.write("{:14}{:}\n".format("nrows", old_hdr.NROWS.values[0]))
            header.write("{:14}{:}\n".format("xllcorner", old_hdr.ULXMAP.values[0]))
            header.write("{:14}{:}\n".format("yllcorner", yllcorner))
            header.write("{:14}{:}\n".format("cellsize", old_hdr.XDIM.values[0]))
            header.write("{:14}{:}\n".format("NODATA_value", old_hdr.NODATA.values[0]))
            header.write("{:14}{:}".format("byteorder", "LSBFIRST"))

        return dem_out

    def create_flt(self):



    def generate(self, omni_source: str, mesh: bool = False, altitude_m: int = 3658, n_tracks: int = 1,
                 density: int = 48, n_itr: int = 3, simplify: int = 100):
        """

        Parameters
        ----------

        altitude_m : int, default 3658 meters (equivalent to 12000 ft)
        n_tracks : int
        """

        # self.crs = NMSIM_bbox_utm(self.study_area)

    # def generate(self, omni_source, altitude_m=12000, density: int = 48, n_iter: int = 3,
    #              simplify: int = 100):
    #     """
    #
    #     Parameters
    #     ----------
    #     n_iter : int, default 3
    #         Number of iterations to run.
    #     """
    #     # TODO: altitude
    #     velocity = 70  # m/s, this doesn't matter actually (?)
    #
    #     active_space = None
    #     total_space = None
    #     xy = None
    #
    #     for i in range(n_iter ):
    #
    #
    #         data = gdal.Open(DEM_tif, GA_ReadOnly)
    #
    #         # extract the (projected) geometric transformation of the NAD83 raster
    #         # glean each of the corners
    #         geoTransform = data.GetGeoTransform()
    #         minx = geoTransform[0]
    #         maxy = geoTransform[3]
    #         maxx = minx + geoTransform[1] * data.RasterXSize
    #         miny = maxy + geoTransform[5] * data.RasterYSize
    #
    #
    #         transformer = pyproj.Transformer.from_crs('epsg:4269', UTM_zone, always_xy=True)
    #         out_points = transformer.transform([minx, maxx], [miny, maxy])
    #
    #         return out_points, UTM_zone
    #
    #
    #         full_extent, utm_zone = get_extent(self.project_dir)
    #
    #         # if this isn't the first iteration, active_space exists
    #         if active_space is not None:
    #             xmax, ymax, _ = active_space.max().astype(float)
    #             xmin, ymin, _ = active_space.min().astype(float)
    #             xpad = 0.2 * (xmax - xmin)  # pad extents by 20% on each side, 40% total
    #             ypad = 0.2 * (ymax - ymin)
    #
    #    pass
