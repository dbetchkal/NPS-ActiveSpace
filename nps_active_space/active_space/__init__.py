import os
import subprocess
from typing import List, Optional, Tuple, TYPE_CHECKING, Union

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal

from nps_active_space.utils import build_src_point_mesh, Microphone, NMSIM_bbox_utm

if TYPE_CHECKING:
    from nps_active_space.utils import Nvspl
    from shapely.geometry import Point


class ActiveSpace:

    def computer_f1(self):
        pass

    def plot(self):
        pass


class ActiveSpaceGenerator:
    """
    A Class that stores general active space generation information and produced individual ActiveSpace objects
    based on customizable parameters.

    Parameters
    ----------
    NMSIM : str
        Absolute path to the NMSIM executable to be used to generate active spaces.
    study_area : gpd.GeoDataFrame
        A gpd.GeoDataFrame of polygon(s) that make up the study area.
    root_dir : str
        Absolute path to a directory where all generated files required for running NMSIM can be stored.
    dem_src : str
        Path to a DEM raster file to be used as NMSIM input.
    ambience_src : Nvspl or str
        an NVSPL object to calculate ambience from or the absolute file path to a raster of ambience.
    quantile : int, default 50
        If using Nvspl data as the ambience source, this quantile of the data will be used to calculate the ambience.
    broadband : bool, default False
        If True and using Nvspl data as the ambience source, quantiles will be calculated from the dBA column
        instead of the 1/3rd octave band columns.
    """
    def __init__(self, NMSIM: str, study_area: gpd.GeoDataFrame, root_dir: str, dem_src: str,
                 ambience_src: Union['Nvspl', str], quantile: int = 50, broadband: bool = False):

        self.root_dir = root_dir
        self.study_area = study_area.to_crs('epsg:4269')  # Convert study area to NAD83 (NMSIM preference)
        self.broadband = broadband
        self.quantile = quantile
        self.NMSIM = NMSIM

        self._make_dir_tree()
        self.dem_tif, self.dem_flt = self._create_dem_files(dem_src, self.study_area)
        # self.ambience = self._set_ambience(ambience_src) # TODO

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

    def _create_dem_files(self, dem_src: str, study_area: gpd.GeoDataFrame, project: bool = True,
                          suffix: str = '') -> Tuple[str, str]:
        """
        Project and mask a DEM raster to be used as input into NMSIM.
        Follows this tutorial: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html

        Parameters
        ----------
        dem_src : str
            Absolute path to a DEM raster to be projected and clipped for use by NMSIM.
        study_area : gpd.GeoDataFrame
            The study area to clip the DEM raster to.
        suffix : str, default ''
        project : bool, default True
           True to project the DEM file to NAD83 before clipping it.

        Returns
        -------
        DEM tif filename for the study area, DEM flt filename for the study area
        """
        # ------ Mask the DEM File ------- #

        # Project the DEM file to NAD83.
        if project:
            dem_projected_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation_NAD83{suffix}.tif"
            gdal.Warp(dem_projected_filename, dem_src, dstSRS='EPSG:4269')
            dem_src = dem_projected_filename

        # Output the study area, in the proper projection, to a shapefile so it can be used for masking.
        study_area_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/study_area_NAD83{suffix}.shp"
        study_area.to_file(study_area_filename)

        # Mask the DEM file with the study area shapefile.
        dem_masked_filename = f"{self.root_dir}\\Input_Data\\01_ELEVATION\\elevation_NAD83_mask{suffix}.tif"
        gdal.Warp(
            dem_masked_filename,
            dem_src,
            cutlineDSName=study_area_filename,
            cropToCutline=True
        )

        # ------ Create .flt DEM File & Header ------- #

        flt_filename = f"{self.root_dir}\\Input_Data\\01_ELEVATION\\elevation_NAD83_mask{suffix}.flt"
        flt_header_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation_NAD83_mask{suffix}.hdr"

        gdal.Translate(flt_filename, dem_masked_filename, options="-ot Float32 -of ehdr -a_nodata -9999")

        # the header file doesn't write correctly... manually overwrite this:
        old_hdr = pd.read_csv(flt_header_filename, header=None, delim_whitespace=True, index_col=0).T

        # compute new lower left corner y-val
        yllcorner = float(old_hdr.ULYMAP) - float(old_hdr.NROWS) * float(old_hdr.XDIM)

        # write a new header file exactly as output by ESRI
        with open(flt_header_filename, 'w') as header:
            header.write("{:14}{:}\n".format("ncols", old_hdr.NCOLS.values[0]))
            header.write("{:14}{:}\n".format("nrows", old_hdr.NROWS.values[0]))
            header.write("{:14}{:}\n".format("xllcorner", old_hdr.ULXMAP.values[0]))
            header.write("{:14}{:}\n".format("yllcorner", yllcorner))
            header.write("{:14}{:}\n".format("cellsize", old_hdr.XDIM.values[0]))
            header.write("{:14}{:}\n".format("NODATA_value", old_hdr.NODATA.values[0]))
            header.write("{:14}{:}".format("byteorder", "LSBFIRST"))

        return dem_masked_filename, flt_filename

    def _create_trajectory_file(self, points: List['Point'], crs: str, filename: str) -> str:
        """
        Create a trajectory file from a list of points.

        Parameters
        ----------
        points : List of shapely Points
            Trajectory points to be written and run through NMSIM.
        crs : str
            The crs the Points should be in before being written to the trajectory file.
        filename : str
            Name of the trajectory file.

        Returns
        -------
        The trajectory file name.
        """
        trajectory_filename = f"{self.root_dir}\\Input_Data\\03_TRAJECTORY\\{filename}.trj"

        trajectory = gpd.GeoDataFrame([], geometry=points, crs=crs)
        trajectory['heading'] = np.random.choice(range(0, 360), size=len(points), replace=True)
        trajectory['climb_angle'] = 0
        trajectory['power'] = 95
        trajectory['rol'] = 0

        velocity = 70   # m/s
        trajectory["knots"] = 1.94384 * velocity    # Convert airspeed to knots

        dist = np.diff(np.array([[x, y, z] for x, y, z in zip(trajectory.geometry.x, trajectory.geometry.y, trajectory.geometry.z)]), axis=0)
        time_elapsed = np.cumsum(np.array([np.linalg.norm(d) for d in dist]) / velocity)
        time_elapsed = np.append(time_elapsed, np.nan)  # last row must be NaN because time is based on differencing
        trajectory["time_elapsed"] = time_elapsed

        with open(trajectory_filename, 'w') as trajectory_file:
            trajectory_file.write("Flight track trajectory variable description:\n")
            trajectory_file.write(" time - time in seconds from the reference time\n")
            trajectory_file.write(" Xpos - x coordinate (UTM)\n")
            trajectory_file.write(" Ypos - y coordinate (UTM)\n")
            trajectory_file.write(" UTM Zone  " + crs + "\n")
            trajectory_file.write(" Zpos - z coordinate in meters MSL\n")
            trajectory_file.write(" heading - aircraft compass bearing in degrees\n")
            trajectory_file.write(" climbANG - aircraft climb angle in degrees\n")
            trajectory_file.write(" vel - aircraft velocity in knots\n")
            trajectory_file.write(" power - % engine power\n")
            trajectory_file.write(" roll - bank angle (right wing down), degrees\n")
            trajectory_file.write("FLIGHT " + filename + "\n")
            trajectory_file.write("TEMP.  59.0\n")
            trajectory_file.write("Humid.  70.0\n")
            trajectory_file.write("\n")
            trajectory_file.write(
                "         time(s)        Xpos           Ypos           Zpos         heading        climbANG       Vel            power          rol\n")

            for ind, point in trajectory.iterrows():
                trajectory_file.write(
                    "{0:15.3f}".format(point.time_elapsed) +
                    "{0:15.3f}".format(point.geometry.x) +
                    "{0:15.3f}".format(point.geometry.y) +
                    "{0:15.3f}".format(point.geometry.z) +
                    "{0:15.3f}".format(point.heading) +
                    "{0:15.3f}".format(point.climb_angle) +
                    "{0:15.3f}".format(point.knots) +
                    "{0:15.3f}".format(point.power) +
                    "{0:15.3f}".format(point.rol) + "\n"
                )

        return trajectory_filename

    def _create_site_file(self, mic: Microphone, dem_file: str) -> str:
        """
        Create a required NMSIM site file representing the location we are testing audibility for.

        Parameters
        ----------
        mic : Microphone
            A Microphone object to create a NMSIM site file for.

        Returns
        -------
        The name of the site file.
        """
        site_filename = f"{self.root_dir}\\Input_Data\\05_SITES\\{mic.name}.sit"
        with open(site_filename, 'w') as site_file:
            site_file.write("    0\n")
            site_file.write("    1\n")
            site_file.write("{0:19.0f}.{1:9.0f}.{2:10.5f} {3:20}\n".format(mic.x, mic.y, mic.z, mic.name))
            site_file.write(f"{dem_file}\n")

        return site_filename

    def _create_instruction_files(self, dem_file: str, site_file: str, trajectory_file: str,
                                  omni_source_file: str) -> str:
        """
        Create the batch.txt and control.nms instructions files needed to run NMSIM.

        Parameters
        ----------
        dem_file : str
            Absolute path to elevation flt file.
        site_file : str
            Absolute path to site file.
        trajectory_file : str
            Absolute path to trajectory file.
        omni_source_file : str
            Absolute path to omni source file.

        Returns
        -------
        Batch file name.
        """
        control_file = f"{self.root_dir}\\control.nms"
        batch_file = f"{self.root_dir}\\batch.txt"
        tis_directory = f"{self.root_dir}\\Output_Data\\TIG_TIS"

        with open(control_file, 'w') as nms:
            nms.write(dem_file + "\n")
            nms.write("-\n")
            nms.write(site_file + "\n")
            nms.write(trajectory_file + "\n")
            nms.write("-\n")
            nms.write("-\n")
            nms.write(omni_source_file + "\n")
            nms.write("{0:11.4f}   \n".format(500.0000))
            nms.write("-\n")
            nms.write("-")

        # write the batch file to create a site-based analysis
        with open(batch_file, 'w') as batch:
            batch.write("open\n")
            batch.write(control_file + "\n")
            batch.write("site\n")
            batch.write(f"{tis_directory}\\{os.path.basename(trajectory_file)[:-4]}" + "\n")
            batch.write("dbf: no\n")
            batch.write("hrs: 0\n")
            batch.write("min: 0\n")
            batch.write("sec: 0.0")

        return batch_file

    def _run_NMSIM(self, batch_file: str):

        process = subprocess.Popen([self.NMSIM, batch_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = process.communicate()
        print(stdout)
        print(stderr)

        # output_messages = stdout.decode("utf-8").split("\r\n")
        # output_messages = [out for out in output_messages if out.strip() != ''])

        # # slightly messier printing for error messages
        # if (stderr != None):
        #     for s in stderr.decode("utf-8").split("\r\n"):
        #         print(s.strip())

    def generate(self, omni_source: str, altitude_m: int = 3658, n_iter: int = 3,
                 simplify: int = 100, src_pt_density: int = 48, mic: Optional[Microphone] = None,
                 mesh: bool = False, mesh_density: int = 10):
        """
        # TODO

        Parameters
        ----------
        omni_source : str
        altitude_m : int, default 3658 meters (equivalent to 12000 ft)
        n_iter : int
        simplify : int
        mesh : bool
        mesh_density : int
        src_pt_density : int
        mic : Microphone, default None
        """
        study_areas = None if mesh else [self.study_area]   # TODO: mesh

        for i, study_area in enumerate(study_areas):

            # Determine the UTM CRS on the western-most edge of the initial study area.
            crs = NMSIM_bbox_utm(study_area)

            # Set the active space to initially be the same as the study area. We will refine it.
            active_space = study_area.to_crs(crs).geometry.iloc[0]

            # If not creating a mesh or if no specific microphone location is passed in, create a microphone at the
            #  center point of the study area.
            if mic and not mesh:
                mic.to_crs(crs, inplace=True)
            else:
                mic = Microphone(
                    name=f"centroid_{i}",
                    lat=active_space.centroid.x,
                    lon=active_space.centroid.y,
                    z=1.60,  # m, average height of human ear
                    crs=crs
                )

            # TODO: does this need to happen every time...?
            dem_filename, flt_filename = self._create_dem_files(self.dem_tif, study_area=study_area, project=False,
                                                    suffix=f'_{mic.name}') if mesh else (self.dem_tif, self.dem_flt)
            site_filename = self._create_site_file(mic, flt_filename)

            for j in range(n_iter):

                source_pt_mesh = build_src_point_mesh(active_space, src_pt_density, altitude_m)
                trajectory_filename = self._create_trajectory_file(source_pt_mesh, crs, f"{mic.name}_{i}")

                batch_file = self._create_instruction_files(flt_filename, site_filename, trajectory_filename, omni_source)

                self._run_NMSIM(batch_file)


                # self._find_audible() # TODO

                # TODO: shrink the study size
