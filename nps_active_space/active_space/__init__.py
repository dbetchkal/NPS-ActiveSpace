import glob
import multiprocessing as mp
import os
import subprocess
from functools import partial
from typing import Iterable, List, Optional, Tuple, Union
from uuid import uuid4

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from osgeo import gdal
from shapely.geometry import Point, Polygon
from shapely.validation import make_valid
from tqdm import tqdm

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.utils import (
    ambience_from_nvspl,
    ambience_from_raster,
    build_src_point_mesh,
    coords_to_utm,
    create_overlapping_mesh,
    Microphone,
    NMSIM_bbox_utm,
    Nvspl,
    project_raster
)


class ActiveSpaceGenerator:
    """
    A class that stores active space generation logic and produces individual active spaces as well as active space
    meshes based on customizable parameters.

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
    quantile : int, default 90
        If using Nvspl data as the ambience source, this quantile of the data will be used to calculate the ambience.
        90 = L90 = 10th quantile
    broadband : bool, default False
        If True and using Nvspl data as the ambience source, quantiles will be calculated from the dBA column
        instead of the 1/3rd octave band columns.
    """
    def __init__(self, NMSIM: str, study_area: gpd.GeoDataFrame, root_dir: str, dem_src: str,
                 ambience_src: Union['Nvspl', str], quantile: int = 90, broadband: bool = False):

        self.study_area = study_area.to_crs('epsg:4269')
        self.root_dir = root_dir
        self.broadband = broadband
        self.quantile = quantile
        self.ambience_src = ambience_src
        self.dem_src = dem_src
        self.NMSIM = NMSIM

        self._dem_file = None
        self._flt_file = None
        self._site_file = None

        self._make_dir_tree()

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

    def _mask_dem_file(self, dem_src: str, study_area: gpd.GeoDataFrame, project: bool = False,
                       buffer: Optional[int] = None, suffix: str = '') -> str:
        """
        Project and mask a DEM .tif raster.
        Follows this tutorial: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html

        Parameters
        ----------
        dem_src : str
            Absolute path to a DEM raster to be projected and clipped for use by NMSIM.
        study_area : gpd.GeoDataFrame
            The study area to clip the DEM raster to.
        suffix : str, default ''
            A suffix to add to the end of the filename to distinguish it from others.
        project : bool, default True
           True to project the DEM file to NAD83 before clipping it.

        Returns
        -------
        The absolute path to the .tif file of the projected and masked raster.
        """
        # Project the DEM file to NAD83 which is also what the study area is in.
        if project:
            dem_projected_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation{suffix}.tif"
            project_raster(dem_src, dem_projected_filename, study_area.crs)
            dem_src = dem_projected_filename

        # Output the study area, in the proper projection, to a shapefile so it can be used for masking.
        study_area_filename_prefix = f"{self.root_dir}/Input_Data/01_ELEVATION/study_area{suffix}_{uuid4()}"
        study_area_filename = f"{study_area_filename_prefix}.shp"
        if buffer:
            equal_area_crs = coords_to_utm(study_area.centroid.iat[0].y, study_area.centroid.iat[0].x)
            study_area_m = study_area.to_crs(equal_area_crs)
            study_area_m = study_area_m.buffer(buffer*1000)
            study_area = study_area_m.to_crs(study_area.crs)

        study_area.index = pd.Index(study_area.index) #, dtype=np.int64)
        study_area.to_file(study_area_filename)

        # Mask the DEM file with the study area shapefile.
        dem_masked_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation_mask{suffix}.tif"
        gdal.Warp(
            dem_masked_filename,
            dem_src,
            cutlineDSName=study_area_filename,
            cropToCutline=True
        )

        # Remove the temporary shapefile (and related files) since they were only needed for masking.
        for filename in glob.glob(f"{study_area_filename_prefix}*"):
            os.remove(filename)

        return dem_masked_filename

    def _create_dem_flt(self, dem_file: str) -> str:
        """
        Convert the DEM .tif to a DEM .flt as input into NMSIM.

        Parameters
        ----------
        dem_file : str
            Absolute path to the DEM .tif file to convert to a .flt file.

        Returns
        -------
        DEM flt filename for the study area.
        """
        flt_filename = dem_file.replace('.tif', '.flt')
        flt_header_filename = dem_file.replace('.tif', '.hdr')

        gdal.Translate(flt_filename, dem_file, options="-ot Float32 -of ehdr -a_nodata -9999")

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

        return flt_filename

    def _create_trajectory_file(self, points: List['Point'], crs: str, filename: str,
                                heading: Optional[int] = None) -> str:
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
        heading : int, default None
            The heading (yaw) to use for all points in the trajectory file. If None, a random heading will be used
            for each point.

        Returns
        -------
        The trajectory file name.
        """
        trajectory_filename = f"{self.root_dir}/Input_Data/03_TRAJECTORY/{filename}.trj"

        trajectory = gpd.GeoDataFrame([], geometry=points, crs=crs)
        trajectory['heading'] = heading if heading is not None else np.random.choice(range(0, 360), size=len(points), replace=True)
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
        dem_file : str
            Absolute path to the .flt DEM file for the study area.

        Returns
        -------
        The name of the site file.
        """
        site_filename = f"{self.root_dir}/Input_Data/05_SITES/{mic.name}.sit"
        with open(site_filename, 'w') as site_file:
            site_file.write("    0\n")
            site_file.write("    1\n")
            site_file.write("{0:19.0f}.{1:9.0f}.{2:10.5f} {3:20}\n".format(mic.x, mic.y, mic.z, mic.name))
            site_file.write(f"{dem_file}\n")

        return site_filename

    def _create_instruction_files(self, flt_file: str, site_file: str, trajectory_file: str,
                                  omni_source_file: str) -> str:
        """
        Create the batch.txt and control.nms instructions files needed to run NMSIM.

        Parameters
        ----------
        flt_file : str
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
        control_file = f"{self.root_dir}/control_{os.path.basename(trajectory_file).replace('.trj', '')}.nms"
        batch_file = f"{self.root_dir}/batch_{os.path.basename(trajectory_file).replace('.trj', '')}.txt"
        tis_directory = f"{self.root_dir}/Output_Data/TIG_TIS"

        with open(control_file, 'w') as nms:
            nms.write(flt_file + "\n")
            nms.write("-\n")
            nms.write(site_file + "\n")
            nms.write(trajectory_file + "\n")
            nms.write(f"{ACTIVE_SPACE_DIR}/data/default.wea" + "\n")
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
            batch.write(f"{tis_directory}/{os.path.basename(trajectory_file)[:-4]}" + "\n")
            batch.write("dbf: no\n")
            batch.write("hrs: 0\n")
            batch.write("min: 0\n")
            batch.write("sec: 0.0")

        return batch_file

    def _find_audible_points(self, trajectory_file: str, tis_file: str, crs: str,
                             ambience: Union[int, Iterable[int]]) -> gpd.GeoDataFrame:
        """
        Determine which points from a trajectory file are audible given the corresponding NMSIM tis output.

        Parameters
        ----------
        trajectory_file : str
            Absolute path to the trajectory file.
        tis_file : str
            Absolute path to the corresponding tis file.
        crs : str
            crs of the trajectory file and of the output GeoDataFrame. In the format 'epsg:XXXX'
        ambience : int or Iterable[int]
            The ambience level(s) at the microphone site.

        Returns
        -------
        total_space ; gpd.GeoDataFrame
            A GeoDataFrame of the NMSIM tested points from the trajectory file with an 'audible' column set to 0 or 1
            depending on the point's audibility.
        """
        # Read in the trajectory input file.
        trajectory_df = pd.read_fwf(trajectory_file, header=14, widths=[16, 14] + [15]*7)
        # trajectory_df = pd.read_fwf(trajectory_file, header=15, widths=[16, 14] + [15]*7)

        # Read in the tis NMSIM output file.
        tis_df = pd.read_fwf(tis_file, header=15, skipfooter=1, widths=[4, 12] + [5]*35)
        tis_df.drop(['SP#', 'F', '42'], axis=1, inplace=True)  # Drop the sensitivity index column
        tis_df.drop(tis_df.head(2).index, inplace=True)        # Drop dividing *****, and unneeded AMBIENT row
        tis_df.reset_index(drop=True, inplace=True)
        tis_df.columns = ["TIME", "A", "10", "12.5", "15.8", "20", "25", "31.5", "40", "50", "63",
                          "80", "100", "125", "160", "200", "250", "315", "400", "500", "630", "800", "1000",
                          "1250", "1600", "2000", "2500", "3150", "4000", "5000", "6300", "8000", "10000", "12500"]
        tis_df.loc[:, 'A':'12500'] = tis_df.loc[:, 'A':'12500'].astype(float) * 0.1  # centibels (cB) to decibels (dB)

        # Check to see if any of the frequency bands are louder than the ambient levels.
        if not self.broadband and type(self.ambience_src) == Nvspl:
            audible_times = (tis_df.loc[:, "12.5":"12500"] > ambience["12.5":"12500"].values).sum(axis=1)
        else:
            audible_times = tis_df.loc[:, "A"] > ambience

        # Determine which aircraft locations produce or do not produce audible sounds.
        audible_pts = (trajectory_df.loc[tis_df[audible_times > 0].index, ["Xpos", "Ypos", "Zpos"]])
        inaudible_pts = (trajectory_df.loc[tis_df[audible_times == 0].index, ["Xpos", "Ypos", "Zpos"]])
        audible_pts["audible"] = 1
        inaudible_pts["audible"] = 0

        active_space = gpd.GeoDataFrame(
            audible_pts,
            crs=crs,
            geometry=gpd.points_from_xy(x=audible_pts.Xpos, y=audible_pts.Ypos, z=audible_pts.Zpos)
        )
        null_space = gpd.GeoDataFrame(
            inaudible_pts,
            crs=crs,
            geometry=gpd.points_from_xy(x=inaudible_pts.Xpos, y=inaudible_pts.Ypos, z=inaudible_pts.Zpos)
        )

        total_space = pd.concat([active_space, null_space], ignore_index=True)
 
        return total_space

    @staticmethod
    def _contour_active_space(total_space: gpd.GeoDataFrame, altitude: int, max_pts: int = 5184) -> List[Point]:
        """
        Use triangulation to select points along the audible/inaudible line of the active space to more precisely
        define the boundaries.

        Parameters
        ----------
        total_space : gpd.GeoDataFrame
            All points that have been tested for audibility so far -- both audible and inaudible.
        altitude : int
            The altitude (in meters) we are calculating the active space at.
        max_pts : int, default 5184
            The number of points to use for triangulation.

        Returns
        -------
        List of Points to pass into NMSIM along the audible/inaudible line help define the active space edge.
        """
        levels = np.linspace(0, 1, 10, endpoint=False)  # these are the contour levels to calculate
        fig, ax = plt.subplots()

        # Uses Delaunay triangulation
        tri = mpl.tri.Triangulation(total_space.geometry.x.tolist(), total_space.geometry.y.tolist())
        cs = ax.tricontour(tri, total_space.audible.tolist(), levels=levels)  # contour with arbitrary point cloud
        plt.close('all')

        # append all contour points to a new xy file
        xyz = np.array([])
        for contour_level in cs.collections:  # iterate through each contour interval
            for contour_path in contour_level.get_paths():  # for each interval, iterate through each path
                x = contour_path.vertices[:, 0]
                y = contour_path.vertices[:, 1]
                z = [altitude] * len(x)
                xyz = np.array([x, y, z]).T if len(xyz) == 0 else np.append(xyz, np.array([x, y, z]).T, axis=0)

        # if we have more contour points than NMSim can handle, down-sample randomly
        if xyz.shape[0] > max_pts:
            random_inds = np.random.randint(0, xyz.shape[0], max_pts)
            xyz = xyz[random_inds, :]

        return list(map(Point, xyz))

    def _run_nmsim(self, job_name: str, source_pts: List[Point], crs: str, flt_file: str, site_file: str,
                   omni_source: str, ambience, heading: Optional[int] = None) -> gpd.GeoDataFrame:
        """
        Execute a single NMSIM job.

        Parameters
        ----------
        job_name : str
            Name of this NMSIM run to use a suffix to input and output files.
        source_pts : List[Point]
            List of shapely Points to test the audibility of in NMSIM.
        crs : str
            The crs of the points. Of the format 'epsg:XXXX..'
        flt_file : str
            Absolute path to the .flt DEM file required to run NMSIM.
        site_file : str
            Absolute path to the .sit file of the receiver point required to run NMSIM.
        heading : int, default None
            The heading (yaw) to use for all points in the trajectory file. If None, a random heading will be used
            for each point.

        Returns
        -------
        new_audibility_pts : gpd.GeoDataFrame
            A GeoDataFrame of points tested during the NMSIM run and their audibility.
        """
        trajectory_filename = self._create_trajectory_file(source_pts, crs, job_name, heading)
        batch_file = self._create_instruction_files(flt_file, site_file, trajectory_filename, omni_source)

        # Run NMSIM.
        process = subprocess.Popen([self.NMSIM, batch_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = process.communicate()

        if stderr:
            for s in stderr.decode("utf-8").split("\r\n"):
                print(s.strip())

        # Determine the audibility of points that were tested during the NMSIM run.
        new_audibility_pts = self._find_audible_points(
            trajectory_filename,
            f"{self.root_dir}/Output_Data/TIG_TIS/{job_name}.tis",
            crs,
            ambience
        )

        return new_audibility_pts

    @staticmethod
    def _build_active_space(total_space: gpd.GeoDataFrame, crs: str) -> gpd.GeoDataFrame:
        """
        Build the final active space polygon given the audibility of all tested points.

        Parameters
        ----------
        total_space : gpd.GeoDataFrame
            A GeoDataFrame of all tested points and their audibility to be used to build the final active space
            polygon from.
        crs : str
            The crs of the points. Of the format 'epsg:XXXX...'

        Returns
        -------
        A GeoDataFrame of the active space.
        """
        fig, ax = plt.subplots()  # need an axis to call tricontour function
        levels = np.linspace(0, 1, 10, endpoint=False)

        # create the triangulated irregular network and contour lines
        tri = mpl.tri.Triangulation(total_space.geometry.x.tolist(), total_space.geometry.y.tolist())
        cs = ax.tricontour(tri, total_space.audible.tolist(), levels=levels)

        # pick some contour level to choose as the active space boundary... 0.5 is somewhat arbitrary
        level_ind = np.where(cs.levels == 0.5)[0][0]
        plt.close('all')  # close triangulation figure

        # iterate through all contour paths in the line collection at level_ind
        active_space_poly = None
        for i, contour_path in enumerate(cs.collections[level_ind].get_paths()):
            x = contour_path.vertices[:, 0]
            y = contour_path.vertices[:, 1]
            new_poly = make_valid(Polygon([(i[0], i[1]) for i in zip(x, y)]))

            # Don't bother with polygons that are smaller than .5 km^2.
            if new_poly.area <= 50000:
                continue
            elif active_space_poly is None:
                active_space_poly = new_poly
            else:
                active_space_poly = active_space_poly.symmetric_difference(new_poly)

        active_space_polys_gdf = gpd.GeoDataFrame(data={'geometry': [active_space_poly]}, geometry='geometry', crs=crs)

        return active_space_polys_gdf

    def _generate(self, study_area: gpd.GeoDataFrame, dem_file: str, omni_source: str, name: str = '',
                  mic: Optional[Microphone] = None, project_dem: bool = True, altitude_m: int = 3658,
                  heading: Optional[int] = None, src_pt_density: int = 48, n_contour: int = 1) -> gpd.GeoDataFrame:
        """
        The main active space generating function. It has been separated from the other generate functions to allow
        for multiprocessing when creating an active space mesh.
        """
        crs = NMSIM_bbox_utm(study_area)  # Determine the UTM CRS on the western-most edge of the study area

        # Initialize a GeoDataFrame of source points that have gone through NMSIM and a GeoDataFrame of the
        #  current active space. The active space will initially be the same as the study area, but will be refined.
        tested_space = gpd.GeoDataFrame(columns=['audible', 'geometry'], geometry='geometry', crs=crs)
        active_space = study_area.to_crs(crs)
        study_area_extent = ([active_space.total_bounds[0], active_space.total_bounds[2]],  # ([minx, maxx],
                             [active_space.total_bounds[1], active_space.total_bounds[3]])  # [miny, maxy])

        if mic:
            mic.to_crs(crs, inplace=True)
        else:
            # Create a microphone at the center point of the study area.
            study_area_wgs84 = study_area.to_crs('epsg:4326')
            mic = Microphone(
                name=f"centroid{name}",
                lat=study_area_wgs84.centroid.iat[0].y,
                lon=study_area_wgs84.centroid.iat[0].x,
                z=1.60,  # m, average height of human ear
                crs=crs
            )

        # If no dem file has been set, create the DEM file now. Also create the flt and site files needed by NMSIM.
        dem_filename = self._dem_file or self._mask_dem_file(dem_file, study_area=study_area, project=project_dem, suffix=f'_{mic.name}')
        flt_filename = self._flt_file or self._create_dem_flt(dem_filename)
        site_filename = self._site_file or self._create_site_file(mic, flt_filename)

        # NOTE: These lines were written with the assumption that using Nvspl data for ambience levels would only
        #  really happen when not running a mesh.
        if type(self.ambience_src) == Nvspl:
            ambience = ambience_from_nvspl(self.ambience_src, self.quantile, self.broadband)
        else:
            ambience = ambience_from_raster(self.ambience_src, mic)

        # Run the point mesh step a maximum of two times.
        for j in range(2):
            source_pts = build_src_point_mesh(active_space, src_pt_density, altitude_m)
            new_audibility_pts = self._run_nmsim(
                f"{mic.name}_mesh{j + 1}",
                source_pts,
                crs,
                flt_filename,
                site_filename,
                omni_source,
                ambience,
                heading
            )

            tested_space = pd.concat([tested_space, new_audibility_pts], ignore_index=True) 
            active_space = tested_space[tested_space.audible == 1]

            # Create a small buffer around the extent of audible points. If the padded active space is less than
            #  30% smaller than the original study area before the for loop is completed, break out of the for
            #  loop and proceed to the contouring step.
            minx, miny, maxx, maxy = active_space.total_bounds
            xpad = 0.2 * (maxx - minx)  # pad extents by 20% on each side, 40% total.
            ypad = 0.2 * (maxy - miny)
            extent = ([minx - xpad, maxx + xpad], [miny - ypad, maxy + ypad])
            shrinkage = np.divide(np.diff(extent) - np.diff(study_area_extent), np.diff(study_area_extent))
            if min(shrinkage) > -0.30:
                break

        # Run triangulation n_counter times to refine the edges of the active space.
        for k in range(n_contour):
            source_pts = self._contour_active_space(tested_space, altitude_m)
            new_audibility_pts = self._run_nmsim(
                f"{mic.name}_contour{k + 1}",
                source_pts,
                crs,
                flt_filename,
                site_filename,
                omni_source,
                ambience,
                heading
            )
            tested_space = pd.concat([tested_space, new_audibility_pts], ignore_index=True)

        active_space = self._build_active_space(tested_space, crs)
        active_space['altitude_m'] = altitude_m
        active_space['mic_name'] = mic.name

        return active_space

    def set_dem(self, mic: Microphone):
        """
        Projecting and masking a DEM file are bottleneck steps in the active space creation process. If active
        space generation is going to be run for the same location just with different parameters like omni source,
        altitude, etc. there is no reason to project and mask the DEM every time.
        This function provides a way to only project and mask the DEM for the study area once. Then, every time the
        generate() function is run, it will use the created DEM file.

        NOTE: This function is only useful when running generate(). Running generate_mesh() will overwrite anything
        set by this function because it's not applicable.

        Parameters
        ----------
        mic : Microphone
            Microphone object acting as the "listener" in NMSIM.
        """
        crs = NMSIM_bbox_utm(self.study_area.iloc[[0]])
        projected_mic = mic.to_crs(crs)

        self._dem_file = self._mask_dem_file(self.dem_src, study_area=self.study_area.iloc[[0]], project=True)
        self._flt_file = self._create_dem_flt(self._dem_file)
        self._site_file = self._create_site_file(projected_mic, self._flt_file)

    def generate(self, omni_source: str, altitude_m: int = 3658, mic: Optional[Microphone] = None,
                 heading: Optional[int] = None, src_pt_density: int = 48, n_contour: int = 1) -> gpd.GeoDataFrame:
        """
        Generate an active space for the study area.

        Parameters
        ----------
        omni_source : str
            Absolute path to the omni source tuning file to use when running NMSIM.
        altitude_m : int, default 3658 meters (equivalent to 12000 ft)
            Single altitude value to use when creating NSMIM trajectories.
        mic : Microphone, default None
            A Microphone object to use as the NMSIM receiver. If no Microphone is passed, the study area centroid
            will be used as the NMSIM receiver location.
        heading : int, default None
            The heading (yaw) to use for all points in the trajectory file. If None, a random heading will be used
            for each point.
        src_pt_density : int
            Density of the point mesh to be used in the first two rounds of active space definition. The point mesh will
            have src_pt_density x src_point_density points.
        n_contour : int, default 1
            Number of rounds of contouring to perform after the two rounds of active space point meshing.

        Returns
        -------
        active_space : gpd.GeoDataFrame
            A GeoDataFrame of the generated active space polygon.
        """
        active_space = self._generate(
            study_area=self.study_area.iloc[[0]],   # Select the study area so that it's a 1 row GeoDataFrame.
            dem_file=self.dem_src,
            omni_source=omni_source,
            mic=mic,
            project_dem=True,
            altitude_m=altitude_m,
            heading=heading,
            src_pt_density=src_pt_density,
            n_contour=n_contour,
        )
        return active_space

    def generate_mesh(self, omni_source: str, altitude_m: int = 3658, heading: Optional[int] = None,
                      src_pt_density: int = 48, n_contour: int = 1, mesh_density: Tuple[int, int] = (1, 25),
                      n_cpus: int = mp.cpu_count() - 1) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """
        Generate multiple active spaces in a mesh pattern for the study area.

        Parameters
        ----------
        omni_source : str
            Absolute path to the omni source tuning file to use when running NMSIM.
        altitude_m : int, default 3658 meters (equivalent to 12000 ft)
            Single altitude value to use when creating NSMIM trajectories.
        heading : int, default None
            The heading (yaw) to use for all points in the trajectory file. If None, a random heading will be used
            for each point.
        src_pt_density : int
            Density of the point mesh to be used in the first two rounds of active space definition. The point mesh will
            have src_pt_density x src_point_density points.
        n_contour : int, default 1
            Number of rounds of contouring to perform after the two rounds of active space point meshing.
        mesh_density : Tuple[int, int], default (1km, 25km)
            Coarseness of the mesh in kilometers. The first value is how far apart the mesh centroids should be and
            the second value is how large the mesh squares around the centroids should be, both in kilometers.
        n_cpus : int, default mp.cpu_count() - 1
            How many cpus to use for multiprocessing the mesh. Defaults to the total number of computer cpus.

        Returns
        -------
        A GeoDataFrame of all generated active space polygons.
        A GeoDataFrame of centroids used to make the mesh.
        """
        self._dem_file = None
        self._flt_file = None
        self._site_file = None

        study_areas, centroids = create_overlapping_mesh(self.study_area, mesh_density[0], mesh_density[1])
        centroids['name'] = centroids.apply(lambda x: f"centroid{x.name+1}", axis=1)

        dem_file = self._mask_dem_file(self.dem_src, self.study_area, project=True, buffer=mesh_density[1] + 1)

        # Since most arguments are the same for each process, create a partial.
        _generate = partial(self._generate, dem_file=dem_file, omni_source=omni_source, project_dem=False,
                            altitude_m=altitude_m, heading=heading, src_pt_density=src_pt_density, n_contour=n_contour)

        pbar = tqdm(desc='Study Area', unit='study area', colour='green', total=study_areas.shape[0], leave=True)
        _update_pbar = lambda _: pbar.update()
        _handle_error = lambda error: print(f'Error: {error}', flush=True)

        with mp.Pool(n_cpus) as pool:
            processes = []
            for i in range(study_areas.shape[0]):
                processes.append(pool.apply_async(_generate,
                                                  kwds={'study_area': study_areas.iloc[[i]], 'name': f'{i+1}'},
                                                  callback=_update_pbar,
                                                  error_callback=_handle_error))
            results = [p.get() for p in processes]
            active_spaces = results.pop()
            for res in results:
                active_spaces = pd.concat([active_spaces, res], ignore_index=True)

        return active_spaces, centroids.to_crs(active_spaces.crs)
