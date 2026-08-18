import glob
import logging
import multiprocessing as mp
import os
import subprocess
from functools import partial
from pathlib import Path
from typing import Optional, Tuple, Union
from uuid import uuid4

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from osgeo import gdal
from pyproj import Transformer
import rasterio
from shapely.geometry import Polygon, box
from shapely.validation import make_valid
from tqdm import tqdm
from warnings import warn

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.setup.site_writer import create_site_dir, write_listener_site_file
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.utils.models import Microphone
from nps_active_space.utils.computation import (
    build_src_point_mesh,
    coords_to_utm,
    create_overlapping_mesh,
    NMSIM_bbox_utm,
    project_raster,
    round_points
)

# dB threshold of human hearing at each 1/3 octave band, from ISO 226:2003
human_hearing_threshold = pd.Series({
    "20": 78.5, "25": 68.7, "31.5": 59.5, "40": 51.1, "50": 44.0, "63": 37.5, "80": 31.5,
    "100": 26.5, "125": 22.1, "160": 17.9, "200": 14.4, "250": 11.4, "315": 8.6, "400": 6.2,
    "500": 4.4, "630": 3.0, "800": 2.2, "1000": 2.4, "1250": 3.5, "1600": 1.7, "2000": -1.3,
    "2500": -4.2, "3150": -6.0, "4000": -5.4, "5000": -1.5, "6300": 6.0, "8000": 12.6,
    "10000": 13.9, "12500": 12.3
})

logger = logging.getLogger(__name__)

class PolygonCreationError(Exception):
    pass

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
        This directory is specific to a single microphone location.
    dem_src : str
        Path to a DEM raster file to be used as NMSIM input.
    ambience : float or pd.Series[float]
        The ambience level(s) at the microphone site.
        If float (broadband ambience), will be compared against the predicted A-weighted broadband level of noises.
        If pd.Series[float], should contain sound levels for the 12.5 to 12500 Hz 1/3 octave bands.
    """
    def __init__(self, NMSIM: str, study_area: gpd.GeoDataFrame, root_dir: str, dem_src: str,
                 ambience: Union[float, pd.Series]):

        assert os.path.exists(NMSIM), "NMSIM not found"
        assert os.path.exists(dem_src), "DEM not found"
        assert os.path.exists(root_dir), "Root directory not found"
        assert isinstance(ambience, (float, int)) or isinstance(ambience, pd.Series), "Improper ambience input"
        if isinstance(ambience, (float, int)):
            warn("Using broadband ambience. This feature has not been maintained and has possible buggy, incorrect, "
                 "or unexpected behavior. Only use if you know what you are doing.", UserWarning)

        self.study_area = study_area.to_crs('epsg:4269')
        self.root_dir = root_dir
        self.ambience = ambience
        self.dem_src = dem_src
        self.NMSIM = NMSIM

        self._dem_file = None
        self._flt_file = None
        self._site_file = None

        create_site_dir(self.root_dir)

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
            equal_area_crs,_ = coords_to_utm(study_area.centroid.iat[0].y, study_area.centroid.iat[0].x)
            study_area_m = study_area.to_crs(equal_area_crs)
            study_area_m = study_area_m.buffer(buffer*1000)
            study_area = study_area_m.to_crs(study_area.crs)

        # sometimes ArcGIS Pro will sneak in an immutable 'FID' column
        # thankfully if it exists, a viable solution is just to drop the column
        if 'FID' in study_area.columns:
            study_area = study_area.drop(columns=['FID'])
        
        # frequently the columns in the `study_area` `gpd.GeoDataFrame` are also immutable
        # to fix this we can convert them to the `object` type which allows us to write the study area
        # to the disk temporarily as we need to here...
        for col in study_area.columns:
            study_area[col] = study_area[col].astype(object)

        # we avoid the deprecated `pd.Int64Index` and an associated AttributeError
        # by simply setting the `index` parameter to False...
        study_area.to_file(study_area_filename, index=False) 

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
        old_hdr = pd.read_csv(flt_header_filename, header=None, sep=r'\s+', index_col=0).T

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

    def _create_trajectory_file(self, points: gpd.GeoDataFrame, crs: str, filename: str,
                                heading: Optional[int] = None) -> str:
        """
        Create a trajectory file from a list of points.

        Parameters
        ----------
        points : gpd.GeoDataFrame
            GeoDataFrame of trajectory points to be written and run through NMSIM.
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

        trajectory = points.to_crs(crs)
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
        return str(write_listener_site_file(
            self.root_dir, mic.name, mic.x, mic.y, mic.z, dem_file
        ))

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

    def _postprocess_trj_tis(self, trajectory_file, tis_file, cleanup: bool = True):
        """
        Reads a job's input trajectory file and output TIS file, and combines them into a DataFrame.
        Optionally removes the .trj and .tis files after doing so.

        Parameters
        ----------
        trajectory_file: str
            Absolute path to the job's trajectory file.
        tis_file: str
            Absolute path to the job's output TIS file.
        cleanup: bool, default True
            If True, delete the .trj and .tis files after extracting their data.

        Returns
        -------
        nmsim_df: pd.DataFrame
            DataFrame representing tested point audibility. Has the columns:
            Xpos, Ypos, Zpos, 
        """
        # Read in the trajectory input file.
        traj_df = pd.read_fwf(trajectory_file, header=14, widths=[16, 14] + [15]*7)
        # drop unneeded fields
        traj_df = traj_df.drop(["time(s)", "heading", "climbANG", "Vel", "power", "rol"], axis=1)

        # Read in the tis NMSIM output file.
        tis_df = pd.read_fwf(tis_file, header=15, skipfooter=1, widths=[4, 12] + [5]*35)
        tis_df.drop(['SP#', 'F', '42'], axis=1, inplace=True)  # Drop the sensitivity index column
        tis_df.drop(tis_df.head(2).index, inplace=True)        # Drop dividing *****, and unneeded AMBIENT row
        tis_df.reset_index(drop=True, inplace=True)
        tis_df.columns = ["TIME", "A", "10", "12.5", "15.8", "20", "25", "31.5", "40", "50", "63",
                            "80", "100", "125", "160", "200", "250", "315", "400", "500", "630", "800", "1000",
                            "1250", "1600", "2000", "2500", "3150", "4000", "5000", "6300", "8000", "10000", "12500"]
        tis_df = tis_df.drop("TIME", axis=1)
        # convert centibels (cB) to decibels (dB) and round to remove floating point precision errors
        # file is less precise than 6 decimal places, so rounding just removes floating point issues
        # all columns are sound levels now so can apply this to the whole df
        tis_df = (tis_df.astype(float) * 0.1).round(6)

        assert not tis_df.empty, "NMSIM didn't run, try reducing max_pts in _preprocess_source_pts() and try again"
        assert len(traj_df) == len(tis_df), f"# trajectory points ({len(traj_df)}) is not equal to # tis points ({len(tis_df)})"
        new_rows = pd.concat([traj_df, tis_df], axis=1)
        
        if cleanup:
            os.remove(trajectory_file)
            os.remove(tis_file)
        
        return new_rows

    def _find_audible_points(self, nmsim_df: pd.DataFrame, crs: str) -> gpd.GeoDataFrame:
        """
        Determine which points from a trajectory file are audible given the corresponding NMSIM tis output.

        Parameters
        ----------
        nmsim_df: pd.DataFrame
            DataFrame listing source point coordinates and NMSIM's predictions in dB of the sound level
            at the microphone for each band. Contains at least columns "Xpos", "Ypos", "Zpos", "A",
            and 1/3 octave bands "20" through "12500", in order.
        crs : str
            crs of the trajectory file and of the output GeoDataFrame. In the format 'epsg:XXXX'

        Returns
        -------
        total_space ; gpd.GeoDataFrame
            A GeoDataFrame of the NMSIM tested points from the trajectory file with an 'audible' column set to 0 or 1
            depending on the point's audibility.
        """
        # Check to see if any of the frequency bands are louder than the ambient levels.
        if isinstance(self.ambience, (float, int)):
            # broadband ambience
            audible = nmsim_df.loc[:, "A"] > self.ambience
        else:
            # spectral ambience
            # the audibility threshold is the maximum of ambience and the limit of human hearing
            # we only have human hearing threshold in bands 20-12500, so just use those
            threshold = np.maximum(self.ambience.loc["20":"12500"], human_hearing_threshold)
            audible = (nmsim_df.loc[:, "20":"12500"] > threshold.values).any(axis=1)

        return gpd.GeoDataFrame(
            {"audible": audible.astype(int)},
            crs=crs,
            geometry=gpd.points_from_xy(nmsim_df["Xpos"], nmsim_df["Ypos"], nmsim_df["Zpos"])
        )

    @staticmethod
    def _contour_active_space(total_space: gpd.GeoDataFrame, altitude: int) -> gpd.GeoDataFrame:
        """
        Use triangulation to select points along the audible/inaudible line of the active space to more precisely
        define the boundaries.

        Parameters
        ----------
        total_space : gpd.GeoDataFrame
            All points that have been tested for audibility so far -- both audible and inaudible.
        altitude : int
            The altitude (in meters) we are calculating the active space at.

        Returns
        -------
        points: gpd.GeoDataFrame
            GeoDataFrame of points to pass into NMSIM along the audible/inaudible line.
        """
        # Uses Delaunay triangulation
        fig, ax = plt.subplots()
        tri = mpl.tri.Triangulation(total_space.geometry.x.tolist(), total_space.geometry.y.tolist())
        cs = ax.tricontour(tri, total_space.audible.tolist(), levels=[0.5])  # contour with arbitrary point cloud
        plt.close(fig)

        contour_path = cs.get_paths()[0]
        x = contour_path.vertices[:, 0]
        y = contour_path.vertices[:, 1]
        pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x, y, altitude), crs=total_space.crs)

        return pts

    def _determine_underground_pts(self, source_pts: gpd.GeoDataFrame) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """Determines which source points are aboveground / underground.
        
        Parameters
        ----------
        source_pts: gpd.GeoDataFrame
            GeoDataFrame of 3D points.
        
        Returns
        -------
        aboveground_pts: gpd.GeoDataFrame
            Subset of source_pts that are above ground.
        underground_pts: gpd.GeoDataFrame
            Subset of source_pts that are below ground.
        """
        crs = source_pts.crs.to_string().lower()

        aboveground_indices = []
        underground_indices = []  # underground or no DEM data
        with rasterio.open(self._dem_file) as dem:
            proj = Transformer.from_crs(crs, dem.crs, always_xy=True)
            xs = source_pts.geometry.x
            ys = source_pts.geometry.y
            zs = source_pts.geometry.z
            xs, ys = proj.transform(xs, ys)
            proj_pts = np.stack([xs, ys], axis=1)
            elevs = list(dem.sample(proj_pts))
            for i in range(len(source_pts)):
                raw = elevs[i]
                if raw is None:
                    underground_indices.append(i)
                    continue
                elev_m = float(raw[0])
                if elev_m == float(dem.nodata) or zs.iloc[i] < elev_m:
                    underground_indices.append(i)
                else:
                    aboveground_indices.append(i)
        
        return source_pts.iloc[aboveground_indices], source_pts.iloc[underground_indices]

    @staticmethod
    def _nmsim_cache_is_readable(csv_filename: str) -> bool:
        """Return True when ``csv_filename`` exists and contains parseable NMSIM cache rows."""
        if not os.path.exists(csv_filename) or os.path.getsize(csv_filename) == 0:
            return False
        try:
            preview = pd.read_csv(csv_filename, nrows=1)
        except pd.errors.EmptyDataError:
            return False
        required = {"Xpos", "Ypos", "A"}
        return not preview.empty and required.issubset(preview.columns)

    @staticmethod
    def load_prev_nmsim_predictions(source_pts: gpd.GeoDataFrame, csv_filename: str, altitude_m: int
                                     ) -> Tuple[pd.DataFrame, pd.DataFrame, gpd.GeoDataFrame]:
        """
        Loads previous NMSIM predictions and compares them against source points we want to compute
        to see if any have been previously computed. This method is static so external scripts can use
        it to examine NMSIM predictions easily.

        Parameters
        ----------
        source_pts: gpd.GeoDataFrame
            GeoDataFrame of 3D source points we wish to get NMSIM predictions for.
        csv_filename: str
            Path to CSV file containing previous NMSIM predictions.
        altitude_m: int
            Altitude of the points in the CSV file.

        Returns
        -------
        nmsim_df_all: pd.DataFrame
            DataFrame containing all past NMSIM predictions.
        nmsim_df: pd.DataFrame
            DataFrame containing past NMSIM predictions that correspond to some source point.
        new_pts: gpd.GeoDataFrame
            Subset of source_pts containing only the points we don't have previous results for.
        """
        if not ActiveSpaceGenerator._nmsim_cache_is_readable(csv_filename):
            if os.path.exists(csv_filename):
                logger.warning(
                    "Ignoring unreadable NMSIM prediction cache (empty or corrupt): %s",
                    csv_filename,
                )
            # no previous predictions, return empty dataframes
            nmsim_df_all = pd.DataFrame()
            nmsim_df = pd.DataFrame()
            new_pts = source_pts
        else:
            # read in previous CSV
            # note: to save disk space / file read-write time, we used NA to represent no sound aka -99.9dB,
            # and stored centibels not decibels to avoid writing the decimal point. So, convert back.
            # Also add back in Zpos field, which we didn't store since it is constant within the file.
            nmsim_df_all = pd.read_csv(csv_filename).fillna(-999).astype("float64")
            if not nmsim_df_all.empty:
                sound_cols = [
                    col for col in nmsim_df_all.columns
                    if col not in {"Xpos", "Ypos", "Zpos"}
                ]
                nmsim_df_all[sound_cols] /= 10
            nmsim_df_all["Zpos"] = altitude_m
            
            # build two multi-indices to quickly determine which points we've run before
            source_idx = pd.MultiIndex.from_frame(pd.DataFrame({
                "Xpos": source_pts.geometry.x,
                "Ypos": source_pts.geometry.y,
            }))
            prev_idx = pd.MultiIndex.from_frame(nmsim_df_all[["Xpos", "Ypos"]])

            nmsim_df = nmsim_df_all[prev_idx.isin(source_idx)].drop_duplicates(["Xpos", "Ypos"])
            new_pts = source_pts[~source_idx.isin(prev_idx)].drop_duplicates("geometry")
            # each source pt should be represented by a row in either new_pts or nmsim_df
            assert len(new_pts) + len(nmsim_df) == len(source_pts)
        
        return nmsim_df_all, nmsim_df, new_pts

    @staticmethod
    def save_nmsim_predictions(nmsim_df_all: pd.DataFrame, csv_filename: str):
        """
        Saves NMSIM predictions to a CSV file for future reference.
        Compresses the data somewhat for better read-write performance and disk space usage.

        Parameters
        ----------
        nmsim_df_all: pd.DataFrame
            DataFrame of NMSIM predictions. Should contain columns Xpos, Ypos, Zpos (optional), A, and 1/3 octave bands
        csv_filename: str
            Path to CSV file to store NMSIM predictions in.
        """
        if nmsim_df_all.empty:
            return

        # remove duplicate points, this can happen sometimes I think and cause counting weirdness
        nmsim_df_all = nmsim_df_all.drop_duplicates(subset=["Xpos", "Ypos"])

        # Save to CSV, saving disk space and read-write time by:
        # - using centibels to avoid writing the decimal point
        # - storing no sound (-99.9dB) as NA
        # - omitting the constant-valued Zpos field
        dB_cols = nmsim_df_all.loc[:,"A":"12500"].columns
        nmsim_df_all[dB_cols] = (nmsim_df_all[dB_cols] * 10).astype("Int64").replace(-999, pd.NA)
        nmsim_df_all.drop("Zpos", axis=1, inplace=True, errors="ignore")
        nmsim_df_all.to_csv(csv_filename, index=False)

    def _run_nmsim(self, job_name: str, source_pts: gpd.GeoDataFrame, flt_file: str, site_file: str,
                   omni_source: str, altitude_m: int, heading: Optional[int] = None) -> gpd.GeoDataFrame:
        """
        Execute a single NMSIM job.

        Parameters
        ----------
        job_name : str
            Name of this NMSIM run to use a suffix to input and output files.
        source_pts : gpd.GeoDataFrame
            GeoDataFrame of points to test the audibility of in NMSIM.
        flt_file : str
            Absolute path to the .flt DEM file required to run NMSIM.
        site_file : str
            Absolute path to the .sit file of the receiver point required to run NMSIM.
        omni_source : str
            Absolute path to the omni source file to use.
        altitude_m : int
            The altitude of the source points, in meters.
        heading : int, default None
            The heading (yaw) to use for all points in the trajectory file. If None, a random heading will be used
            for each point.

        Returns
        -------
        new_audibility_pts : gpd.GeoDataFrame
            A GeoDataFrame of points tested during the NMSIM run and their audibility.
        """
        assert len(source_pts) > 0, "Trying to run NMSIM on zero source points"
        source_pts = source_pts.drop_duplicates("geometry")  # duplicates cause problems when storing prev nmsim predictions
        crs = source_pts.crs.to_string().lower()

        # Mark any underground points as inaudible and don't pass them to NMSIM
        aboveground_pts, underground_pts = self._determine_underground_pts(source_pts)
        
        # mark underground points as inaudible
        audibility_pts = underground_pts
        audibility_pts["audible"] = 0
        if len(aboveground_pts) == 0:
            return audibility_pts

        # Check if we've run any of the aboveground points through NMSIM before
        # The csv filename is important - we assume all gains, altitudes, and headings in a csv are the same,
        # and so omit this information inside the csv to save space / read-write time
        omni_str = os.path.splitext(os.path.basename(omni_source))[0]
        csv_filename = f"{self.root_dir}/Output_Data/TIG_TIS/{altitude_m}m_{omni_str}_{heading}deg.csv"
        nmsim_df_all, nmsim_df, new_pts = ActiveSpaceGenerator.load_prev_nmsim_predictions(
            aboveground_pts, csv_filename, altitude_m)
        # print(f"{job_name} n={len(aboveground_pts)}, old={len(nmsim_df)}, new={len(new_pts)}")

        if len(new_pts) == 0:
            # no need to run NMSIM, we have all the predictions already
            new_nmsim_df = pd.DataFrame()
        else:
            # Prepare NMSIM instruction files
            trajectory_filename = self._create_trajectory_file(new_pts, crs, job_name, heading)
            batch_file = self._create_instruction_files(flt_file, site_file, trajectory_filename, omni_source)

            # Run NMSIM.
            process = subprocess.Popen([self.NMSIM, batch_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout, stderr = process.communicate()
            if stderr:
                for s in stderr.decode("utf-8").split("\r\n"):
                    logger.error(s.strip())

            # Read NMSIM outputs and remove unneeded .trj and .tis file
            new_nmsim_df = self._postprocess_trj_tis(
                trajectory_filename,
                f"{self.root_dir}/Output_Data/TIG_TIS/{job_name}.tis",
                cleanup=True)
        
            # Combine with results from previous NMSIM runs
            nmsim_df = pd.concat([nmsim_df, new_nmsim_df], ignore_index=True)
            nmsim_df = nmsim_df.drop_duplicates(subset=["Xpos", "Ypos"])
        
        # Combine new NMSIM predictions with ALL previous predictions (not just ones matching source_pts),
        # and save back to the csv file
        nmsim_df_all = pd.concat([nmsim_df_all, new_nmsim_df], ignore_index=True)
        ActiveSpaceGenerator.save_nmsim_predictions(nmsim_df_all, csv_filename)

        # Determine the audibility of points that were tested.
        nmsim_audibility_pts = self._find_audible_points(nmsim_df, crs)
        audibility_pts = pd.concat([audibility_pts, nmsim_audibility_pts], ignore_index=True)

        return audibility_pts

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
        # augment total space with inaudible points around the boundary, to avoid contour artifacts
        minx, miny, maxx, maxy = total_space.total_bounds
        region = gpd.GeoDataFrame({"geometry": [box(minx, miny, maxx, maxy)]}, crs=total_space.crs)
        region.geometry = region.buffer(100)  # small buffer so that the active space boundary occurs just outside the study area
        new_points = build_src_point_mesh(region)
        new_points = new_points[~new_points.geometry.within(region.union_all())]
        new_points["audible"] = 0
        total_space = pd.concat([total_space, new_points], ignore_index=True)

        # create the triangulated irregular network and contour lines
        fig, ax = plt.subplots()  # need an axis to call tricontour function
        levels = np.linspace(0, 1, 10, endpoint=False)
        tri = mpl.tri.Triangulation(total_space.geometry.x.tolist(), total_space.geometry.y.tolist())
        cs = ax.tricontour(tri, total_space.audible.tolist(), levels=levels)

        # pick some contour level to choose as the active space boundary... 0.5 is somewhat arbitrary
        level_ind = np.where(cs.levels == 0.5)[0][0]
        plt.close(fig)  # close triangulation figure

        # iterate through all contour paths in the `TriContourSet` at level_ind
        # https://matplotlib.org/stable/api/tri_api.html#matplotlib.tri.TriContourSet
        contour_path = cs.get_paths()[level_ind] # in recent versions of `matplotlib` there is 1:1 correspondence `cs.levels` : `Path`
        polygons = [Polygon(P) for P in contour_path.to_polygons()] # convert to `shapely.Polygon`

        # mark inner polygons as holes - otherwise when we dissolve or union_all(), they will disappear
        outer_polys = []
        for i, outer in enumerate(polygons):
            # skip if this is contained in another polygon
            if any(outer.within(other) for j, other in enumerate(polygons) if i != j):
                continue
            # find polygons inside this one
            inner_polys = [inner for j, inner in enumerate(polygons) if j != i and inner.within(outer)]
            holes = [inner.exterior.coords for inner in inner_polys]
            # construct polygon with holes
            poly = Polygon(shell=outer.exterior.coords, holes=holes)
            outer_polys.append(poly)
        
        # ensure valid geometries
        outer_polys = [make_valid(poly) if not poly.is_valid else poly for poly in outer_polys]

        return gpd.GeoDataFrame(data={'geometry': outer_polys}, geometry='geometry', crs=crs)

    def _preprocess_source_points(self, source_pts: gpd.GeoDataFrame, valid_query_region: gpd.GeoDataFrame,
                                  tested_pts: gpd.GeoDataFrame, max_pts: int = 4000):
        """        
        Filter a set of source points, removing points that:
        1) Are outside the valid query region - meaning outside the study area, or too close
           to the study area boundary (since DEM info isn't present outside the study area)
        2) We already know are audible / inaudible
        3) More than what NMSIM can handle (~4000 points at once)

        Parameters
        ----------
        source_pts: gpd.GeoDataFrame
            GeoDataFrame of 3D points.
        valid_query_region: gpd.GeoDataFrame
            GeoDataFrame containing geometry that defines the valid query region. Any source points outside
            of this region will be filtered out and not tested. The valid query region is the study area,
            minus a padding region around the edge so that points are tested too close to where DEM data is missing.
        tested_pts: gpd.GeoDataFrame
            A GeoDataFrame of 3D points, with a field "audible" = 0 or 1, representing points that have been
            tested already by NMSIM. There is no need to retest a point that has already been tested.
        max_pts: int, default 4000
            The maximum number of source points to give to NMSIM in one trajectory file. If len(source_pts)
            exceeds this, points will be dropped at random (but with a fixed seed) to match this.
            Reason: If too many points are given to NMSIM, NMSIM won't calculate anything and will leave
            the output TIS file blank. Annoyingly, this behavior is not deterministic / there is not a hard threshold;
            it becomes more likely with more points (potentially a memory allocation thing that fails silently?).
        
        Returns
        -------
        source_pts: gpd.GeoDataFrame
            A copy of source_pts filtered appropriately.
        """
        # only query points inside the study area and far enough from the study area boundary
        source_pts = source_pts[source_pts.within(valid_query_region)]

        # don't query points we already know the answer for
        source_pts = source_pts[~source_pts.geometry.isin(tested_pts.geometry)]

        # if too many points for NMSIM, randomly downsample
        if source_pts.shape[0] > max_pts:
            source_pts = source_pts.sample(max_pts, random_state=5)       

        return source_pts
    
    def _generate(self, study_area: gpd.GeoDataFrame, dem_file: str, omni_source: str, name: str = '',
                  mic: Optional[Microphone] = None, project_dem: bool = True, altitude_m: int = 3658,
                  heading: Optional[int] = None, src_pt_density: int = 48, n_contour: int = 1,
                  predetermined_audibility_pts: Optional[gpd.GeoDataFrame] = None
                  ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """
        The main active space generating function. It has been separated from the other generate functions to allow
        for multiprocessing when creating an active space mesh.
        """
        crs = NMSIM_bbox_utm(study_area)  # Determine the UTM CRS on the western-most edge of the study area

        # Initialize a GeoDataFrame of source points that have gone through NMSIM
        if predetermined_audibility_pts is not None:
            tested_pts = predetermined_audibility_pts
        else:
            # start with empty tested space
            tested_pts = gpd.GeoDataFrame(columns=['audible', 'geometry'], geometry='geometry', crs=crs)

        # Initialize a GeoDataFrame of the current active space. The active space will initially
        # be the same as the study area, but will be refined.
        study_area = study_area.to_crs(crs)
        valid_query_region = study_area.to_crs(crs).union_all().buffer(-100)  # require points to not be right on the boundary
        study_area_extent = ([study_area.total_bounds[0], study_area.total_bounds[2]],  # ([minx, maxx],
                             [study_area.total_bounds[1], study_area.total_bounds[3]])  # [miny, maxy])

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
        
        # Use meter elevation artifacts from project_setup (required before generate_active_space).
        if self._dem_file is None or self._flt_file is None:
            tif_path, flt_path = get_project_setup_elevation(self.root_dir)
            dem_filename = str(tif_path)
            flt_filename = str(flt_path)
        else:
            dem_filename = self._dem_file
            flt_filename = self._flt_file
        site_filename = self._site_file or self._create_site_file(mic, flt_filename)

        # Prepare a coarse and fine grid to use for the 1st and 2nd point mesh steps
        # we end up rounding the source_pts coords to the nearest 0.001m later, so do this now
        # to make comparisons with the output of past runs work properly
        coarse_grid = build_src_point_mesh(study_area, src_pt_density, altitude_m)
        fine_grid = build_src_point_mesh(study_area, 2*src_pt_density-1, altitude_m)
        round_points(coarse_grid, 3)
        round_points(fine_grid, 3)

        # Run the point mesh step two times.
        for j in range(2):
            if j == 0:
                source_pts = coarse_grid
            elif j == 1:
                # use the fine grid, but only near the boundary to be efficient.
                # if you are within a short distance of an audible and an inaudible point,
                # you are near the boundary - can use .buffer() to figure this out
                audible_pts = tested_pts[tested_pts["audible"] == 1]
                inaudible_pts = tested_pts[tested_pts["audible"] != 1]
                # Buffer - reduce the buffering memory load by only buffering from the coarse grid.
                # Not reducing memory load can cause out-of-memory crashes
                audible_pts = gpd.sjoin(audible_pts, coarse_grid, how="inner", predicate='intersects')
                inaudible_pts = gpd.sjoin(inaudible_pts, coarse_grid, how="inner", predicate='intersects')
                # buffer by about 1.5 coarse grid cells, emprically is enough
                buffer_amt = 1.5 * np.diff(np.array(study_area_extent)).max() / src_pt_density
                # Also use cap_style=3 for square buffers to further reduce memory load.
                near_audible = audible_pts.union_all().buffer(buffer_amt, cap_style=3)
                near_inaudible = inaudible_pts.union_all().buffer(buffer_amt, cap_style=3)

                boundary_zone = near_audible.intersection(near_inaudible)
                source_pts = fine_grid[fine_grid.within(boundary_zone)]
            
            source_pts = self._preprocess_source_points(source_pts, valid_query_region, tested_pts)
            if source_pts.empty:
                logger.info(f"Mesh step {j+1}: no source points, skipping")
                break
            
            new_audibility_pts = self._run_nmsim(
                f"{mic.name}_{altitude_m}m_mesh{j + 1}",
                source_pts,
                flt_filename,
                site_filename,
                omni_source,
                altitude_m,
                heading
            )
            tested_pts = pd.concat([tested_pts, new_audibility_pts], ignore_index=True) 

        # Run triangulation n_contour times to refine the edges of the active space.
        for k in range(n_contour):
            source_pts = self._contour_active_space(tested_pts, altitude_m)
            # we end up rounding the source_pts coords to the nearest 0.001m later, so do this now
            # to make comparisons with the output of past runs work properly
            round_points(source_pts, 3)
            source_pts = self._preprocess_source_points(source_pts, valid_query_region, tested_pts)
            if source_pts.empty:
                # print(f"Refine step {k+1}: no source points, skipping")
                break
            new_audibility_pts = self._run_nmsim(
                f"{mic.name}_{altitude_m}m_contour{k + 1}",
                source_pts,
                flt_filename,
                site_filename,
                omni_source,
                altitude_m,
                heading
            )
            tested_pts = pd.concat([tested_pts, new_audibility_pts], ignore_index=True)

        active_space = self._build_active_space(tested_pts, crs)
        active_space['altitude_m'] = altitude_m
        active_space['mic_name'] = mic.name

        return active_space, tested_pts

    def set_dem(self, mic: Microphone):
        """
        Cache project_setup elevation artifacts and listener ``.sit`` for repeated ``generate()`` calls.

        Requires ``elevation_m_nad83_utm*.tif`` / ``.flt`` / ``.hdr`` from ``project_setup`` under
        ``Input_Data/01_ELEVATION``. Run ``project_setup`` for the site before ``generate_active_space``.

        NOTE: This function is only useful when running generate(). Running generate_mesh() uses a
        separate buffered DEM path.

        Parameters
        ----------
        mic : Microphone
            Microphone object acting as the "listener" in NMSIM.
        """
        crs = NMSIM_bbox_utm(self.study_area.iloc[[0]])
        projected_mic = mic.to_crs(crs)

        tif_path, flt_path = get_project_setup_elevation(self.root_dir)
        self._dem_file = str(tif_path)
        self._flt_file = str(flt_path)
        flt_for_sit = Path(flt_path)
        if flt_for_sit.is_relative_to(Path(self.root_dir)):
            flt_for_sit = flt_for_sit.relative_to(Path(self.root_dir))
        self._site_file = str(write_listener_site_file(
            self.root_dir,
            projected_mic.name,
            projected_mic.x,
            projected_mic.y,
            projected_mic.z,
            flt_for_sit,
        ))
        logger.info(
            "Using project_setup elevation artifacts for NMSIM: %s",
            tif_path.name,
        )

    def generate(self, omni_source: str, altitude_m: int = 3658, mic: Optional[Microphone] = None,
                 heading: Optional[int] = None, src_pt_density: int = 48, n_contour: int = 1,
                 predetermined_audibility_pts: Optional[gpd.GeoDataFrame] = None
                 )-> Union[gpd.GeoDataFrame, Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]]:
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
        predetermined_audibility_pts: gpd.GeoDataFrame, default None
            A GeoDataFrame of points we already know are audible/inaudible for this omni source and heading.
            It's geometry is 3D points, and contains an "audible" field = 0 or 1.
            Use case - if we previously tested a quieter omni source, anywhere it was audible will also be audible
            for a louder omni source. A similar thing holds for louder sources being inaudible at certain points.

        Returns
        -------
        active_space : gpd.GeoDataFrame
            A GeoDataFrame of the generated active space polygon.
        tested_points : gpd.GeoDataFrame
            A GeoDataFrame of 3D points that were tested for audibility, with an "audible" field (0 or 1)
            listing whether they were determined to be audible or not.
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
            predetermined_audibility_pts=predetermined_audibility_pts
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
        _handle_error = lambda error: logger.error(f'Error: {error}')

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
