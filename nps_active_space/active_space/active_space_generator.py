import logging
import os
from pathlib import Path
from typing import Optional, Tuple, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pyproj import Transformer
import rasterio
from shapely.geometry import Polygon, box
from shapely.validation import make_valid
from tqdm import tqdm
from warnings import warn

from nps_active_space.active_space.propagation_model import PropagationModel, prediction_cache_csv_path
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.setup.site_writer import create_site_dir
from nps_active_space.utils.models import Microphone
from nps_active_space.utils.computation import (
    build_src_point_mesh,
    NMSIM_bbox_utm,
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
    A class that stores active space generation logic and produces individual active spaces.

    Parameters
    ----------
    NMSIM : str
        Absolute path to the NMSIM executable to be used to generate active spaces.
    study_area : gpd.GeoDataFrame
        A gpd.GeoDataFrame of polygon(s) that make up the study area.
    root_dir : str
        Absolute path to a directory where all generated files required for running NMSIM can be stored.
        This directory is specific to a single microphone location. Elevation comes from
        ``project_setup`` artifacts under ``Input_Data/01_ELEVATION`` (see ``set_dem``).
    ambience : float or pd.Series[float]
        The ambience level(s) at the microphone site.
        If float (broadband ambience), will be compared against the predicted A-weighted broadband level of noises.
        If pd.Series[float], should contain sound levels for the 12.5 to 12500 Hz 1/3 octave bands.
    """
    def __init__(
        self,
        NMSIM: str | None,
        study_area: gpd.GeoDataFrame,
        root_dir: str,
        ambience: Union[float, pd.Series],
        propagation_model: PropagationModel | None = None,
    ):
        assert os.path.exists(root_dir), "Root directory not found"
        assert isinstance(ambience, (float, int)) or isinstance(ambience, pd.Series), "Improper ambience input"
        if isinstance(ambience, (float, int)):
            warn("Using broadband ambience. This feature has not been maintained and has possible buggy, incorrect, "
                 "or unexpected behavior. Only use if you know what you are doing.", UserWarning)

        if propagation_model is None:
            assert NMSIM is not None and os.path.exists(NMSIM), "NMSIM not found"
            from nps_active_space.active_space.nmsim_propagation_model import NmsimPropagationModel
            propagation_model = NmsimPropagationModel(NMSIM, root_dir)

        self.study_area = study_area.to_crs('epsg:4269')
        self.root_dir = root_dir
        self.ambience = ambience
        self.NMSIM = NMSIM
        self.propagation_model = propagation_model
        self._site_context = None
        self._dem_file = None
        self._flt_file = None

        create_site_dir(self.root_dir)

    def _resolve_site_elevation(self) -> str:
        """Return the project_setup GeoTIFF path, caching tif/flt paths when needed."""
        if self._dem_file is not None:
            return self._dem_file
        tif_path, flt_path = get_project_setup_elevation(self.root_dir)
        self._dem_file = str(tif_path)
        self._flt_file = str(flt_path)
        return self._dem_file

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
        assert self._site_context is not None, "prepare_site must run before propagation predict"
        dem_path = self._site_context.dem_file

        aboveground_indices = []
        underground_indices = []  # underground or no DEM data
        with rasterio.open(dem_path) as dem:
            proj = Transformer.from_crs(crs, dem.crs, always_xy=True)
            xs = source_pts.geometry.x
            ys = source_pts.geometry.y
            zs = source_pts.geometry.z
            xs, ys = proj.transform(xs, ys)
            proj_pts = np.stack([xs, ys], axis=1)
            elevs = list(dem.sample(proj_pts))
            for i in range(len(source_pts)):
                if elevs[i] is None or elevs[i] == dem.nodata or zs.iloc[i] < elevs[i]:
                    underground_indices.append(i)
                else:
                    aboveground_indices.append(i)
        
        return source_pts.iloc[aboveground_indices], source_pts.iloc[underground_indices]

    @staticmethod
    def _nmsim_cache_failure_reason(csv_filename: str) -> str | None:
        """
        Return why an on-disk cache cannot be used, or None if the file is readable.

        Returns None when ``csv_filename`` does not exist (no cache yet) or when the
        file contains valid NMSIM prediction rows.
        """
        if not os.path.exists(csv_filename):
            return None
        if os.path.getsize(csv_filename) == 0:
            return "file is empty (0 bytes)"
        try:
            preview = pd.read_csv(csv_filename, nrows=1)
        except pd.errors.EmptyDataError:
            return "file has no parseable CSV content"
        if preview.empty:
            return "file has a header but no data rows"
        required = {"Xpos", "Ypos", "A"}
        missing = sorted(required - set(preview.columns))
        if missing:
            return f"missing required columns: {missing}"
        return None

    @staticmethod
    def _nmsim_cache_is_readable(csv_filename: str) -> bool:
        """Return True when ``csv_filename`` exists and contains parseable NMSIM cache rows."""
        return (
            os.path.exists(csv_filename)
            and ActiveSpaceGenerator._nmsim_cache_failure_reason(csv_filename) is None
        )

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
            cache_failure = ActiveSpaceGenerator._nmsim_cache_failure_reason(csv_filename)
            if cache_failure is not None:
                logger.warning(
                    "Ignoring unreadable NMSIM prediction cache at %s (%s). "
                    "Removing file; points will be recomputed and cache rewritten after NMSIM.",
                    csv_filename,
                    cache_failure,
                )
                try:
                    os.remove(csv_filename)
                except OSError as exc:
                    logger.warning(
                        "Could not remove unreadable NMSIM prediction cache %s: %s",
                        csv_filename,
                        exc,
                    )
            nmsim_df_all = pd.DataFrame()
            nmsim_df = pd.DataFrame()
            new_pts = source_pts
        else:
            # Stored as centibels with NA for no sound (-99.9 dB); convert back to dB.
            # Zpos is omitted in the file because it is constant within a cache.
            nmsim_df_all = pd.read_csv(csv_filename).fillna(-999).astype("float64")
            if not nmsim_df_all.empty:
                sound_cols = [
                    col for col in nmsim_df_all.columns
                    if col not in {"Xpos", "Ypos", "Zpos"}
                ]
                nmsim_df_all[sound_cols] /= 10
            nmsim_df_all["Zpos"] = altitude_m

            source_idx = pd.MultiIndex.from_frame(pd.DataFrame({
                "Xpos": source_pts.geometry.x,
                "Ypos": source_pts.geometry.y,
            }))
            prev_idx = pd.MultiIndex.from_frame(nmsim_df_all[["Xpos", "Ypos"]])

            nmsim_df = nmsim_df_all[prev_idx.isin(source_idx)].drop_duplicates(["Xpos", "Ypos"])
            new_pts = source_pts[~source_idx.isin(prev_idx)].drop_duplicates("geometry")
            assert len(new_pts) + len(nmsim_df) == len(source_pts)

        return nmsim_df_all, nmsim_df, new_pts

    @staticmethod
    def _source_pts_missing_predictions(
        source_pts: gpd.GeoDataFrame,
        pred_df: pd.DataFrame,
    ) -> gpd.GeoDataFrame:
        """Points sent to predict() that have no row in the returned prediction DataFrame."""
        if len(source_pts) == 0:
            return source_pts.iloc[0:0]
        if len(pred_df) == 0:
            return source_pts

        source_idx = pd.MultiIndex.from_frame(pd.DataFrame({
            "Xpos": source_pts.geometry.x,
            "Ypos": source_pts.geometry.y,
        }))
        pred_idx = pd.MultiIndex.from_frame(pred_df[["Xpos", "Ypos"]].drop_duplicates())
        return source_pts[~source_idx.isin(pred_idx)].drop_duplicates("geometry")

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
        nmsim_df_all[dB_cols] = (
            (nmsim_df_all[dB_cols] * 10).round().astype("Int64").replace(-999, pd.NA)
        )
        nmsim_df_all.drop("Zpos", axis=1, inplace=True, errors="ignore")
        nmsim_df_all.to_csv(csv_filename, index=False)

    def _run_propagation_model(self, job_name: str, source_pts: gpd.GeoDataFrame,
                               omni_source: str, altitude_m: int,
                               heading: Optional[int] = None) -> gpd.GeoDataFrame:
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

        # AAM uses its own ELV grid; GDAL DEM samples can disagree and break whole batches.
        filter_below = getattr(self.propagation_model, "filter_below_terrain", None)
        if filter_below is not None:
            aboveground_pts, below_aam_pts = filter_below(
                self._site_context,
                aboveground_pts,
                job_name=job_name,
            )
            if len(below_aam_pts) > 0:
                below_aam_pts = below_aam_pts.copy()
                below_aam_pts["audible"] = 0
                underground_pts = pd.concat(
                    [underground_pts, below_aam_pts],
                    ignore_index=True,
                )

        # mark underground points as inaudible
        audibility_pts = underground_pts
        if len(audibility_pts) > 0:
            audibility_pts["audible"] = 0
        if len(aboveground_pts) == 0:
            return audibility_pts

        # Check if we've run any of the aboveground points through NMSIM before
        # The csv filename is important - we assume all gains, altitudes, and headings in a csv are the same,
        # and so omit this information inside the csv to save space / read-write time
        omni_str = os.path.splitext(os.path.basename(omni_source))[0]
        csv_filename = prediction_cache_csv_path(
            self.root_dir,
            self.propagation_model.predictions_subdir,
            altitude_m,
            omni_str,
            heading,
        )
        nmsim_df_all, nmsim_df, new_pts = ActiveSpaceGenerator.load_prev_nmsim_predictions(
            aboveground_pts, csv_filename, altitude_m)
        # print(f"{job_name} n={len(aboveground_pts)}, old={len(nmsim_df)}, new={len(new_pts)}")

        if len(new_pts) == 0:
            # no need to run propagation model, we have all the predictions already
            new_nmsim_df = pd.DataFrame()
        else:
            new_nmsim_df = self.propagation_model.predict(
                self._site_context,
                new_pts,
                omni_source,
                altitude_m,
                job_name,
                heading,
            )
            failed_pts = self._source_pts_missing_predictions(new_pts, new_nmsim_df)
            if len(failed_pts) > 0:
                logger.warning(
                    "%s: marking %d point(s) inaudible after predict failure/skip",
                    job_name,
                    len(failed_pts),
                )
                failed_audibility = failed_pts.copy()
                failed_audibility["audible"] = 0
                audibility_pts = pd.concat(
                    [audibility_pts, failed_audibility],
                    ignore_index=True,
                )
            nmsim_df = pd.concat([nmsim_df, new_nmsim_df], ignore_index=True)
            nmsim_df = nmsim_df.drop_duplicates(subset=["Xpos", "Ypos"])

        # Combine new predictions with ALL previous predictions (not just ones matching source_pts),
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
    
    def _generate(self, study_area: gpd.GeoDataFrame, omni_source: str, name: str = '',
                  mic: Optional[Microphone] = None, altitude_m: int = 3658,
                  heading: Optional[int] = None, src_pt_density: int = 48, n_contour: int = 1,
                  predetermined_audibility_pts: Optional[gpd.GeoDataFrame] = None
                  ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """
        The main active space generating function. It has been separated from the other generate functions to allow
        for multiprocessing when generating multiple active spaces in parallel.
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
        
        if self._site_context is None:
            self._site_context = self.propagation_model.prepare_site(
                self._resolve_site_elevation(),
                study_area,
                mic,
                project_dem=False,
                suffix=f"_{mic.name}",
            )

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
            
            source_pts = self._preprocess_source_points(
                source_pts, valid_query_region, tested_pts,
                max_pts=self.propagation_model.max_points_per_run,
            )
            if source_pts.empty:
                logger.info(f"Mesh step {j+1}: no source points, skipping")
                break
            
            new_audibility_pts = self._run_propagation_model(
                f"{mic.name}_{altitude_m}m_mesh{j + 1}",
                source_pts,
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
            source_pts = self._preprocess_source_points(
                source_pts, valid_query_region, tested_pts,
                max_pts=self.propagation_model.max_points_per_run,
            )
            if source_pts.empty:
                # print(f"Refine step {k+1}: no source points, skipping")
                break
            new_audibility_pts = self._run_propagation_model(
                f"{mic.name}_{altitude_m}m_contour{k + 1}",
                source_pts,
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
        Cache project_setup elevation and run model ``prepare_site`` once per microphone.

        Requires ``elevation_m_nad83_utm*.tif`` / ``.flt`` / ``.hdr`` from ``project_setup`` under
        ``Input_Data/01_ELEVATION``. Run ``project_setup`` for the site before active space generation.

        Parameters
        ----------
        mic : Microphone
            Microphone object acting as the "listener" in NMSIM.
        """
        crs = NMSIM_bbox_utm(self.study_area.iloc[[0]])
        projected_mic = mic.to_crs(crs)
        dem_file = self._resolve_site_elevation()
        self._site_context = self.propagation_model.prepare_site(
            dem_file,
            self.study_area.iloc[[0]],
            projected_mic,
            project_dem=False,
            suffix=f"_{projected_mic.name}",
        )
        logger.info("Using project_setup elevation artifacts: %s", Path(dem_file).name)

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
            omni_source=omni_source,
            mic=mic,
            altitude_m=altitude_m,
            heading=heading,
            src_pt_density=src_pt_density,
            n_contour=n_contour,
            predetermined_audibility_pts=predetermined_audibility_pts
        )
        return active_space
