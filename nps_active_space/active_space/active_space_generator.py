import logging
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer
import rasterio
from warnings import warn

from nps_active_space.active_space.active_space_geometry import (
    build_active_space,
    contour_active_space,
)
from nps_active_space.active_space.audibility import audible_points_gdf
from nps_active_space.active_space.prediction_cache import (
    predict_with_cache,
    prediction_cache_csv_path,
)
from nps_active_space.propagation_model.protocol import PropagationModel
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.setup.site_writer import create_site_dir
from nps_active_space.utils.models import Microphone
from nps_active_space.utils.computation import (
    build_src_point_mesh,
    study_area_utm_crs,
    round_points
)

logger = logging.getLogger(__name__)

__all__ = ["ActiveSpaceGenerator"]


class ActiveSpaceGenerator:
    """
    A class that stores active space generation logic and produces individual active spaces.

    Parameters
    ----------
    NMSIM : str, optional
        Absolute path to the NMSIM executable. Required unless ``propagation_model`` is passed.
    study_area : gpd.GeoDataFrame
        A gpd.GeoDataFrame of polygon(s) that make up the study area.
    root_dir : str
        Absolute path to a directory where all generated files required for running the propagation model can be stored.
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
        ambience: float | int | pd.Series,
        propagation_model: PropagationModel | None = None,
    ):
        assert os.path.exists(root_dir), "Root directory not found"
        assert isinstance(ambience, (float, int)) or isinstance(ambience, pd.Series), "Improper ambience input"
        if isinstance(ambience, (float, int)):
            warn("Using broadband ambience. This feature has not been maintained and has possible buggy, incorrect, "
                 "or unexpected behavior. Only use if you know what you are doing.", UserWarning)

        if propagation_model is None:
            assert NMSIM is not None and os.path.exists(NMSIM), "NMSIM not found"
            from nps_active_space.propagation_model.nmsim.model import NmsimPropagationModel
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

    def set_dem(self, mic: Microphone):
        """
        Cache project_setup elevation and run model ``prepare_site`` once per microphone.

        Requires ``elevation_m_nad83_utm*.tif`` / ``.flt`` / ``.hdr`` from ``project_setup`` under
        ``Input_Data/01_ELEVATION``. Run ``project_setup`` for the site before active space generation.

        Parameters
        ----------
        mic : Microphone
            Microphone object acting as the receiver for the propagation model.
        """
        crs = study_area_utm_crs(self.study_area.iloc[[0]])
        projected_mic = mic.to_crs(crs)
        self._prepare_site_context(self.study_area.iloc[[0]], projected_mic)
        dem_file = self._resolve_site_elevation()
        logger.info("Using project_setup elevation artifacts: %s", Path(dem_file).name)

    def generate(
        self,
        omni_source: str,
        altitude_m: int = 3658,
        mic: Microphone | None = None,
        heading: int | None = None,
        src_pt_density: int = 48,
        n_contour: int = 1,
        predetermined_audibility_pts: gpd.GeoDataFrame | None = None,
    ) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """
        Generate an active space for the study area.

        Parameters
        ----------
        omni_source : str
            Absolute path to the omni source tuning file to use when running the propagation model.
        altitude_m : int, default 3658 meters (equivalent to 12000 ft)
            Single altitude value to use when creating source-point trajectories.
        mic : Microphone, default None
            A Microphone object to use as the propagation model receiver. If no Microphone is passed, the study area centroid
            will be used as the receiver location.
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
        return self._generate(
            study_area=self.study_area.iloc[[0]],   # Select the study area so that it's a 1 row GeoDataFrame.
            omni_source=omni_source,
            mic=mic,
            altitude_m=altitude_m,
            heading=heading,
            src_pt_density=src_pt_density,
            n_contour=n_contour,
            predetermined_audibility_pts=predetermined_audibility_pts
        )

    def _generate(
        self,
        study_area: gpd.GeoDataFrame,
        omni_source: str,
        name: str = '',
        mic: Microphone | None = None,
        altitude_m: int = 3658,
        heading: int | None = None,
        src_pt_density: int = 48,
        n_contour: int = 1,
        predetermined_audibility_pts: gpd.GeoDataFrame | None = None,
    ) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """
        The main active space generating function. It has been separated from the other generate functions to allow
        for multiprocessing when generating multiple active spaces in parallel.
        """
        crs = study_area_utm_crs(study_area)  # Determine the UTM CRS on the western-most edge of the study area

        # Initialize a GeoDataFrame of source points that have been tested by the propagation model
        if predetermined_audibility_pts is not None:
            tested_pts = predetermined_audibility_pts
        else:
            # start with empty tested space
            tested_pts = gpd.GeoDataFrame(columns=['audible', 'geometry'], geometry='geometry', crs=crs)

        study_area = study_area.to_crs(crs)
        valid_query_region = study_area.to_crs(crs).union_all().buffer(-100)  # require points to not be right on the boundary
        study_area_extent = ([study_area.total_bounds[0], study_area.total_bounds[2]],  # ([minx, maxx],
                             [study_area.total_bounds[1], study_area.total_bounds[3]])  # [miny, maxy])

        mic = self._receiver_microphone(study_area, mic, name, crs)

        if self._site_context is None:
            self._prepare_site_context(study_area, mic)

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
            else:
                source_pts = self._source_pts_near_audibility_boundary(
                    tested_pts, coarse_grid, fine_grid, study_area_extent, src_pt_density,
                )

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
            source_pts = contour_active_space(tested_pts, altitude_m)
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

        active_space = build_active_space(tested_pts, crs)
        active_space['altitude_m'] = altitude_m
        active_space['mic_name'] = mic.name

        return active_space, tested_pts

    def _run_propagation_model(
        self,
        job_name: str,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        heading: int | None = None,
    ) -> gpd.GeoDataFrame:
        """
        Execute a single propagation model prediction run.

        Parameters
        ----------
        job_name : str
            Name of this propagation model run to use as a suffix for input and output files.
        source_pts : gpd.GeoDataFrame
            GeoDataFrame of points to test the audibility of with the propagation model.
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
            A GeoDataFrame of points tested during the propagation model run and their audibility.
        """
        assert len(source_pts) > 0, "Trying to run propagation model on zero source points"
        source_pts = source_pts.drop_duplicates("geometry")  # duplicates cause problems when storing prev predictions
        crs = source_pts.crs.to_string().lower()

        # Mark any underground points as inaudible and don't pass them to the propagation model
        aboveground_pts, underground_pts = self._determine_underground_pts(source_pts)

        # AAM uses its own ELV grid; GDAL DEM samples can disagree and break whole batches.
        aboveground_pts, below_model_pts = self.propagation_model.filter_below_terrain(
            self._site_context,
            aboveground_pts,
            job_name=job_name,
        )
        if len(below_model_pts) > 0:
            underground_pts = pd.concat(
                [underground_pts, self._as_inaudible(below_model_pts)],
                ignore_index=True,
            )

        audibility_pts = self._as_inaudible(underground_pts)
        if len(aboveground_pts) == 0:
            return audibility_pts

        # Cache lookup / predict / rewrite lives in prediction_cache; this method
        # only filters geometry and turns spectra into audibility.
        omni_str = os.path.splitext(os.path.basename(omni_source))[0]
        csv_filename = prediction_cache_csv_path(
            self.root_dir,
            self.propagation_model.predictions_subdir,
            altitude_m,
            omni_str,
            heading,
        )
        pred_df, failed_pts = predict_with_cache(
            lambda pts: self.propagation_model.predict(
                self._site_context,
                pts,
                omni_source,
                altitude_m,
                job_name,
                heading,
            ),
            aboveground_pts,
            csv_filename,
            altitude_m,
            job_name,
        )
        if len(failed_pts) > 0:
            audibility_pts = pd.concat(
                [audibility_pts, self._as_inaudible(failed_pts)],
                ignore_index=True,
            )

        predicted_audibility_pts = audible_points_gdf(pred_df, self.ambience, crs)
        audibility_pts = pd.concat([audibility_pts, predicted_audibility_pts], ignore_index=True)

        return audibility_pts

    def _preprocess_source_points(
        self,
        source_pts: gpd.GeoDataFrame,
        valid_query_region: gpd.GeoDataFrame,
        tested_pts: gpd.GeoDataFrame,
        max_pts: int,
    ) -> gpd.GeoDataFrame:
        """
        Filter a set of source points, removing points that:
        1) Are outside the valid query region - meaning outside the study area, or too close
           to the study area boundary (since DEM info isn't present outside the study area)
        2) We already know are audible / inaudible
        3) More than ``propagation_model.max_points_per_run`` allows in one batch

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
            tested already by the propagation model. There is no need to retest a point that has already been tested.
        max_pts: int
            Maximum source points per batch (``propagation_model.max_points_per_run``). If len(source_pts)
            exceeds this, points will be dropped at random (but with a fixed seed) to match this.
            NMSim tolerates up to ~4000 per trajectory; AAM ``ONE TRACK`` is capped at 400.
            If too many points are given in one batch, the model may fail silently and leave prediction output blank.

        Returns
        -------
        source_pts: gpd.GeoDataFrame
            A copy of source_pts filtered appropriately.
        """
        # only query points inside the study area and far enough from the study area boundary
        source_pts = source_pts[source_pts.within(valid_query_region)]

        # don't query points we already know the answer for
        source_pts = source_pts[~source_pts.geometry.isin(tested_pts.geometry)]

        # if too many points for the propagation model, randomly downsample
        if source_pts.shape[0] > max_pts:
            source_pts = source_pts.sample(max_pts, random_state=5)

        return source_pts

    def _determine_underground_pts(
        self, source_pts: gpd.GeoDataFrame,
    ) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
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
    def _source_pts_near_audibility_boundary(
        tested_pts: gpd.GeoDataFrame,
        coarse_grid: gpd.GeoDataFrame,
        fine_grid: gpd.GeoDataFrame,
        study_area_extent: tuple[list[float], list[float]],
        src_pt_density: int,
    ) -> gpd.GeoDataFrame:
        """Fine-grid points near the audible/inaudible boundary of the coarse mesh."""
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
        return fine_grid[fine_grid.within(boundary_zone)]

    @staticmethod
    def _as_inaudible(pts: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """Copy of ``pts`` with ``audible=0``, or ``pts`` unchanged if empty."""
        if pts.empty:
            return pts
        inaudible = pts.copy()
        inaudible["audible"] = 0
        return inaudible

    def _prepare_site_context(self, study_area: gpd.GeoDataFrame, mic: Microphone) -> None:
        """Run model ``prepare_site`` and cache the result on ``_site_context``."""
        self._site_context = self.propagation_model.prepare_site(
            self._resolve_site_elevation(),
            study_area,
            mic,
            project_dem=False,
            suffix=f"_{mic.name}",
        )

    @staticmethod
    def _receiver_microphone(
        study_area: gpd.GeoDataFrame,
        mic: Microphone | None,
        name: str,
        crs: str,
    ) -> Microphone:
        """Return ``mic`` in ``crs``, or a 1.60 m centroid microphone if none was given."""
        if mic:
            mic.to_crs(crs, inplace=True)
            return mic
        study_area_wgs84 = study_area.to_crs('epsg:4326')
        return Microphone(
            name=f"centroid{name}",
            lat=study_area_wgs84.centroid.iat[0].y,
            lon=study_area_wgs84.centroid.iat[0].x,
            z=1.60,  # m, average height of human ear
            crs=crs
        )

    def _resolve_site_elevation(self) -> str:
        """Return the project_setup GeoTIFF path, caching tif/flt paths when needed."""
        if self._dem_file is not None:
            return self._dem_file
        tif_path, flt_path = get_project_setup_elevation(self.root_dir)
        self._dem_file = str(tif_path)
        self._flt_file = str(flt_path)
        return self._dem_file
