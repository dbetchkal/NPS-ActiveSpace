"""AAM implementation of :class:`PropagationModel`."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import pandas as pd

from aam_translator import hop_speed_kn
from aam_translator.context import TerrainResult

from nps_active_space.propagation_model.aam.config import resolve_aam_chunk_size
from nps_active_space.propagation_model.aam.run_batch import (
    execute_aam_batch,
    geo_source_to_track,
    mic_to_pois,
)
from nps_active_space.propagation_model.aam.run_log import (
    aam_log,
    configure_aam_run_log,
    FORTRAN_FPA_SUBSCRIPT_ERROR,
    is_fortran_fpa_subscript_error,
    short_aam_work_dir_name,
    summarize_aam_error,
)
from nps_active_space.propagation_model.aam.source import (
    aam_template_nc_path,
    ensure_aam_nc_for_source,
    stage_run_ncfiles,
)
from nps_active_space.propagation_model.aam.terrain import (
    AamTerrainParams,
    ensure_aam_terrain,
    log_terrain_summary,
    resolve_dem_for_aam,
    split_below_aam_terrain,
    split_safe_aam_track_runs,
    terrain_dir_for_site,
)
from nps_active_space.propagation_model.aam.track import (
    order_source_pts_for_track,
    pad_single_point_track,
)
from nps_active_space.propagation_model.protocol import DEFAULT_MAX_POINTS_PER_PREDICT
from nps_active_space.utils.models import Microphone
from nps_active_space.utils.paths import AAM_PREDICTIONS_SUBDIR, AAM_RUNS_SUBDIR


@dataclass(frozen=True)
class AamSiteContext:
    terrain: TerrainResult
    dem_file: str
    mic: Microphone


class AamPropagationModel:
    max_points_per_run = DEFAULT_MAX_POINTS_PER_PREDICT
    predictions_subdir = AAM_PREDICTIONS_SUBDIR

    def __init__(
        self,
        root_dir: str,
        aam_shim: str = "/usr/local/bin/aam",
        *,
        receiver_agl_m: float | None = None,
    ) -> None:
        self.root_dir = root_dir
        self.aam_shim = aam_shim
        self.receiver_agl_m = receiver_agl_m
        self._root = Path(root_dir).resolve()
        self._runs_dir = Path(root_dir) / AAM_RUNS_SUBDIR
        self._runs_dir.mkdir(parents=True, exist_ok=True)
        configure_aam_run_log(self._root)

    def __setstate__(self, state: dict) -> None:
        self.__dict__.update(state)
        configure_aam_run_log(self._root)

    def filter_below_terrain(
        self,
        site: AamSiteContext,
        source_pts: gpd.GeoDataFrame,
        *,
        job_name: str = "",
    ) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """Drop points below the AAM ELV surface before building a track batch."""
        return split_below_aam_terrain(
            site.terrain, source_pts, job_name=job_name,
        )

    def prepare_site(
        self,
        dem_src: str,
        study_area: gpd.GeoDataFrame,
        mic: Microphone,
        *,
        project_dem: bool = True,
        suffix: str = "",
    ) -> AamSiteContext:
        dem_file = resolve_dem_for_aam(
            dem_src,
            study_area,
            self.root_dir,
            project_dem=project_dem,
            suffix=suffix,
        )
        terrain_dir = terrain_dir_for_site(self.root_dir, suffix)

        mic_wgs84 = mic.to_crs("EPSG:4326")
        receiver_agl_m = self._receiver_agl_m(mic_wgs84)
        aoi = study_area.to_crs("EPSG:4326").union_all()
        params = AamTerrainParams.from_receiver_agl(receiver_agl_m, aoi)

        terrain = ensure_aam_terrain(
            self.root_dir,
            terrain_dir,
            dem_file,
            aoi,
            params,
        )
        log_terrain_summary(terrain)

        return AamSiteContext(
            terrain=terrain,
            dem_file=dem_file,
            mic=mic_wgs84,
        )

    def predict(
        self,
        site: AamSiteContext,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None = None,
    ) -> pd.DataFrame:
        """Predict spectra at the mic for each source point.

        Pipeline: snake-order mesh points → ELV vertex filter → pack clear-hop
        ``ONE TRACK`` runs → chunk (``resolve_aam_chunk_size``) → pad singleton
        tracks → subprocess + POI read. Failed below-ground batches are skipped
        (not bisected). Fortran FPA subscript errors may be retried by halving
        the track. An empty frame means the caller should mark points inaudible.
        """
        ordered = order_source_pts_for_track(source_pts)
        above_pts, _below_pts = self.filter_below_terrain(
            site, ordered, job_name=job_name,
        )
        if len(above_pts) == 0:
            return pd.DataFrame()

        frames, run_idx = self._predict_chunked_runs(
            site, above_pts, omni_source, altitude_m, job_name, heading,
        )
        if not frames:
            aam_log(
                "predict",
                f"no predictions for {job_name}: 0/{run_idx} batch(es) "
                f"succeeded for {len(above_pts)} above-ground point(s); "
                "caller will mark inaudible. Inspect "
                f"Output_Data/aam/runs/{job_name}_r*/scenario.txt "
                "and aam_stderr.txt (terrain, NCfiles, below-ground).",
            )
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)

    def _predict_chunked_runs(
        self,
        site: AamSiteContext,
        above_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None,
    ) -> tuple[list[pd.DataFrame], int]:
        chunk_size = resolve_aam_chunk_size()
        frames: list[pd.DataFrame] = []
        runs = split_safe_aam_track_runs(
            site.terrain, above_pts, job_name=job_name,
        )
        run_idx = 0
        for run_pts in runs:
            for start in range(0, len(run_pts), chunk_size):
                chunk_pts = run_pts.iloc[start : start + chunk_size]
                chunk_job = f"{job_name}_r{run_idx:03d}"
                run_idx += 1
                try:
                    chunk_frame = self._predict_batch_with_fpa_split(
                        site,
                        chunk_pts,
                        omni_source,
                        altitude_m,
                        chunk_job,
                        heading,
                    )
                except Exception as exc:
                    aam_log(
                        "predict",
                        f"skipped {chunk_job} n={len(chunk_pts)} "
                        f"({summarize_aam_error(str(exc))})",
                    )
                    continue
                if len(chunk_frame) > 0:
                    frames.append(chunk_frame)
        return frames, run_idx

    def _predict_batch_with_fpa_split(
        self,
        site: AamSiteContext,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None,
        *,
        split_depth: int = 0,
    ) -> pd.DataFrame:
        try:
            return self._predict_batch(
                site,
                source_pts,
                omni_source,
                altitude_m,
                job_name,
                heading,
            )
        except Exception as exc:
            if not is_fortran_fpa_subscript_error(exc) or len(source_pts) <= 2:
                raise
            mid = len(source_pts) // 2
            aam_log(
                "predict",
                f"split {job_name} n={len(source_pts)} after {FORTRAN_FPA_SUBSCRIPT_ERROR}",
            )
            left = self._predict_batch_with_fpa_split(
                site,
                source_pts.iloc[:mid],
                omni_source,
                altitude_m,
                f"{job_name}_sa{split_depth}",
                heading,
                split_depth=split_depth + 1,
            )
            right = self._predict_batch_with_fpa_split(
                site,
                source_pts.iloc[mid:],
                omni_source,
                altitude_m,
                f"{job_name}_sb{split_depth}",
                heading,
                split_depth=split_depth + 1,
            )
            return pd.concat([left, right], ignore_index=True)

    def _predict_batch(
        self,
        site: AamSiteContext,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None = None,
    ) -> pd.DataFrame:
        work_dir_name = short_aam_work_dir_name(job_name)
        work_dir = self._runs_dir / work_dir_name
        work_dir.mkdir(parents=True, exist_ok=True)
        aam_log("run-dir", f"{work_dir_name} <- {job_name}", to_console=False)

        above_pts, below_pts = self.filter_below_terrain(
            site, source_pts, job_name=job_name,
        )
        if len(below_pts) > 0 and len(above_pts) == 0:
            raise RuntimeError(
                f"all {len(below_pts)} source points below AAM terrain for {job_name}",
            )
        if len(above_pts) == 0:
            return pd.DataFrame()

        track = pad_single_point_track(geo_source_to_track(above_pts))
        receiver_agl_m = self._receiver_agl_m(site.mic)
        pois = mic_to_pois(site.mic, receiver_agl_m)
        template_nc = aam_template_nc_path(self.aam_shim)
        source_id, cached_nc = ensure_aam_nc_for_source(
            omni_source,
            self.root_dir,
            template_nc,
        )
        run_nc_dir = stage_run_ncfiles(work_dir, cached_nc)
        heading_deg = float(heading if heading is not None else 90.0)
        speed_kn = hop_speed_kn(track, site.terrain)

        return execute_aam_batch(
            root=self._root,
            aam_shim=self.aam_shim,
            work_dir=work_dir,
            terrain=site.terrain,
            track=track,
            pois=pois,
            source_pts=above_pts,
            run_nc_dir=run_nc_dir,
            source_id=source_id,
            job_name=job_name,
            heading_deg=heading_deg,
            speed_kn=speed_kn,
        )

    def _receiver_agl_m(self, mic: Microphone) -> float:
        if self.receiver_agl_m is not None:
            return self.receiver_agl_m
        return mic.z


__all__ = [
    "AamPropagationModel",
    "AamSiteContext",
    "resolve_aam_chunk_size",
]
