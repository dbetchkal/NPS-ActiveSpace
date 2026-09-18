"""AAM implementation of :class:`PropagationModel`."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import dataclass, replace
from pathlib import Path

import geopandas as gpd
import pandas as pd

from aam_translator import (
    assert_track_alignment,
    hop_speed_kn,
    read_poi,
    read_run_log,
    write_inp,
)
from aam_translator.context import TerrainResult
from aam_translator.write_inp import PoiPoint, TrackPoint

from nps_active_space.propagation_model.aam.output import poi_history_to_predictions_df
from nps_active_space.propagation_model.aam.run_log import (
    aam_log,
    configure_aam_run_log,
    FORTRAN_FPA_SUBSCRIPT_ERROR,
    is_fortran_fpa_subscript_error,
    log_run_batch,
    short_aam_work_dir_name,
    summarize_aam_cli_output,
    summarize_aam_error,
)
from nps_active_space.propagation_model.aam.source import (
    aam_subprocess_env,
    aam_template_nc_path,
    ensure_aam_nc_for_source,
    stage_run_ncfiles,
)
from nps_active_space.propagation_model.aam.terrain import (
    AAM_INP_BASENAME,
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
from nps_active_space.utils.paths import AAM_PREDICTIONS_SUBDIR, AAM_RUNS_SUBDIR, display_path

AAM_RUN_TIMEOUT_S = 600
DEFAULT_AAM_CHUNK_SIZE = 400


def resolve_aam_chunk_size() -> int:
    """Points per AAM process (Windows native or Docker+Wine).

    Default is ``DEFAULT_AAM_CHUNK_SIZE`` (400), AAM's ``ONE TRACK`` cap /
    ``aam_translator.MAX_TRACK_POINTS``. Override with env ``AAM_CHUNK_SIZE``.
    """
    return max(1, int(os.environ.get("AAM_CHUNK_SIZE", str(DEFAULT_AAM_CHUNK_SIZE))))


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

        Avoids AAM below-ground aborts rather than retrying: filter vertices
        against the ELV grid, snake the lattice, pack clear-hop tracks, then
        chunk at ``resolve_aam_chunk_size()``. Pad a leftover singleton (~1 m).
        Skip a failed below-ground chunk (do not bisect). Fortran FPA-bounds
        errors on a long high-altitude track may be retried by halving via
        ``_predict_batch_with_fpa_split``. When every batch fails, returns an
        empty frame so the caller can mark those points inaudible.
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
        """Run one batch; if AAM exits with a Fortran subscript error on internal array FPA (see ``summarize_aam_error``), halve the track and retry. Not used for below-ground failures."""
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
        start = time.perf_counter()
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
        ordered_pts = above_pts

        track = pad_single_point_track(self._build_track(ordered_pts))
        pois = self._build_pois(site)
        template_nc = aam_template_nc_path(self.aam_shim)
        source_id, cached_nc = ensure_aam_nc_for_source(
            omni_source,
            self.root_dir,
            template_nc,
        )
        run_nc_dir = stage_run_ncfiles(work_dir, cached_nc)
        heading_deg = float(heading if heading is not None else 90.0)
        speed_kn = hop_speed_kn(track, site.terrain)
        inp_path = work_dir / f"{AAM_INP_BASENAME}.inp"
        aam_log_path = work_dir / f"{AAM_INP_BASENAME}.txt"

        try:
            self._stage_run_dir(
                work_dir,
                site,
                track,
                pois,
                source_id,
                job_name,
                heading_deg,
                speed_kn,
            )
            self._run_aam(inp_path, work_dir, run_nc_dir)
            frame = self._read_run_predictions(
                work_dir, site, track, ordered_pts, omni_source, job_name,
            )
        except Exception as exc:
            log_run_batch(
                self._root,
                job_name=job_name,
                n_track=len(track),
                source_id=source_id,
                heading_deg=heading_deg,
                speed_kn=speed_kn,
                elapsed_s=time.perf_counter() - start,
                inp_path=inp_path,
                ok=False,
                error=str(exc),
                to_console=False,
            )
            raise

        log_run_batch(
            self._root,
            job_name=job_name,
            n_track=len(track),
            source_id=source_id,
            heading_deg=heading_deg,
            speed_kn=speed_kn,
            elapsed_s=time.perf_counter() - start,
            inp_path=inp_path,
            aam_log_path=aam_log_path if aam_log_path.is_file() else None,
            ok=True,
            to_console=False,
        )
        return frame

    def _build_track(self, source_pts: gpd.GeoDataFrame) -> list[TrackPoint]:
        wgs84_pts = source_pts.to_crs("EPSG:4326")
        return [
            TrackPoint(
                lon=float(row.geometry.x),
                lat=float(row.geometry.y),
                alt_m=float(row.geometry.z),
            )
            for _, row in wgs84_pts.iterrows()
        ]

    def _build_pois(self, site: AamSiteContext) -> list[PoiPoint]:
        receiver_agl_m = self._receiver_agl_m(site.mic)
        return [
            PoiPoint(
                name=site.mic.name,
                lon=float(site.mic.lon),
                lat=float(site.mic.lat),
                agl_m=receiver_agl_m,
            ),
        ]

    def _stage_run_dir(
        self,
        work_dir: Path,
        site: AamSiteContext,
        track: list[TrackPoint],
        pois: list[PoiPoint],
        source_id: str,
        job_name: str,
        heading_deg: float,
        speed_kn: float,
    ) -> None:
        shutil.copy2(site.terrain.elv_path, work_dir / "scenario.elv")
        if site.terrain.imp_path:
            shutil.copy2(site.terrain.imp_path, work_dir / "scenario.imp")

        inp_path = work_dir / f"{AAM_INP_BASENAME}.inp"
        write_inp(
            site.terrain,
            inp_path,
            track=track,
            pois=pois,
            source_id=source_id,
            track_name=job_name,
            speed_kn=speed_kn,
            heading_deg=heading_deg,
            elv_basename="scenario.elv",
            imp_basename="scenario.imp",
            remark=f"ActiveSpace {job_name}",
        )

    def _read_run_predictions(
        self,
        work_dir: Path,
        site: AamSiteContext,
        track: list[TrackPoint],
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        job_name: str,
    ) -> pd.DataFrame:
        poi_path = work_dir / f"{AAM_INP_BASENAME}.POI"
        log_path = work_dir / f"{AAM_INP_BASENAME}.txt"
        run_log = read_run_log(log_path)
        if not run_log.ok:
            raise RuntimeError(
                f"AAM run failed for {job_name}: {run_log.read_error}",
            )

        histories = read_poi(poi_path)
        if not histories:
            raise RuntimeError(f"AAM produced no POI zones for {job_name}")

        history = histories[0]
        assert_track_alignment(
            history=history,
            track=track,
            terrain=site.terrain,
            run_log=run_log,
        )
        n_real = len(source_pts)
        if history.n_samples > n_real:
            history = replace(
                history,
                time_s=history.time_s[:n_real],
                broadband_db=history.broadband_db[:n_real],
                band_levels_db=history.band_levels_db[:n_real],
            )
        return poi_history_to_predictions_df(history, source_pts)

    def _run_aam(self, inp_path: Path, work_dir: Path, nc_root: Path) -> None:
        if not os.path.isfile(self.aam_shim):
            raise FileNotFoundError(
                f"AAM executable not found at {display_path(self.aam_shim)}; "
                "set [project] aam in your .config (path to AAM_3.0.0.exe, "
                "or /usr/local/bin/aam in Docker)",
            )
        proc = subprocess.run(
            [self.aam_shim, inp_path.name],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=AAM_RUN_TIMEOUT_S,
            env=aam_subprocess_env(self.aam_shim, nc_root),
        )
        combined = "\n".join(
            part for part in (proc.stderr, proc.stdout) if part
        ).strip()
        if combined:
            (work_dir / "aam_stderr.txt").write_text(
                combined, encoding="utf-8", errors="replace",
            )
        if proc.returncode != 0:
            raise RuntimeError(
                f"AAM exited {proc.returncode}: {summarize_aam_cli_output(combined)}",
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
