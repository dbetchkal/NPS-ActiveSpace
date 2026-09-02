"""AAM implementation of :class:`PropagationModel`."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import dataclass
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

from nps_active_space.active_space.aam_run_log import (
    aam_log,
    configure_aam_run_log,
    log_run_batch,
    summarize_aam_cli_output,
    summarize_aam_error,
)
from nps_active_space.active_space.aam_output import (
    aam_source_id_from_omni,
    apply_omni_gain_offset,
    poi_history_to_predictions_df,
)
from nps_active_space.active_space.aam_terrain import (
    AAM_INP_BASENAME,
    AamTerrainParams,
    ensure_aam_terrain,
    log_terrain_summary,
    resolve_dem_for_aam,
    split_below_aam_terrain,
    terrain_dir_for_site,
)
from nps_active_space.active_space.propagation_model import (
    AAM_PREDICTIONS_SUBDIR,
    AAM_RUNS_SUBDIR,
)
from nps_active_space.utils.models import Microphone

AAM_RUN_TIMEOUT_S = 600
DEFAULT_AAM_CHUNK_SIZE = 50


def _resolve_aam_ncfiles_dir(aam_exe: str | Path) -> Path:
    """Locate AAM NetCDF noise database for a native Windows ``.exe`` launch."""
    override = os.environ.get("AAM_NC", "").strip()
    if override:
        nc_root = Path(override)
        if nc_root.is_dir():
            return nc_root
        raise FileNotFoundError(
            f"AAM_NC is set but not a directory: {nc_root}",
        )

    exe = Path(aam_exe)
    candidates = [
        exe.parent / "NCfiles",
        exe.parent.parent / "NCfiles",
    ]
    for nc_root in candidates:
        if nc_root.is_dir():
            return nc_root

    tried = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"AAM NCfiles/ not found for {exe}; tried {tried}. "
        "Set AAM_NC to the NCfiles directory, or place NCfiles next to the exe "
        "(typical layouts: ...\\Bin\\NCfiles or ...\\AAM\\NCfiles).",
    )


def _aam_subprocess_env(aam_exe: str | Path) -> dict[str, str]:
    """Env for one AAM subprocess. Native Windows exe needs NCfiles next to it."""
    env = os.environ.copy()
    exe = Path(aam_exe)
    if exe.suffix.lower() != ".exe":
        return env
    nc_root = _resolve_aam_ncfiles_dir(exe)
    nc = str(nc_root.resolve()) + os.sep
    env["ROTOR_NOISE"] = nc
    env["FWING_NOISE"] = nc
    env["QUARRY_NOISE"] = nc
    return env

__all__ = [
    "AAM_PREDICTIONS_SUBDIR",
    "AamPropagationModel",
    "AamSiteContext",
    "aam_source_id_from_omni",
    "poi_history_to_predictions_df",
    "split_below_aam_terrain",
]


def _runs_dir_for_site(root_dir: str | Path) -> Path:
    runs_dir = Path(root_dir) / AAM_RUNS_SUBDIR
    runs_dir.mkdir(parents=True, exist_ok=True)
    return runs_dir


def resolve_aam_chunk_size() -> int:
    """Points per Wine/AAM invocation (see docs/aam_integration_notes.md)."""
    return max(1, int(os.environ.get("AAM_CHUNK_SIZE", str(DEFAULT_AAM_CHUNK_SIZE))))


def _order_source_pts_for_track(source_pts: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Sort mesh points for ``ONE TRACK`` so consecutive hops stay local."""
    if len(source_pts) <= 1:
        return source_pts
    ordered = source_pts.copy()
    ordered["_sort_x"] = ordered.geometry.x
    ordered["_sort_y"] = ordered.geometry.y
    return ordered.sort_values(["_sort_x", "_sort_y"]).drop(
        columns=["_sort_x", "_sort_y"],
    )


@dataclass(frozen=True)
class AamSiteContext:
    terrain: TerrainResult
    dem_file: str
    mic: Microphone


class AamPropagationModel:
    max_points_per_run = 400
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
        self._runs_dir = _runs_dir_for_site(root_dir)
        configure_aam_run_log(self._root)

    def __setstate__(self, state: dict) -> None:
        self.__dict__.update(state)
        configure_aam_run_log(self._root)

    def _receiver_agl_m(self, mic: Microphone) -> float:
        if self.receiver_agl_m is not None:
            return self.receiver_agl_m
        return mic.z

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
        chunk_size = resolve_aam_chunk_size()
        if len(source_pts) <= chunk_size:
            return self._predict_batch_resilient(
                site, source_pts, omni_source, altitude_m, job_name, heading,
            )

        frames: list[pd.DataFrame] = []
        ordered = _order_source_pts_for_track(source_pts)
        for chunk_idx, start in enumerate(range(0, len(ordered), chunk_size)):
            chunk_pts = ordered.iloc[start : start + chunk_size]
            chunk_job = f"{job_name}_c{chunk_idx:03d}"
            chunk_frame = self._predict_batch_resilient(
                site,
                chunk_pts,
                omni_source,
                altitude_m,
                chunk_job,
                heading,
            )
            if len(chunk_frame) > 0:
                frames.append(chunk_frame)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)

    def _predict_batch_resilient(
        self,
        site: AamSiteContext,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None = None,
    ) -> pd.DataFrame:
        """Run one batch; on failure split in half until bad single points are skipped."""
        if len(source_pts) == 0:
            return pd.DataFrame()
        try:
            return self._predict_batch(
                site, source_pts, omni_source, altitude_m, job_name, heading,
            )
        except Exception as exc:
            reason = summarize_aam_error(str(exc))
            if len(source_pts) <= 1:
                aam_log(
                    "predict",
                    f"skipped 1 point in {job_name} ({reason})",
                )
                return pd.DataFrame()

            mid = len(source_pts) // 2
            aam_log(
                "predict",
                f"isolating {job_name} n={len(source_pts)} ({reason})",
            )
            frames: list[pd.DataFrame] = []
            for suffix, sub_pts in (
                ("_L", source_pts.iloc[:mid]),
                ("_R", source_pts.iloc[mid:]),
            ):
                sub_frame = self._predict_batch_resilient(
                    site,
                    sub_pts,
                    omni_source,
                    altitude_m,
                    f"{job_name}{suffix}",
                    heading,
                )
                if len(sub_frame) > 0:
                    frames.append(sub_frame)
            if not frames:
                return pd.DataFrame()
            return pd.concat(frames, ignore_index=True)

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
        work_dir = self._runs_dir / job_name
        work_dir.mkdir(parents=True, exist_ok=True)

        ordered_pts = _order_source_pts_for_track(source_pts)
        above_pts, below_pts = self.filter_below_terrain(
            site, ordered_pts, job_name=job_name,
        )
        if len(below_pts) > 0 and len(above_pts) == 0:
            raise RuntimeError(
                f"all {len(below_pts)} source points below AAM terrain for {job_name}",
            )
        if len(above_pts) == 0:
            return pd.DataFrame()
        ordered_pts = above_pts

        track = self._build_track(ordered_pts)
        pois = self._build_pois(site)
        source_id = aam_source_id_from_omni(omni_source)
        heading_deg = float(heading if heading is not None else 90.0)
        speed_kn = hop_speed_kn(track, site.terrain) if len(track) > 1 else 0.0
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
            self._run_aam(inp_path, work_dir)
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

        assert_track_alignment(
            history=histories[0],
            track=track,
            terrain=site.terrain,
            run_log=run_log,
        )
        frame = poi_history_to_predictions_df(histories[0], source_pts)
        return apply_omni_gain_offset(frame, omni_source)

    def _run_aam(self, inp_path: Path, work_dir: Path) -> None:
        if not os.path.isfile(self.aam_shim):
            raise FileNotFoundError(
                f"AAM executable not found at {self.aam_shim}; "
                "set [project] aam in your .config (path to AAM_3.0.0.exe, "
                "or /usr/local/bin/aam in Docker)",
            )
        proc = subprocess.run(
            [self.aam_shim, inp_path.name],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=AAM_RUN_TIMEOUT_S,
            env=_aam_subprocess_env(self.aam_shim),
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
