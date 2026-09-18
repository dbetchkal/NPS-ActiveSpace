"""Stage one AAM subprocess, run the executable, and read POI predictions."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import replace
from pathlib import Path

import geopandas as gpd
import pandas as pd

from aam_translator import (
    assert_track_alignment,
    read_poi,
    read_run_log,
    write_inp,
)
from aam_translator.context import TerrainResult
from aam_translator.write_inp import PoiPoint, TrackPoint

from nps_active_space.propagation_model.aam.config import AAM_RUN_TIMEOUT_S
from nps_active_space.propagation_model.aam.output import poi_history_to_predictions_df
from nps_active_space.propagation_model.aam.run_log import (
    log_run_batch,
    summarize_aam_cli_output,
)
from nps_active_space.propagation_model.aam.source import aam_subprocess_env
from nps_active_space.propagation_model.aam.terrain_cache import AAM_INP_BASENAME
from nps_active_space.utils.models import Microphone
from nps_active_space.utils.paths import display_path

__all__ = [
    "execute_aam_batch",
    "geo_source_to_track",
    "mic_to_pois",
    "read_aam_batch_predictions",
    "run_aam_subprocess",
    "stage_aam_run_directory",
]


def geo_source_to_track(source_pts: gpd.GeoDataFrame) -> list[TrackPoint]:
    wgs84_pts = source_pts.to_crs("EPSG:4326")
    return [
        TrackPoint(
            lon=float(row.geometry.x),
            lat=float(row.geometry.y),
            alt_m=float(row.geometry.z),
        )
        for _, row in wgs84_pts.iterrows()
    ]


def mic_to_pois(mic: Microphone, receiver_agl_m: float) -> list[PoiPoint]:
    return [
        PoiPoint(
            name=mic.name,
            lon=float(mic.lon),
            lat=float(mic.lat),
            agl_m=receiver_agl_m,
        ),
    ]


def stage_aam_run_directory(
    work_dir: Path,
    terrain: TerrainResult,
    track: list[TrackPoint],
    pois: list[PoiPoint],
    *,
    source_id: str,
    job_name: str,
    heading_deg: float,
    speed_kn: float,
) -> tuple[Path, Path]:
    """Copy ELV/IMP and write ``scenario.inp``; return inp and expected AAM log paths."""
    shutil.copy2(terrain.elv_path, work_dir / "scenario.elv")
    if terrain.imp_path:
        shutil.copy2(terrain.imp_path, work_dir / "scenario.imp")

    inp_path = work_dir / f"{AAM_INP_BASENAME}.inp"
    aam_log_path = work_dir / f"{AAM_INP_BASENAME}.txt"
    write_inp(
        terrain,
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
    return inp_path, aam_log_path


def run_aam_subprocess(
    aam_shim: str,
    inp_path: Path,
    work_dir: Path,
    nc_root: Path,
    *,
    timeout_s: int = AAM_RUN_TIMEOUT_S,
) -> None:
    if not os.path.isfile(aam_shim):
        raise FileNotFoundError(
            f"AAM executable not found at {display_path(aam_shim)}; "
            "set [project] aam in your .config (path to AAM_3.0.0.exe, "
            "or /usr/local/bin/aam in Docker)",
        )
    proc = subprocess.run(
        [aam_shim, inp_path.name],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=timeout_s,
        env=aam_subprocess_env(aam_shim, nc_root),
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


def read_aam_batch_predictions(
    work_dir: Path,
    terrain: TerrainResult,
    track: list[TrackPoint],
    source_pts: gpd.GeoDataFrame,
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
        terrain=terrain,
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


def execute_aam_batch(
    *,
    root: Path,
    aam_shim: str,
    work_dir: Path,
    terrain: TerrainResult,
    track: list[TrackPoint],
    pois: list[PoiPoint],
    source_pts: gpd.GeoDataFrame,
    run_nc_dir: Path,
    source_id: str,
    job_name: str,
    heading_deg: float,
    speed_kn: float,
) -> pd.DataFrame:
    """Run one staged AAM batch and append success/failure lines to the site log."""
    start = time.perf_counter()
    inp_path, aam_log_path = stage_aam_run_directory(
        work_dir,
        terrain,
        track,
        pois,
        source_id=source_id,
        job_name=job_name,
        heading_deg=heading_deg,
        speed_kn=speed_kn,
    )
    try:
        run_aam_subprocess(aam_shim, inp_path, work_dir, run_nc_dir)
        frame = read_aam_batch_predictions(
            work_dir, terrain, track, source_pts, job_name,
        )
    except Exception as exc:
        log_run_batch(
            root,
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
        root,
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
