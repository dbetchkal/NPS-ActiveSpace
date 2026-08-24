"""AAM implementation of :class:`PropagationModel`."""

from __future__ import annotations

import logging
import json
import re
import shutil
import os
import subprocess
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import CRS

from aam_translator import (
    assert_track_alignment,
    hop_speed_kn,
    read_poi,
    read_run_log,
    write_inp,
    write_terrain,
)
from aam_translator.bands import band_number_for_frequency
from aam_translator.context import TerrainResult
from aam_translator.grid_spec import GridSpec
from aam_translator.nmbgf_io import read_nmbgf_header
from aam_translator.read_poi import PoiTimeHistory
from aam_translator.write_inp import PoiPoint, TrackPoint

from nps_active_space.active_space.propagation_model import NMSIM_BAND_COLUMNS
from nps_active_space.utils.computation import project_raster
from nps_active_space.utils.helpers import omni_to_gain
from nps_active_space.utils.models import Microphone

logger = logging.getLogger(__name__)

AAM_WORK_SUBDIR = "Input_Data/AAM"
AAM_INP_BASENAME = "scenario"
AAM_RUN_TIMEOUT_S = 600
AAM_ELEVATION_MASK = "elevation_mask.tif"
AAM_TERRAIN_CACHE_META = "terrain_cache.json"
AAM_DEFAULT_FLOW_RESISTIVITY = 200.0
AOI_BOUNDS_TOLERANCE_DEG = 1e-4


def _terrain_log(msg: str) -> None:
    """Visible in Docker validate runs (logging may be unset)."""
    line = f"[aam-terrain] {msg}"
    logger.info(line)
    print(line, flush=True)


@contextmanager
def _timed_terrain_step(label: str):
    _terrain_log(f"{label}...")
    start = time.perf_counter()
    try:
        yield
    finally:
        _terrain_log(f"{label} done ({time.perf_counter() - start:.1f}s)")


def _dem_raster_summary(path: str) -> str:
    with rasterio.open(path) as src:
        res_x, res_y = src.res
        return (
            f"{Path(path).name}: {src.width}×{src.height} px, "
            f"CRS={src.crs}, res≈{abs(res_x):.6f}×{abs(res_y):.6f}"
        )


def _crs_matches_dem(study_area_crs: str | object, dem_path: str) -> bool:
    with rasterio.open(dem_path) as src:
        if src.crs is None:
            return False
        return CRS.from_user_input(study_area_crs).equals(CRS.from_user_input(src.crs))


def _aoi_bounds_deg(aoi_wgs84) -> str:
    xmin, ymin, xmax, ymax = aoi_wgs84.bounds
    return (
        f"lon [{xmin:.5f}, {xmax:.5f}], "
        f"lat [{ymin:.5f}, {ymax:.5f}]"
    )


def _terrain_grid_summary(terrain: TerrainResult) -> str:
    spec = terrain.spec
    return (
        f"AEQD grid {spec.cell_count_x}×{spec.cell_count_y} cells "
        f"at {spec.cell_dx_m:.1f}×{spec.cell_dy_m:.1f} m "
        f"({spec.cell_count_x * spec.cell_count_y:,} cells)"
    )


def _aoi_bounds_tuple(aoi_wgs84) -> tuple[float, float, float, float]:
    return tuple(aoi_wgs84.bounds)


def _bounds_close(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
    tol: float = AOI_BOUNDS_TOLERANCE_DEG,
) -> bool:
    return all(abs(x - y) <= tol for x, y in zip(a, b, strict=True))


def _terrain_artifact_paths(terrain_dir: Path) -> tuple[Path, Path, Path]:
    elv = terrain_dir / f"{AAM_INP_BASENAME}.elv"
    imp = terrain_dir / f"{AAM_INP_BASENAME}.imp"
    clip = terrain_dir / f"{AAM_INP_BASENAME}_clip.tif"
    return elv, imp, clip


def _terrain_from_disk(
    elv_path: Path,
    imp_path: Path,
    clip_path: Path,
    *,
    grid_agl_ft: float,
    flow_resistivity: float,
) -> TerrainResult:
    if not clip_path.is_file():
        raise FileNotFoundError(f"clip GeoTIFF missing: {clip_path}")

    with rasterio.open(clip_path) as src:
        transform = src.transform
        cell_dx_m = abs(transform.a)
        cell_dy_m = abs(transform.e)
        spec = GridSpec(
            cell_count_x=src.width,
            cell_count_y=src.height,
            cell_dx_m=cell_dx_m,
            cell_dy_m=cell_dy_m,
            grid_origin_x_m=transform.c,
            grid_origin_y_m=transform.f + src.height * transform.e,
        )
        aeqd_crs = CRS.from_user_input(src.crs)

    elv_header = read_nmbgf_header(elv_path)
    elv_header_feet = elv_header.units.strip().upper() == "FEET"

    return TerrainResult(
        spec=spec,
        aeqd_crs=aeqd_crs,
        elv_header_feet=elv_header_feet,
        elv_path=str(elv_path),
        imp_path=str(imp_path),
        clip_tif_path=str(clip_path),
        grid_agl_ft=grid_agl_ft,
        flow_resistivity=flow_resistivity,
    )


def _write_terrain_cache_meta(
    terrain_dir: Path,
    *,
    dem_file: str,
    receiver_agl_m: float,
    flow_resistivity: float,
    grid_agl_ft: float,
    aoi_bounds: tuple[float, float, float, float],
) -> None:
    meta = {
        "dem_path": str(Path(dem_file).resolve()),
        "dem_mtime": os.path.getmtime(dem_file),
        "receiver_agl_m": receiver_agl_m,
        "flow_resistivity": flow_resistivity,
        "grid_agl_ft": grid_agl_ft,
        "aoi_bounds_wgs84": list(aoi_bounds),
    }
    (terrain_dir / AAM_TERRAIN_CACHE_META).write_text(
        json.dumps(meta, indent=2),
        encoding="utf-8",
    )


def _terrain_cache_valid(
    terrain_dir: Path,
    dem_file: str,
    receiver_agl_m: float,
    flow_resistivity: float,
    grid_agl_ft: float,
    aoi_bounds: tuple[float, float, float, float],
) -> bool:
    elv_path, imp_path, clip_path = _terrain_artifact_paths(terrain_dir)
    if not elv_path.is_file() or not imp_path.is_file():
        return False

    dem_resolved = str(Path(dem_file).resolve())
    if not Path(dem_file).is_file():
        return False
    if os.path.getmtime(dem_file) > os.path.getmtime(elv_path):
        return False

    meta_path = terrain_dir / AAM_TERRAIN_CACHE_META
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("dem_path") != dem_resolved:
            return False
        if meta.get("dem_mtime", 0) < os.path.getmtime(dem_file):
            return False
        if meta.get("receiver_agl_m") != receiver_agl_m:
            return False
        if meta.get("flow_resistivity") != flow_resistivity:
            return False
        if meta.get("grid_agl_ft") != grid_agl_ft:
            return False
        cached_bounds = tuple(meta.get("aoi_bounds_wgs84", []))
        if len(cached_bounds) != 4 or not _bounds_close(
            aoi_bounds, cached_bounds,
        ):
            return False
        return True

    # Legacy cache (ELV/IMP from a prior run without terrain_cache.json).
    return clip_path.is_file()


def _try_load_cached_terrain(
    terrain_dir: Path,
    dem_file: str,
    aoi_wgs84,
    receiver_agl_m: float,
    flow_resistivity: float,
    grid_agl_ft: float,
) -> TerrainResult | None:
    aoi_bounds = _aoi_bounds_tuple(aoi_wgs84)
    if not _terrain_cache_valid(
        terrain_dir,
        dem_file,
        receiver_agl_m,
        flow_resistivity,
        grid_agl_ft,
        aoi_bounds,
    ):
        return None

    elv_path, imp_path, clip_path = _terrain_artifact_paths(terrain_dir)
    with _timed_terrain_step("load cached terrain from disk"):
        terrain = _terrain_from_disk(
            elv_path,
            imp_path,
            clip_path,
            grid_agl_ft=grid_agl_ft,
            flow_resistivity=flow_resistivity,
        )
    _terrain_log(
        f"reusing cached terrain in {terrain_dir} "
        f"(ELV not older than {Path(dem_file).name})"
    )
    return terrain

# NMSim omni tuning basename -> AAM NetCDF source id (FLATO*.nc in NCfiles/).
OMNI_TO_AAM_SOURCE_ID: dict[str, str] = {
    "O_+200": "FLATO200",
}


def aam_source_id_from_omni(omni_source: str) -> str:
    stem = Path(omni_source).stem
    if stem in OMNI_TO_AAM_SOURCE_ID:
        return OMNI_TO_AAM_SOURCE_ID[stem]
    if stem.startswith("FLATO"):
        return stem
    # NMSim omni tuning files (O_+005 = +0.5 dB, etc.) map to FLATO200; gain applied in predict().
    if re.fullmatch(r"O_[+-]\d{3}", stem):
        return "FLATO200"
    raise ValueError(
        f"no AAM source_id mapping for omni source {omni_source!r}; "
        f"known keys: {sorted(OMNI_TO_AAM_SOURCE_ID)}",
    )


def apply_omni_gain_offset(frame: pd.DataFrame, omni_source: str) -> pd.DataFrame:
    """Shift prediction levels by the NMSim omni tuning gain (dB)."""
    offset_db = omni_to_gain(omni_source)
    if offset_db == 0:
        return frame
    out = frame.copy()
    out["A"] = out["A"] + offset_db
    for col in NMSIM_BAND_COLUMNS:
        if col in out.columns:
            out[col] = out[col] + offset_db
    return out


def _band_number_for_nmsim_column(column: str) -> int:
    return band_number_for_frequency(float(column))



def poi_history_to_predictions_df(
    history: PoiTimeHistory,
    source_pts: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """Map one POI zone (single receiver) to the NMSim TIS-shaped DataFrame."""
    n = history.n_samples
    if n != len(source_pts):
        raise ValueError(
            f"POI has {n} rows but source_pts has {len(source_pts)} points",
        )

    coords = source_pts.to_crs(source_pts.crs)
    frame = pd.DataFrame({
        "Xpos": coords.geometry.x.values,
        "Ypos": coords.geometry.y.values,
        "Zpos": coords.geometry.z.values,
        "A": history.broadband("dBA"),
    })

    band_index = {bn: i for i, bn in enumerate(history.band_numbers)}
    for col in NMSIM_BAND_COLUMNS:
        band_num = _band_number_for_nmsim_column(col)
        if band_num in band_index:
            frame[col] = history.band_levels_db[:, band_index[band_num]]
        else:
            frame[col] = np.nan

    return frame


@dataclass(frozen=True)
class AamSiteContext:
    terrain: TerrainResult
    dem_file: str
    mic: Microphone


class AamPropagationModel:
    max_points_per_run = 400

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
        self._aam_work_dir = Path(root_dir) / AAM_WORK_SUBDIR
        self._aam_work_dir.mkdir(parents=True, exist_ok=True)

    def prepare_site(
        self,
        dem_src: str,
        study_area: gpd.GeoDataFrame,
        mic: Microphone,
        *,
        project_dem: bool = True,
        suffix: str = "",
    ) -> AamSiteContext:
        elevation_dir = Path(self.root_dir) / "Input_Data/01_ELEVATION"
        masked_dem = elevation_dir / AAM_ELEVATION_MASK

        _terrain_log(
            f"prepare_site: parent DEM {_dem_raster_summary(dem_src)}; "
            f"study_area CRS={study_area.crs}"
        )
        if masked_dem.is_file():
            _terrain_log(
                f"NMSim study-area clipped DEM exists "
                f"({_dem_raster_summary(str(masked_dem))}) — "
                "AAM path still uses config [data] dem unless changed"
            )

        dem_file = dem_src
        if project_dem:
            dem_projected = str(elevation_dir / f"elevation_aam{suffix}.tif")
            if _crs_matches_dem(study_area.crs, dem_src):
                _terrain_log(
                    "skipping GDAL warp: DEM CRS already matches study_area "
                    f"(would have written {Path(dem_projected).name})"
                )
            else:
                with _timed_terrain_step(
                    f"GDAL warp to study_area CRS -> {Path(dem_projected).name}"
                ):
                    project_raster(dem_src, dem_projected, study_area.crs)
                dem_file = dem_projected
                _terrain_log(f"warped DEM {_dem_raster_summary(dem_file)}")

        terrain_dir = self._aam_work_dir / f"terrain{suffix}"
        terrain_dir.mkdir(parents=True, exist_ok=True)

        mic_wgs84 = mic.to_crs("EPSG:4326")
        receiver_agl_m = self.receiver_agl_m if self.receiver_agl_m is not None else mic_wgs84.z

        aoi = study_area.to_crs("EPSG:4326").union_all()
        aoi_bounds = _aoi_bounds_tuple(aoi)
        grid_agl_ft = receiver_agl_m * 3.28084
        flow_resistivity = AAM_DEFAULT_FLOW_RESISTIVITY

        _terrain_log(f"AOI clip envelope (WGS84): {_aoi_bounds_deg(aoi)}")

        terrain = _try_load_cached_terrain(
            terrain_dir,
            dem_file,
            aoi,
            receiver_agl_m,
            flow_resistivity,
            grid_agl_ft,
        )
        if terrain is None:
            _terrain_log(
                f"write_terrain -> {terrain_dir}/scenario.elv "
                f"(native DEM posting, receiver AGL {receiver_agl_m:.2f} m)"
            )
            with _timed_terrain_step("write_terrain (AEQD resample + ELV/IMP)"):
                terrain = write_terrain(
                    dem_file,
                    aoi,
                    terrain_dir,
                    crs_in="EPSG:4326",
                    elv_basename="scenario.elv",
                    imp_basename="scenario.imp",
                    flow_resistivity=flow_resistivity,
                    grid_agl_ft=grid_agl_ft,
                )
            _write_terrain_cache_meta(
                terrain_dir,
                dem_file=dem_file,
                receiver_agl_m=receiver_agl_m,
                flow_resistivity=flow_resistivity,
                grid_agl_ft=grid_agl_ft,
                aoi_bounds=aoi_bounds,
            )

        _terrain_log(_terrain_grid_summary(terrain))
        if terrain.clip_tif_path:
            _terrain_log(f"clip sidecar: {Path(terrain.clip_tif_path).name}")

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
        work_dir = self._aam_work_dir / job_name
        work_dir.mkdir(parents=True, exist_ok=True)

        wgs84_pts = source_pts.to_crs("EPSG:4326")
        track = [
            TrackPoint(
                lon=float(row.geometry.x),
                lat=float(row.geometry.y),
                alt_m=float(row.geometry.z),
            )
            for _, row in wgs84_pts.iterrows()
        ]

        receiver_agl_m = (
            self.receiver_agl_m if self.receiver_agl_m is not None else site.mic.z
        )
        pois = [
            PoiPoint(
                name=site.mic.name,
                lon=float(site.mic.lon),
                lat=float(site.mic.lat),
                agl_m=receiver_agl_m,
            ),
        ]


        shutil.copy2(site.terrain.elv_path, work_dir / "scenario.elv")
        if site.terrain.imp_path:
            shutil.copy2(site.terrain.imp_path, work_dir / "scenario.imp")

        speed_kn = hop_speed_kn(track, site.terrain) if len(track) > 1 else 0.0
        inp_path = work_dir / f"{AAM_INP_BASENAME}.inp"
        write_inp(
            site.terrain,
            inp_path,
            track=track,
            pois=pois,
            source_id=aam_source_id_from_omni(omni_source),
            track_name=job_name,
            speed_kn=speed_kn,
            heading_deg=float(heading if heading is not None else 90.0),
            elv_basename="scenario.elv",
            imp_basename="scenario.imp",
            remark=f"ActiveSpace {job_name}",
        )

        self._run_aam(inp_path, work_dir)

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
                f"AAM shim not found at {self.aam_shim}; "
                "use docker/run_activespace.sh -m aam or set aam_shim",
            )
        proc = subprocess.run(
            [self.aam_shim, inp_path.name],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=AAM_RUN_TIMEOUT_S,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"AAM exited {proc.returncode}: {proc.stderr or proc.stdout}",
            )
