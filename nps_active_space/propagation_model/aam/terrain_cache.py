"""AAM terrain cache, DEM preparation, and ELV/IMP generation."""

from __future__ import annotations

import json
import os
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import rasterio
from pyproj import CRS

from aam_translator import load_terrain, write_terrain
from aam_translator.context import TerrainResult

from nps_active_space.propagation_model.aam.run_log import aam_log
from nps_active_space.utils.computation import project_raster
from nps_active_space.utils.paths import AAM_TERRAIN_SUBDIR, display_path

AAM_INP_BASENAME = "scenario"
AAM_TERRAIN_CACHE_META = "terrain_cache.json"
AAM_DEFAULT_FLOW_RESISTIVITY = 200.0
AOI_BOUNDS_TOLERANCE_DEG = 1e-4


@dataclass(frozen=True)
class AamTerrainParams:
    receiver_agl_m: float
    flow_resistivity: float
    grid_agl_ft: float
    aoi_bounds: tuple[float, float, float, float]

    @classmethod
    def from_receiver_agl(
        cls,
        receiver_agl_m: float,
        aoi_wgs84,
    ) -> AamTerrainParams:
        return cls(
            receiver_agl_m=receiver_agl_m,
            flow_resistivity=AAM_DEFAULT_FLOW_RESISTIVITY,
            grid_agl_ft=receiver_agl_m * 3.28084,
            aoi_bounds=tuple(aoi_wgs84.bounds),
        )


@contextmanager
def timed_terrain_step(label: str):
    start = time.perf_counter()
    try:
        yield
    finally:
        aam_log("terrain", f"{label} ({time.perf_counter() - start:.1f}s)")


def terrain_dir_for_site(root_dir: str | Path, suffix: str) -> Path:
    """Return ``Input_Data/aam/terrain/{mic}/``, creating it when missing."""
    root = Path(root_dir)
    mic_key = suffix.removeprefix("_") or "default"
    terrain_dir = root / AAM_TERRAIN_SUBDIR / mic_key
    terrain_dir.mkdir(parents=True, exist_ok=True)
    return terrain_dir


def resolve_dem_for_aam(
    dem_src: str,
    study_area: gpd.GeoDataFrame,
    root_dir: str | Path,
    *,
    project_dem: bool,
    suffix: str,
) -> str:
    """Return the DEM path AAM should resample into ELV/IMP."""
    aam_log(
        "terrain",
        f"prepare_site: DEM {_dem_raster_summary(dem_src)}; "
        f"study_area CRS={study_area.crs}",
    )
    if not project_dem:
        return dem_src
    if _crs_matches_dem(study_area.crs, dem_src):
        aam_log(
            "terrain",
            "skipping GDAL warp: DEM CRS already matches study_area",
        )
        return dem_src

    elevation_dir = Path(root_dir) / "Input_Data/01_ELEVATION"
    dem_projected = str(elevation_dir / f"elevation_aam{suffix}.tif")
    with timed_terrain_step(
        f"GDAL warp to study_area CRS -> {Path(dem_projected).name}"
    ):
        project_raster(dem_src, dem_projected, study_area.crs)
    aam_log("terrain", f"warped DEM {_dem_raster_summary(dem_projected)}")
    return dem_projected


def ensure_aam_terrain(
    root_dir: str | Path,
    terrain_dir: Path,
    dem_file: str,
    aoi_wgs84,
    params: AamTerrainParams,
) -> TerrainResult:
    """Load cached ELV/IMP or build fresh AAM terrain for the AOI."""
    root = Path(root_dir)
    aam_log("terrain", f"AOI clip envelope (WGS84): {_aoi_bounds_deg(aoi_wgs84)}")

    terrain = _try_load_cached_terrain(terrain_dir, root, dem_file, params)
    if terrain is not None:
        return terrain

    aam_log(
        "terrain",
        f"write_terrain -> {display_path(terrain_dir)}/scenario.elv "
        f"(receiver AGL {params.receiver_agl_m:.2f} m)",
    )
    with timed_terrain_step("write_terrain (AEQD resample + ELV/IMP)"):
        terrain = write_terrain(
            dem_file,
            aoi_wgs84,
            terrain_dir,
            crs_in="EPSG:4326",
            elv_basename="scenario.elv",
            imp_basename="scenario.imp",
            flow_resistivity=params.flow_resistivity,
            grid_agl_ft=params.grid_agl_ft,
        )
    _write_terrain_cache_meta(terrain_dir, root, dem_file=dem_file, params=params)
    return terrain


def log_terrain_summary(terrain: TerrainResult) -> None:
    aam_log("terrain", terrain_grid_summary(terrain))
    if terrain.clip_tif_path:
        aam_log("terrain", f"clip sidecar: {Path(terrain.clip_tif_path).name}")


def terrain_grid_summary(terrain: TerrainResult) -> str:
    spec = terrain.spec
    return (
        f"AEQD grid {spec.cell_count_x}×{spec.cell_count_y} cells "
        f"at {spec.cell_dx_m:.1f}×{spec.cell_dy_m:.1f} m "
        f"({spec.cell_count_x * spec.cell_count_y:,} cells)"
    )


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
    params: AamTerrainParams,
) -> TerrainResult:
    return load_terrain(
        elv_path,
        imp_path=imp_path,
        clip_tif_path=clip_path,
        grid_agl_ft=params.grid_agl_ft,
        flow_resistivity=params.flow_resistivity,
    )


def _dem_cache_rel(root_dir: Path, dem_file: str) -> str:
    """Site-relative DEM path for cache metadata (portable across host vs /repo)."""
    dem = Path(dem_file).resolve()
    try:
        return str(dem.relative_to(root_dir.resolve()))
    except ValueError:
        return dem.name


def _dem_cache_matches(meta_dem_path: str, root_dir: Path, dem_file: str) -> bool:
    expected = _dem_cache_rel(root_dir, dem_file)
    if meta_dem_path == expected:
        return True
    if meta_dem_path == str(Path(dem_file).resolve()):
        return True
    return Path(meta_dem_path).name == Path(dem_file).name


def _write_terrain_cache_meta(
    terrain_dir: Path,
    root_dir: Path,
    *,
    dem_file: str,
    params: AamTerrainParams,
) -> None:
    meta = {
        "dem_path": _dem_cache_rel(root_dir, dem_file),
        "dem_mtime": os.path.getmtime(dem_file),
        "receiver_agl_m": params.receiver_agl_m,
        "flow_resistivity": params.flow_resistivity,
        "grid_agl_ft": params.grid_agl_ft,
        "aoi_bounds_wgs84": list(params.aoi_bounds),
    }
    (terrain_dir / AAM_TERRAIN_CACHE_META).write_text(
        json.dumps(meta, indent=2),
        encoding="utf-8",
    )


def _terrain_cache_valid(
    terrain_dir: Path,
    root_dir: Path,
    dem_file: str,
    params: AamTerrainParams,
) -> bool:
    elv_path, imp_path, clip_path = _terrain_artifact_paths(terrain_dir)
    if not elv_path.is_file() or not imp_path.is_file():
        return False

    if not Path(dem_file).is_file():
        return False
    if os.path.getmtime(dem_file) > os.path.getmtime(elv_path):
        return False

    meta_path = terrain_dir / AAM_TERRAIN_CACHE_META
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if not _dem_cache_matches(meta.get("dem_path", ""), root_dir, dem_file):
            return False
        if meta.get("dem_mtime", 0) < os.path.getmtime(dem_file):
            return False
        if meta.get("receiver_agl_m") != params.receiver_agl_m:
            return False
        if meta.get("flow_resistivity") != params.flow_resistivity:
            return False
        if meta.get("grid_agl_ft") != params.grid_agl_ft:
            return False
        cached_bounds = tuple(meta.get("aoi_bounds_wgs84", []))
        if len(cached_bounds) != 4 or not _bounds_close(
            params.aoi_bounds, cached_bounds,
        ):
            return False
        return True

    return clip_path.is_file()


def _try_load_cached_terrain(
    terrain_dir: Path,
    root_dir: Path,
    dem_file: str,
    params: AamTerrainParams,
) -> TerrainResult | None:
    if not _terrain_cache_valid(terrain_dir, root_dir, dem_file, params):
        return None

    elv_path, imp_path, clip_path = _terrain_artifact_paths(terrain_dir)
    with timed_terrain_step("load cached terrain from disk"):
        terrain = _terrain_from_disk(
            elv_path,
            imp_path,
            clip_path,
            params=params,
        )
    aam_log(
        "terrain",
        f"reusing cached terrain in {display_path(terrain_dir)} "
        f"(ELV not older than {Path(dem_file).name})",
    )
    return terrain
