"""AAM terrain cache, DEM preparation, and below-ground point filtering."""

from __future__ import annotations

import json
import os
import time
from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from pyproj import CRS, Transformer

from aam_translator import load_terrain, read_nmbgf_grid, write_terrain
from aam_translator.constants import FT_PER_M
from aam_translator.context import TerrainResult

from nps_active_space.active_space.aam_run_log import aam_log
from nps_active_space.active_space.propagation_model import (
    AAM_INPUT_SUBDIR,
    AAM_INPUT_SUBDIR_LEGACY,
    AAM_TERRAIN_SUBDIR,
)
from nps_active_space.utils.computation import project_raster

AAM_INP_BASENAME = "scenario"
AAM_TERRAIN_CACHE_META = "terrain_cache.json"
AAM_DEFAULT_FLOW_RESISTIVITY = 200.0
AOI_BOUNDS_TOLERANCE_DEG = 1e-4
# AAM bilinear-interpolates ELV and rejects the whole batch when z <= surface.
# Filter uses the same ELV grid; tolerance is float noise only, not a clearance margin.
AAM_BELOW_SURFACE_TOLERANCE_M = 0.01


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


def terrain_grid_summary(terrain: TerrainResult) -> str:
    spec = terrain.spec
    return (
        f"AEQD grid {spec.cell_count_x}×{spec.cell_count_y} cells "
        f"at {spec.cell_dx_m:.1f}×{spec.cell_dy_m:.1f} m "
        f"({spec.cell_count_x * spec.cell_count_y:,} cells)"
    )


def terrain_dir_for_site(root_dir: str | Path, suffix: str) -> Path:
    """Return terrain cache dir, preferring ``Input_Data/aam/terrain/{mic}/``."""
    root = Path(root_dir)
    mic_key = suffix.removeprefix("_") or "default"
    canonical = root / AAM_TERRAIN_SUBDIR / mic_key
    if canonical.is_dir():
        return canonical

    legacy_candidates = (
        root / AAM_INPUT_SUBDIR_LEGACY / "terrain" / mic_key,
        root / AAM_INPUT_SUBDIR_LEGACY / f"terrain{suffix}",
        root / AAM_INPUT_SUBDIR / f"terrain{suffix}",
    )
    for legacy_dir in legacy_candidates:
        if legacy_dir.is_dir():
            return legacy_dir

    canonical.mkdir(parents=True, exist_ok=True)
    return canonical


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

    # Legacy cache (ELV/IMP from a prior run without terrain_cache.json).
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
        f"reusing cached terrain in {terrain_dir} "
        f"(ELV not older than {Path(dem_file).name})",
    )
    return terrain


def resolve_dem_for_aam(
    dem_src: str,
    study_area: gpd.GeoDataFrame,
    root_dir: str | Path,
    *,
    project_dem: bool,
    suffix: str,
) -> str:
    """Return the DEM path AAM should resample into ELV/IMP.

    Callers pass the ``project_setup`` GeoTIFF (already clipped to the study area).
    A full-extent warp is only used if ``project_dem`` is true and CRS still differs.
    """
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
        f"write_terrain -> {terrain_dir}/scenario.elv "
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


@lru_cache(maxsize=16)
def _cached_elv_values(elv_path: str, mtime_ns: int) -> np.ndarray:
    """Load ELV payload; cache keyed by path and mtime."""
    return read_nmbgf_grid(elv_path).values


def _elv_grid_values(elv_path: Path) -> np.ndarray:
    resolved = elv_path.resolve()
    return _cached_elv_values(str(resolved), resolved.stat().st_mtime_ns)


def _model_ij_from_aeqd_m(
    terrain: TerrainResult,
    aeqd_x_m: np.ndarray,
    aeqd_y_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Fractional ELV column/row indices from AEQD plane coordinates."""
    spec = terrain.spec
    col_i = (aeqd_x_m - spec.grid_origin_x_m) / spec.cell_dx_m
    row_j = (aeqd_y_m - spec.grid_origin_y_m) / spec.cell_dy_m
    return col_i, row_j


def _bilinear_sample_grid(
    data: np.ndarray,
    col_i: np.ndarray,
    row_j: np.ndarray,
) -> np.ndarray:
    """Bilinear sample a north-up 2D grid at fractional column/row indices."""
    nrows, ncols = data.shape
    samples = np.full(col_i.shape, np.nan, dtype=np.float64)

    c0 = np.floor(col_i).astype(np.int64)
    r0 = np.floor(row_j).astype(np.int64)
    dc = col_i - c0
    dr = row_j - r0
    valid = (c0 >= 0) & (r0 >= 0) & (c0 < ncols - 1) & (r0 < nrows - 1)
    if not np.any(valid):
        return samples

    c0v = c0[valid]
    r0v = r0[valid]
    c1v = c0v + 1
    r1v = r0v + 1
    dcv = dc[valid]
    drv = dr[valid]

    v00 = data[r0v, c0v]
    v01 = data[r0v, c1v]
    v10 = data[r1v, c0v]
    v11 = data[r1v, c1v]
    samples[valid] = (
        (1.0 - drv) * (1.0 - dcv) * v00
        + (1.0 - drv) * dcv * v01
        + drv * (1.0 - dcv) * v10
        + drv * dcv * v11
    )
    return samples


def _terrain_surface_elevation_m(
    source_pts: gpd.GeoDataFrame,
    terrain: TerrainResult,
) -> np.ndarray:
    """Sample AAM terrain MSL (meters) at each source point from the ELV grid."""
    elv_path = terrain.elv_path
    if not elv_path or not Path(elv_path).is_file():
        raise FileNotFoundError(f"AAM ELV grid missing: {elv_path}")

    values = _elv_grid_values(Path(elv_path))
    to_aeqd = Transformer.from_crs(source_pts.crs, terrain.aeqd_crs, always_xy=True)
    xs = source_pts.geometry.x.to_numpy(dtype=np.float64)
    ys = source_pts.geometry.y.to_numpy(dtype=np.float64)
    if xs.size == 1:
        aeqd_x_m, aeqd_y_m = to_aeqd.transform(float(xs[0]), float(ys[0]))
        aeqd_x_m = np.asarray([aeqd_x_m], dtype=np.float64)
        aeqd_y_m = np.asarray([aeqd_y_m], dtype=np.float64)
    else:
        aeqd_x_m, aeqd_y_m = to_aeqd.transform(xs, ys)
    col_i, row_j = _model_ij_from_aeqd_m(terrain, aeqd_x_m, aeqd_y_m)
    raw = _bilinear_sample_grid(values, col_i, row_j)
    if terrain.elv_header_feet:
        return raw / FT_PER_M
    return raw


def split_below_aam_terrain(
    terrain: TerrainResult,
    source_pts: gpd.GeoDataFrame,
    *,
    job_name: str = "",
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Split points that are at or below the AAM ELV terrain surface.

    Samples the same ``.ELV`` grid AAM uses (bilinear, model indices). AAM aborts
    the whole batch when any track point is below that surface; filter them here
    instead of relying on clip-TIF resampling or a large clearance margin.
    """
    if len(source_pts) == 0:
        return source_pts, source_pts.iloc[0:0]

    surface_m = _terrain_surface_elevation_m(source_pts, terrain)
    z_m = source_pts.geometry.z.to_numpy()
    agl_m = z_m - surface_m
    invalid = np.isnan(surface_m)
    below = invalid | (agl_m <= AAM_BELOW_SURFACE_TOLERANCE_M)

    n_below = int(below.sum())
    if n_below > 0:
        deficits_m = -agl_m[below & ~invalid]
        label = f"{job_name}: " if job_name else ""
        if len(deficits_m) > 0:
            aam_log(
                "filter",
                f"{label}filtered {n_below}/{len(source_pts)} below AAM terrain "
                f"(deficit min={deficits_m.min():.2f}m "
                f"max={deficits_m.max():.2f}m mean={deficits_m.mean():.2f}m)",
            )
        else:
            aam_log(
                "filter",
                f"{label}filtered {n_below}/{len(source_pts)} "
                "(no AAM terrain sample / nodata)",
            )

    above = source_pts.loc[~below].copy()
    below_pts = source_pts.loc[below].copy()
    return above, below_pts


def _hop_segment_below_terrain(
    terrain: TerrainResult,
    start_pt,
    end_pt,
    source_crs,
    to_aeqd: Transformer,
    from_aeqd: Transformer,
) -> bool:
    """True if the 3D hop between two vertices intersects the ELV surface."""
    spec = terrain.spec
    step_m = max(min(spec.cell_dx_m, spec.cell_dy_m) / 2.0, 1.0)
    ax0, ay0 = to_aeqd.transform(float(start_pt.x), float(start_pt.y))
    ax1, ay1 = to_aeqd.transform(float(end_pt.x), float(end_pt.y))
    dist_m = float(np.hypot(ax1 - ax0, ay1 - ay0))
    n = max(3, int(np.ceil(dist_m / step_m)) + 1)
    t = np.linspace(0.0, 1.0, n)
    ax = ax0 + t * (ax1 - ax0)
    ay = ay0 + t * (ay1 - ay0)
    xs, ys = from_aeqd.transform(ax, ay)
    z0 = float(start_pt.z)
    z1 = float(end_pt.z)
    zs = z0 + t * (z1 - z0)
    samples = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(xs, ys, zs),
        crs=source_crs,
    )
    surface_m = _terrain_surface_elevation_m(samples, terrain)
    agl_m = zs - surface_m
    invalid = np.isnan(surface_m)
    # Vertices are already ELV-filtered; only the interpolated interior can clip.
    interior = slice(1, -1)
    return bool(np.any(invalid[interior] | (agl_m[interior] <= AAM_BELOW_SURFACE_TOLERANCE_M)))


def split_safe_aam_track_runs(
    terrain: TerrainResult,
    source_pts: gpd.GeoDataFrame,
    *,
    job_name: str = "",
) -> list[gpd.GeoDataFrame]:
    """Break an ordered mesh into ``ONE TRACK`` runs whose hops stay above ELV.

    AAM interpolates between consecutive vertices and aborts the whole track if
    any interpolated sample is below ground. Vertex filtering is not enough.
    """
    if len(source_pts) == 0:
        return []
    if len(source_pts) == 1:
        return [source_pts]

    to_aeqd = Transformer.from_crs(source_pts.crs, terrain.aeqd_crs, always_xy=True)
    from_aeqd = Transformer.from_crs(terrain.aeqd_crs, source_pts.crs, always_xy=True)
    runs: list[list[int]] = [[0]]
    n_cut = 0
    for i in range(len(source_pts) - 1):
        start = source_pts.geometry.iloc[i]
        end = source_pts.geometry.iloc[i + 1]
        if _hop_segment_below_terrain(
            terrain, start, end, source_pts.crs, to_aeqd, from_aeqd,
        ):
            n_cut += 1
            runs.append([i + 1])
        else:
            runs[-1].append(i + 1)

    if n_cut > 0:
        label = f"{job_name}: " if job_name else ""
        aam_log(
            "filter",
            f"{label}split ONE TRACK into {len(runs)} runs "
            f"({n_cut} hops below AAM terrain)",
        )
    return [source_pts.iloc[idx].copy() for idx in runs]
