"""ELV grid sampling, below-ground filtering, and clear-hop track packing."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
from pyproj import Transformer

from aam_translator import read_nmbgf_grid
from aam_translator.constants import FT_PER_M
from aam_translator.context import TerrainResult

from nps_active_space.propagation_model.aam.run_log import aam_log
from nps_active_space.utils.paths import display_path

# AAM bilinear-samples ELV under the flight track (ONE TRACK vertices and the
# linear hops between them; SW-origin j) and aborts that track if vehicle z is
# below ground. POIs are not checked. read_nmbgf_grid is north-up; convert
# model row_j before sampling. Tolerance is float noise only, not a clearance.
AAM_BELOW_SURFACE_TOLERANCE_M = 0.01


def split_below_aam_terrain(
    terrain: TerrainResult,
    source_pts: gpd.GeoDataFrame,
    *,
    job_name: str = "",
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Split points that are at or below the AAM ELV terrain surface."""
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


def split_safe_aam_track_runs(
    terrain: TerrainResult,
    source_pts: gpd.GeoDataFrame,
    *,
    job_name: str = "",
) -> list[gpd.GeoDataFrame]:
    """Pack points into ``ONE TRACK`` runs whose hops stay above ELV."""
    if len(source_pts) == 0:
        return []
    if len(source_pts) == 1:
        return [source_pts]

    to_aeqd = Transformer.from_crs(source_pts.crs, terrain.aeqd_crs, always_xy=True)
    from_aeqd = Transformer.from_crs(terrain.aeqd_crs, source_pts.crs, always_xy=True)
    index_runs = _pack_clear_hop_runs(terrain, source_pts, to_aeqd, from_aeqd)
    if len(index_runs) > 1:
        label = f"{job_name}: " if job_name else ""
        aam_log(
            "filter",
            f"{label}packed {len(source_pts)} points into {len(index_runs)} "
            "AAM tracks (clear hops)",
        )
    return [source_pts.iloc[idx].copy() for idx in index_runs]


def _split_sequential_hop_runs(
    terrain: TerrainResult,
    source_pts: gpd.GeoDataFrame,
) -> list[list[int]]:
    """Cut the given order at every clipping hop (legacy snake split)."""
    if len(source_pts) == 0:
        return []
    if len(source_pts) == 1:
        return [[0]]

    to_aeqd = Transformer.from_crs(source_pts.crs, terrain.aeqd_crs, always_xy=True)
    from_aeqd = Transformer.from_crs(terrain.aeqd_crs, source_pts.crs, always_xy=True)
    geoms = source_pts.geometry
    runs: list[list[int]] = [[0]]
    for i in range(len(source_pts) - 1):
        if _hop_segment_below_terrain(
            terrain, geoms.iloc[i], geoms.iloc[i + 1],
            source_pts.crs, to_aeqd, from_aeqd,
        ):
            runs.append([i + 1])
        else:
            runs[-1].append(i + 1)
    return runs


def _pack_clear_hop_runs(
    terrain: TerrainResult,
    source_pts: gpd.GeoDataFrame,
    to_aeqd: Transformer,
    from_aeqd: Transformer,
) -> list[list[int]]:
    """Greedy nearest-clear-neighbor path cover over vertex-filtered points."""
    n = len(source_pts)
    geoms = source_pts.geometry
    xs = geoms.x.to_numpy(dtype=np.float64)
    ys = geoms.y.to_numpy(dtype=np.float64)
    if n == 1:
        aeqd_x, aeqd_y = to_aeqd.transform(float(xs[0]), float(ys[0]))
        xy = np.array([[aeqd_x, aeqd_y]], dtype=np.float64)
    else:
        aeqd_x, aeqd_y = to_aeqd.transform(xs, ys)
        xy = np.column_stack(
            [np.asarray(aeqd_x, dtype=np.float64), np.asarray(aeqd_y, dtype=np.float64)]
        )

    hop_clips: dict[tuple[int, int], bool] = {}

    def clips(i: int, j: int) -> bool:
        key = (i, j) if i < j else (j, i)
        cached = hop_clips.get(key)
        if cached is None:
            cached = _hop_segment_below_terrain(
                terrain, geoms.iloc[i], geoms.iloc[j],
                source_pts.crs, to_aeqd, from_aeqd,
            )
            hop_clips[key] = cached
        return cached

    unused = set(range(n))
    runs: list[list[int]] = []
    while unused:
        start = min(unused)
        track = [start]
        unused.remove(start)
        while unused:
            current = track[-1]
            remaining = np.fromiter(unused, dtype=np.int64)
            dist_m = np.hypot(
                xy[remaining, 0] - xy[current, 0],
                xy[remaining, 1] - xy[current, 1],
            )
            nearest = remaining[np.lexsort((remaining, dist_m))]
            pick: int | None = None
            for cand in nearest:
                cand_i = int(cand)
                if not clips(current, cand_i):
                    pick = cand_i
                    break
            if pick is None:
                break
            track.append(pick)
            unused.remove(pick)
        runs.append(track)
    return runs


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


def _northup_row_from_model_j(row_south: np.ndarray, nrows: int) -> np.ndarray:
    """Map SW-origin model ``row_j`` onto a north-up ELV array (row 0 = north)."""
    return (nrows - 1) - row_south


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
        raise FileNotFoundError(f"AAM ELV grid missing: {display_path(elv_path)}")

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
    col_i, row_south = _model_ij_from_aeqd_m(terrain, aeqd_x_m, aeqd_y_m)
    row_north = _northup_row_from_model_j(row_south, values.shape[0])
    raw = _bilinear_sample_grid(values, col_i, row_north)
    if terrain.elv_header_feet:
        return raw / FT_PER_M
    return raw


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
    interior = slice(1, -1)
    return bool(np.any(invalid[interior] | (agl_m[interior] <= AAM_BELOW_SURFACE_TOLERANCE_M)))
