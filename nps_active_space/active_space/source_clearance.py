"""Lift near-surface sources just above DEM/ELV so they are not treated as underground."""

from __future__ import annotations

import geopandas as gpd
import numpy as np

# Stack / hull height: just above local ground or water surface.
SOURCE_SURFACE_CLEARANCE_M = 2.0
# DEM/ELV noise allowed below the requested MSL z before a point is truly underground.
SOURCE_SURFACE_LIFT_MAX_BELOW_M = 1.0


def apply_surface_clearance_m(
    z_m: np.ndarray,
    surface_m: np.ndarray,
    *,
    clearance_m: float = SOURCE_SURFACE_CLEARANCE_M,
    max_below_m: float = SOURCE_SURFACE_LIFT_MAX_BELOW_M,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(z_cleared_m, is_underground)``.

    Airborne points keep their MSL z. Points at or slightly below the local
    surface (sea level, ground, DEM noise) are lifted to ``surface + clearance``.
    Points deeper than ``max_below_m`` stay underground (e.g. land cells on a
    0 m vessel layer).
    """
    z_m = np.asarray(z_m, dtype=float)
    surface_m = np.asarray(surface_m, dtype=float)
    invalid = ~np.isfinite(surface_m)
    agl_m = z_m - surface_m
    is_underground = invalid | (agl_m < -max_below_m)
    z_cleared_m = z_m.copy()
    near_surface = ~is_underground & (agl_m < clearance_m)
    z_cleared_m[near_surface] = surface_m[near_surface] + clearance_m
    return z_cleared_m, is_underground


def with_point_z(points: gpd.GeoDataFrame, z_m: np.ndarray) -> gpd.GeoDataFrame:
    """Copy ``points`` with geometry z replaced (x/y unchanged)."""
    out = points.copy()
    out.geometry = gpd.points_from_xy(
        points.geometry.x.to_numpy(),
        points.geometry.y.to_numpy(),
        np.asarray(z_m, dtype=float),
    )
    return out
