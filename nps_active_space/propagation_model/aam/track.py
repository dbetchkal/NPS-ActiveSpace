"""Order mesh points and pad tracks for AAM ``ONE TRACK`` batches."""

from __future__ import annotations

import math

import geopandas as gpd
import pandas as pd
from aam_translator.write_inp import TrackPoint

# AAM 3.0.0 aborts when ``ONE TRACK`` has only one vertex (native Windows or Wine;
# Fortran subscript error on internal array FPA; empty ``.POI`` — see
# ``run_log.summarize_aam_error``). Chunking can leave a lone mesh point as a singleton
# track; pad ~1 m east so AAM always gets two vertices. aam-translator notes:
# https://github.com/elliott-ruebush/aam-translator/blob/main/docs/reading_aam_output.md
# https://github.com/elliott-ruebush/aam-translator/blob/main/references/notes/aam_inp_format.md
SINGLE_TRACK_PAD_M = 1.0

__all__ = [
    "SINGLE_TRACK_PAD_M",
    "order_source_pts_for_track",
    "pad_single_point_track",
]


def pad_single_point_track(track: list[TrackPoint]) -> list[TrackPoint]:
    """Duplicate a lone vertex ~1 m east so AAM can interpolate a track."""
    if len(track) != 1:
        return track
    point = track[0]
    cos_lat = math.cos(math.radians(point.lat))
    # Spherical-Earth scale for this pad only (~111.32 km/° lat); not a package geodesy constant.
    meters_per_deg_lon = 111_320.0 * max(abs(cos_lat), 1e-6)
    pad = TrackPoint(
        lon=point.lon + SINGLE_TRACK_PAD_M / meters_per_deg_lon,
        lat=point.lat,
        alt_m=point.alt_m,
    )
    return [point, pad]


def order_source_pts_for_track(source_pts: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Order mesh points so consecutive ``ONE TRACK`` hops stay spatially local.

    Lattice meshes sorted by ``(x, y)`` walk a column then jump from the last
    row of column *i* to the first row of column *i+1* — a domain-width hop
    that often clips terrain even when both endpoints are above ground.
    Snaking *y* each column keeps that wrap to one cell.
    """
    if len(source_pts) <= 1:
        return source_pts
    ordered = source_pts.copy()
    ordered["_sort_x"] = ordered.geometry.x
    ordered["_sort_y"] = ordered.geometry.y
    columns: list[pd.DataFrame] = []
    for col_i, (_, column) in enumerate(ordered.groupby("_sort_x", sort=True)):
        columns.append(
            column.sort_values("_sort_y", ascending=(col_i % 2 == 0)),
        )
    return gpd.GeoDataFrame(
        pd.concat(columns),
        crs=source_pts.crs,
    ).drop(columns=["_sort_x", "_sort_y"])
