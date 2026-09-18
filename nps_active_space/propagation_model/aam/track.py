"""Order mesh points and pad tracks for AAM ``ONE TRACK`` batches."""

from __future__ import annotations

import math

import geopandas as gpd
import pandas as pd
from aam_translator.write_inp import TrackPoint

# AAM 3.0.0 crashes on a 1-vertex ONE TRACK (Wine exit 152; related Fortran crash
# whose stderr often mentions the internal array FPA; empty .POI). Pad ~1 m so a
# leftover singleton stays two vertices. See aam-translator
# docs/reading_aam_output.md and references/notes/aam_inp_format.md.
SINGLE_TRACK_PAD_M = 1.0
METERS_PER_DEG_LAT = 111_320.0

__all__ = [
    "METERS_PER_DEG_LAT",
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
    meters_per_deg_lon = METERS_PER_DEG_LAT * max(abs(cos_lat), 1e-6)
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
