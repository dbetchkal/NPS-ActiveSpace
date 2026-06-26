from __future__ import annotations

import geopandas as gpd
import numpy as np
import pyvista as pv
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon


def polygon_to_mesh(polygon: Polygon, z: float) -> pv.PolyData:
    exterior = np.array(polygon.exterior.coords)
    points = np.c_[exterior[:, 0], exterior[:, 1], np.full(len(exterior), z)]
    poly = pv.PolyData(points)
    poly.faces = np.hstack([[len(points)], np.arange(len(points))])
    return poly.triangulate()


def active_to_polys(active: gpd.GeoDataFrame) -> list[Polygon]:
    geometry = active.geometry.iloc[0]
    if isinstance(geometry, Polygon):
        polys = [geometry]
    elif isinstance(geometry, MultiPolygon):
        polys = geometry.geoms
    return polys


def active_to_linestrings(active: gpd.GeoDataFrame) -> list[LineString]:
    geometry = active.boundary.iloc[0]
    if isinstance(geometry, MultiLineString):
        linestrings = geometry.geoms
    elif isinstance(geometry, LineString):
        linestrings = [geometry]
    return linestrings


def _point_as_line(point: Point, span_m: float = 10.0) -> LineString:
    """Convert a single-point geometry into a short renderable segment."""
    x, y = point.x, point.y
    return LineString([(x, y), (x + span_m, y)])


def iter_plot_linestrings(geometry) -> list[LineString]:
    """Normalize annotation geometries to LineStrings PyVista can draw."""
    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, Point):
        return [_point_as_line(geometry)]
    if isinstance(geometry, LineString):
        if len(geometry.coords) < 2:
            return [_point_as_line(Point(geometry.coords[0]))]
        return [geometry]
    if isinstance(geometry, MultiLineString):
        lines: list[LineString] = []
        for part in geometry.geoms:
            lines.extend(iter_plot_linestrings(part))
        return lines
    return []


def densify_linestring(linestring: LineString, step_m: float) -> LineString:
    """Insert vertices along a line so elevation can be sampled between sparse AIS points."""
    if linestring.is_empty or linestring.length <= step_m:
        return linestring
    distances = np.arange(0, linestring.length, step_m)
    if distances[-1] < linestring.length:
        distances = np.append(distances, linestring.length)
    return LineString([linestring.interpolate(float(d)) for d in distances])


def create_polyline_3d(
    linestring: LineString,
    z: float | np.ndarray | None = None,
) -> pv.PolyData:
    coords = np.array(linestring.coords)
    xy = coords[:, :2]
    if z is not None:
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        if z_arr.size != coords.shape[0]:
            raise ValueError(
                f"z must be scalar or have one value per vertex ({coords.shape[0]}), got {z_arr.size}"
            )
        coords = np.column_stack((xy, z_arr))
    elif coords.shape[1] == 2:
        coords = np.column_stack((xy, np.zeros(coords.shape[0])))
    assert coords.shape[1] == 3
    coords[:, 2] = np.clip(np.nan_to_num(coords[:, 2], nan=0.0), 0, 10000)
    n_points = coords.shape[0]
    lines = np.hstack([[n_points], np.arange(n_points)])
    poly = pv.PolyData(coords)
    poly.lines = lines
    return poly


def flat_sea_surface_polyline(linestring: LineString, offset_m: float) -> pv.PolyData:
    """Sea-level track at a constant offset — no DEM resampling."""
    z = np.full(len(linestring.coords), offset_m)
    return create_polyline_3d(linestring, z=z)


def track_points_to_linestring(
    track_points: gpd.GeoDataFrame,
    *,
    include_z: bool = False,
) -> LineString:
    """Connect ordered track points into a line in the plot CRS."""
    if include_z and "z" in track_points.columns:
        coords = [
            (geom.x, geom.y, float(z) if np.isfinite(z) else 0.0)
            for geom, z in zip(track_points.geometry, track_points["z"], strict=True)
        ]
    else:
        coords = [(geom.x, geom.y) for geom in track_points.geometry]
    if len(coords) < 2:
        return LineString(coords)
    return LineString(coords)
