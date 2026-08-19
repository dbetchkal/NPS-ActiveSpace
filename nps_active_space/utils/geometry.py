from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

from shapely.geometry import LineString, Point

Coordinate = tuple[float, ...]
SinglePointMode = Literal["segment", "empty"]


def linestring_from_coords(
    coords: Sequence[Coordinate],
    *,
    on_single_point: SinglePointMode = "segment",
    segment_span_m: float = 10.0,
) -> LineString:
    """Build a valid Shapely LineString from coordinate sequences.

    Shapely rejects LineStrings with exactly one vertex. This helper returns an
    empty LineString for no coordinates and normalizes single-point inputs
    according to ``on_single_point``.

    ``segment_span_m`` is in the same planar units as ``coords`` (meters in UTM,
    degrees in WGS84), not necessarily meters.
    """
    if not coords:
        return LineString()
    if len(coords) == 1:
        if on_single_point == "empty":
            return LineString()
        return LineString(_single_point_segment(coords[0], segment_span_m))
    return LineString(coords)


def linestring_from_point(
    point: Point,
    *,
    segment_span_m: float = 10.0,
) -> LineString:
    """Convert a Point into a short renderable LineString segment."""
    return linestring_from_coords(
        [point.coords[0]],
        segment_span_m=segment_span_m,
    )


def _single_point_segment(coord: Coordinate, span_m: float) -> list[Coordinate]:
    x, y = coord[0], coord[1]
    if len(coord) >= 3:
        z = coord[2]
        return [(x, y, z), (x + span_m, y, z)]
    return [(x, y), (x + span_m, y)]
