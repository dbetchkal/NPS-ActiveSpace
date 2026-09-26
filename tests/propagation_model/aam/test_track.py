"""Tests for AAM mesh ordering and single-point track padding."""

from __future__ import annotations

import math

import geopandas as gpd
import pytest
from shapely.geometry import Point

pytest.importorskip("aam_translator")

from aam_translator.write_inp import TrackPoint

from nps_active_space.propagation_model.aam.track import (
    SINGLE_TRACK_PAD_M,
    order_source_pts_for_track,
    pad_single_point_track,
)


class TestOrderSourcePtsForTrack:
    def test_order_source_pts_sorts_by_xy(self) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2]},
            geometry=[Point(2, 0, 100), Point(0, 0, 100), Point(1, 0, 100)],
            crs="EPSG:32606",
        )
        ordered = order_source_pts_for_track(pts)
        assert ordered["id"].tolist() == [1, 2, 0]

    def test_order_source_pts_snakes_grid_columns(self) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2, 3]},
            geometry=[
                Point(0, 0, 100),
                Point(0, 1, 100),
                Point(1, 0, 100),
                Point(1, 1, 100),
            ],
            crs="EPSG:32606",
        )
        ordered = order_source_pts_for_track(pts)
        assert ordered["id"].tolist() == [0, 1, 3, 2]


class TestPadSinglePointTrack:
    def test_leaves_multi_point_track_unchanged(self) -> None:
        track = [
            TrackPoint(lon=-148.87, lat=63.66, alt_m=1500.0),
            TrackPoint(lon=-148.86, lat=63.66, alt_m=1500.0),
        ]
        assert pad_single_point_track(track) == track

    def test_pads_one_vertex_about_one_meter_east(self) -> None:
        point = TrackPoint(lon=-148.87, lat=63.66, alt_m=1500.0)
        padded = pad_single_point_track([point])
        assert len(padded) == 2
        assert padded[0] == point
        assert padded[1].lat == point.lat
        assert padded[1].alt_m == point.alt_m
        meters_per_deg_lon = 111_320.0 * abs(math.cos(math.radians(point.lat)))
        east_m = (padded[1].lon - padded[0].lon) * meters_per_deg_lon
        assert east_m == pytest.approx(SINGLE_TRACK_PAD_M, rel=1e-4)
