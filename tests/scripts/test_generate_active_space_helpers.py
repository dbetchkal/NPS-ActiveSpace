"""Unit tests for small helpers extracted from generate_active_space.py."""

import geopandas as gpd

import nps_active_space.scripts.generate_active_space as generate_active_space
from shapely.geometry import Polygon


def _empty_active_space() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(geometry=[Polygon()], crs="EPSG:4326")


def _nonempty_active_space() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs="EPSG:4326",
    )


class TestNonemptyActiveSpaceCount:
    def test_counts_only_layers_with_geometry(self):
        results = [
            ("O_+010", _nonempty_active_space()),
            ("O_+020", gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")),
            ("O_+030", _empty_active_space()),
        ]

        assert generate_active_space._nonempty_active_space_count(results) == 1
