"""Tests for active-space tricontour geometry helpers."""

import matplotlib

matplotlib.use("Agg")

import geopandas as gpd
import numpy as np
from shapely.geometry import Point

from nps_active_space.active_space.active_space_geometry import (
    build_active_space,
    contour_active_space,
)

_UTM_CRS = "EPSG:32606"
_AUDIBLE_RADIUS_M = 30.0
_GRID_HALF_EXTENT_M = 80.0
_GRID_SPACING_M = 10.0
_CENTER_X = 500_000.0
_CENTER_Y = 6_000_000.0


def _audibility_point_cloud(crs: str) -> gpd.GeoDataFrame:
    """Audible disk inside an inaudible square grid."""
    coords = np.arange(
        -_GRID_HALF_EXTENT_M,
        _GRID_HALF_EXTENT_M + _GRID_SPACING_M,
        _GRID_SPACING_M,
    )
    xs, ys = np.meshgrid(coords, coords)
    xs = xs.ravel() + _CENTER_X
    ys = ys.ravel() + _CENTER_Y
    audible = (xs - _CENTER_X) ** 2 + (ys - _CENTER_Y) ** 2 < _AUDIBLE_RADIUS_M**2
    return gpd.GeoDataFrame(
        {"audible": audible.astype(int)},
        geometry=gpd.points_from_xy(xs, ys),
        crs=crs,
    )


class TestContourActiveSpace:
    def test_contour_points_at_altitude_with_preserved_crs(self) -> None:
        total_space = _audibility_point_cloud(_UTM_CRS)
        altitude_m = 3_658

        contour_pts = contour_active_space(total_space, altitude_m)

        assert not contour_pts.empty
        assert contour_pts.crs == total_space.crs
        assert np.allclose(contour_pts.geometry.z, altitude_m)

        distances = np.hypot(
            contour_pts.geometry.x - _CENTER_X,
            contour_pts.geometry.y - _CENTER_Y,
        )
        inner_extent = _AUDIBLE_RADIUS_M * 0.6
        outer_extent = _AUDIBLE_RADIUS_M * 1.8
        assert distances.min() > inner_extent
        assert distances.max() < outer_extent


class TestBuildActiveSpace:
    def test_build_active_space_polygon_contains_audible_center(self) -> None:
        total_space = _audibility_point_cloud(_UTM_CRS)

        active_space = build_active_space(total_space, _UTM_CRS)

        assert not active_space.empty
        assert active_space.crs.to_string() == _UTM_CRS

        union = active_space.union_all()
        center = Point(_CENTER_X, _CENTER_Y)
        far_corner = Point(
            _CENTER_X + _GRID_HALF_EXTENT_M,
            _CENTER_Y + _GRID_HALF_EXTENT_M,
        )

        assert union.contains(center)
        assert not union.contains(far_corner)
