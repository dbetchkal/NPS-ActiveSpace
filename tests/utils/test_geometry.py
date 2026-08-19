import pytest
from shapely.geometry import LineString, Point

from nps_active_space.utils.geometry import linestring_from_coords, linestring_from_point


class TestLinestringFromCoords:
    def test_empty_coords_returns_empty_linestring(self):
        line = linestring_from_coords([])
        assert line.is_empty

    def test_two_points_returns_linestring(self):
        line = linestring_from_coords([(0.0, 0.0), (1.0, 0.0)])
        assert len(line.coords) == 2

    def test_single_point_defaults_to_short_segment(self):
        line = linestring_from_coords([(1.0, 2.0)])
        assert len(line.coords) == 2
        assert line.coords[0] == pytest.approx((1.0, 2.0))
        assert line.coords[1][0] == pytest.approx(11.0)
        assert line.coords[1][1] == pytest.approx(2.0)

    def test_single_point_preserves_z(self):
        line = linestring_from_coords([(1.0, 2.0, 100.0)])
        assert line.coords[0] == pytest.approx((1.0, 2.0, 100.0))
        assert line.coords[1] == pytest.approx((11.0, 2.0, 100.0))

    def test_single_point_empty_mode(self):
        line = linestring_from_coords([(1.0, 2.0)], on_single_point="empty")
        assert line.is_empty

    def test_single_point_does_not_raise_geos_error(self):
        line = linestring_from_coords([(500_000.0, 6_000_000.0)])
        assert not line.is_empty


class TestLinestringFromPoint:
    def test_point_becomes_short_line(self):
        line = linestring_from_point(Point(1.0, 2.0))
        assert len(line.coords) == 2
