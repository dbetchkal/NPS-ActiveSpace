import numpy as np
import pytest
from shapely.geometry import LineString, Point

from nps_active_space.scripts.viz import (
    create_polyline_3d,
    densify_linestring,
    sea_surface_z_profile,
    track_points_to_linestring,
)


class TestDensifyLinestring:
    def test_short_line_unchanged(self):
        line = LineString([(0, 0), (50, 0)])
        assert densify_linestring(line, step_m=100) == line

    def test_long_line_gets_intermediate_vertices(self):
        line = LineString([(0, 0), (250, 0)])
        dense = densify_linestring(line, step_m=100)
        assert len(dense.coords) == 4
        assert dense.coords[1] == pytest.approx((100.0, 0.0))


class TestSeaSurfaceZProfile:
    def test_water_pixels_use_offset_above_zero(self, monkeypatch):
        line = LineString([(0, 0), (250, 0)])

        def fake_elevation(_dem, _lon, _lat):
            return 0.0

        monkeypatch.setattr(
            "nps_active_space.scripts.viz.get_elevation",
            fake_elevation,
        )
        dense, z_vals = sea_surface_z_profile(
            line,
            dem=object(),
            crs="epsg:32608",
            offset_m=5.0,
            densify_step_m=100.0,
        )
        assert len(dense.coords) == 4
        assert z_vals.tolist() == [5.0, 5.0, 5.0, 5.0]

    def test_negative_dem_clamped_to_sea_level(self, monkeypatch):
        line = LineString([(0, 0), (1, 0)])

        monkeypatch.setattr(
            "nps_active_space.scripts.viz.get_elevation",
            lambda _dem, _lon, _lat: -12.0,
        )
        _, z_vals = sea_surface_z_profile(
            line,
            dem=object(),
            crs="epsg:32608",
            offset_m=2.0,
            densify_step_m=100.0,
        )
        assert z_vals.tolist() == [2.0, 2.0]


class TestCreatePolyline3d:
    def test_rejects_mismatched_z_length(self):
        line = LineString([(0, 0), (1, 1), (2, 2)])
        with pytest.raises(ValueError, match="one value per vertex"):
            create_polyline_3d(line, z=np.array([1.0, 2.0]))


class TestTrackPointsToLinestring:
    def test_builds_line_from_points(self):
        import geopandas as gpd

        gdf = gpd.GeoDataFrame(
            geometry=[Point(0, 0), Point(1, 0), Point(2, 0)],
            crs="epsg:32608",
        )
        line = track_points_to_linestring(gdf)
        assert len(line.coords) == 3
