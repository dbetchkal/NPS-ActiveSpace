import argparse
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import LineString, Point

from nps_active_space.scripts.viz import (
    create_polyline_3d,
    densify_linestring,
    parse_deployment,
    parse_existing_file,
    parse_iso_date,
    parse_max_tracks,
    sea_surface_z_profile,
    track_points_to_linestring,
)


class TestParseDeployment:
    def test_valid_deployment(self):
        assert parse_deployment("DENATRLA2024") == ("DENA", "TRLA", 2024)

    def test_rejects_too_short(self):
        with pytest.raises(argparse.ArgumentTypeError, match="at least 9 characters"):
            parse_deployment("DENA2024")

    def test_rejects_unit_year_only(self):
        with pytest.raises(argparse.ArgumentTypeError, match="at least 9 characters"):
            parse_deployment("DENA2024")

    def test_rejects_non_digit_year(self):
        with pytest.raises(argparse.ArgumentTypeError, match="4 digits"):
            parse_deployment("DENATRLA20XX")


class TestParseIsoDate:
    def test_valid_date(self):
        assert parse_iso_date("2024-06-18", arg_name="--start-date") == "2024-06-18"

    def test_rejects_bad_format(self):
        with pytest.raises(argparse.ArgumentTypeError, match="YYYY-MM-DD"):
            parse_iso_date("06/18/2024", arg_name="--start-date")

    def test_includes_arg_name_in_error(self):
        with pytest.raises(argparse.ArgumentTypeError, match="--end-date"):
            parse_iso_date("not-a-date", arg_name="--end-date")


class TestParseExistingFile:
    def test_existing_file(self, tmp_path: Path):
        file_path = tmp_path / "annotations.geojson"
        file_path.write_text("{}")
        assert parse_existing_file(str(file_path), arg_name="--annotation-file") == str(
            file_path
        )

    def test_missing_file(self, tmp_path: Path):
        with pytest.raises(argparse.ArgumentTypeError, match="file not found"):
            parse_existing_file(str(tmp_path / "missing.geojson"), arg_name="--transits-pkl")


class TestParseMaxTracks:
    def test_valid_positive_int(self):
        assert parse_max_tracks("500") == 500

    def test_rejects_zero(self):
        with pytest.raises(argparse.ArgumentTypeError, match="positive integer"):
            parse_max_tracks("0")

    def test_rejects_non_integer(self):
        with pytest.raises(argparse.ArgumentTypeError, match="positive integer"):
            parse_max_tracks("abc")


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
