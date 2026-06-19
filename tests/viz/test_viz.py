import argparse
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import LineString, MultiLineString, Point

from nps_active_space.viz import (
    Visualizer,
    annotation_z_profile,
    create_polyline_3d,
    densify_linestring,
    format_annotation_summary,
    is_surface_track,
    iter_plot_linestrings,
    parse_deployment,
    parse_existing_file,
    parse_iso_date,
    parse_max_tracks,
    resolve_viz_plot_flags,
    sea_surface_z_profile,
    track_points_to_linestring,
    utm_orientation_axes_kwargs,
    vertex_z_from_coord,
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


class TestResolveVizPlotFlags:
    def test_annotation_file_implies_annotations(self):
        flags = resolve_viz_plot_flags(annotation_file="/tmp/a.geojson")
        assert flags == (False, True, False, False)

    def test_transits_pkl_implies_transits(self):
        flags = resolve_viz_plot_flags(transits_pkl="/tmp/t.pkl")
        assert flags == (False, False, True, False)

    def test_all_enables_standard_layers(self):
        flags = resolve_viz_plot_flags(plot_all=True)
        assert flags == (True, True, True, False)


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
            "nps_active_space.viz.elevation.get_elevation",
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
            "nps_active_space.viz.elevation.get_elevation",
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


class TestFormatAnnotationSummary:
    def test_empty(self):
        import geopandas as gpd

        assert format_annotation_summary(gpd.GeoDataFrame()) == "0 segments"

    def test_includes_track_and_surface_counts(self):
        import geopandas as gpd

        gdf = gpd.GeoDataFrame(
            {
                "_id": ["T1", "T1"],
                "audible": [True, False],
            },
            geometry=[
                LineString([(-136.0, 58.0), (-136.1, 58.1)]),
                LineString([(0.0, 0.0, 100.0), (1.0, 1.0, 200.0)]),
            ],
            crs="EPSG:4326",
        )
        summary = format_annotation_summary(gdf)
        assert "2 segments" in summary
        assert "1 tracks" in summary
        assert "1 audible" in summary
        assert "1 sea-surface" in summary
        assert "1 elevated" in summary


class TestAnnotationZProfile:
    def test_airborne_uses_stored_z_without_dem(self):
        line = LineString([(0, 0, 1500.0), (1, 1, 2000.0), (2, 2, 1800.0)])

        class _FailSampler:
            def sample_utm_many(self, x, y):
                raise AssertionError("DEM sampling should be skipped for airborne tracks")

        z = annotation_z_profile(line, _FailSampler())
        assert z.tolist() == [1500.0, 2000.0, 1800.0]

    def test_surface_vertices_clamped_to_dem(self):
        line = LineString([(0, 0, 0.0), (1, 1, 250.0)])

        class _FakeSampler:
            def sample_utm_many(self, x, y):
                return np.full(len(x), 100.0)

        z = annotation_z_profile(line, _FakeSampler(), offset_m=2.0)
        assert z[0] == pytest.approx(102.0)
        assert z[1] == pytest.approx(250.0)


class TestOrientationWidgets:
    def test_window_title_constant(self):
        from nps_active_space.viz.markers import WINDOW_TITLE

        assert WINDOW_TITLE == "NPS ActiveSpace Visualization"

    def test_bundled_waypoint_icon_exists(self):
        from nps_active_space.viz.markers import default_window_icon_path

        assert default_window_icon_path().is_file()

    def test_utm_axes_use_east_north_up_labels(self):
        kwargs = utm_orientation_axes_kwargs()
        assert kwargs["xlabel"] == "E"
        assert kwargs["ylabel"] == "N"
        assert kwargs["zlabel"] == "Z"
        assert kwargs["viewport"] == (0, 0, 0.2, 0.2)


class TestFlatSeaSurfacePolyline:
    def test_uses_constant_offset_without_dem(self):
        line = LineString([(0, 0), (1000, 0), (2000, 0)])
        poly = Visualizer._flat_sea_surface_polyline(line, offset_m=5.0)
        assert poly.n_points == 3
        assert np.all(poly.points[:, 2] == pytest.approx(5.0))


class TestVertexZAndSurfaceTrack:
    def test_missing_z_is_sea_level(self):
        assert vertex_z_from_coord((-136.0, 58.0)) == 0.0

    def test_nan_z_is_sea_level(self):
        assert vertex_z_from_coord((-136.0, 58.0, float("nan"))) == 0.0

    def test_nan_z_counts_as_surface_track(self):
        line = LineString([(-136.0, 58.0, float("nan")), (-136.1, 58.1, float("nan"))])
        assert is_surface_track(np.array(line.coords))

    def test_positive_z_is_not_surface_track(self):
        line = LineString([(0, 0, 100.0), (1, 1, 200.0)])
        assert not is_surface_track(np.array(line.coords))


class TestIterPlotLinestrings:
    def test_point_becomes_short_line(self):
        lines = iter_plot_linestrings(Point(1.0, 2.0))
        assert len(lines) == 1
        assert len(lines[0].coords) == 2

    def test_multilinestring_splits(self):
        geom = MultiLineString([LineString([(0, 0), (1, 1)]), LineString([(2, 2), (3, 3)])])
        assert len(iter_plot_linestrings(geom)) == 2


class TestCreatePolyline3d:
    def test_nan_z_becomes_zero(self):
        line = LineString([(0, 0, float("nan")), (1, 0, float("nan"))])
        poly = create_polyline_3d(line, z=np.array([float("nan"), 5.0]))
        assert poly.points[0, 2] == pytest.approx(0.0)
        assert poly.points[1, 2] == pytest.approx(5.0)

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


class TestScriptsVizShim:
    def test_scripts_viz_delegates_to_viz_cli(self):
        from nps_active_space.scripts import viz as scripts_viz
        from nps_active_space.viz.cli import main

        assert scripts_viz.main is main
