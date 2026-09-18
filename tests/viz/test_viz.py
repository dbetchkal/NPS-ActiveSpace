import argparse
from configparser import NoOptionError
from pathlib import Path
from unittest.mock import MagicMock

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, MultiLineString, Point, box

from nps_active_space.ground_truthing.load_tracks import LoadedTracks
from nps_active_space.utils.enums import AcousticModel, TrackSource
from nps_active_space.utils.models import Tracks
from nps_active_space.viz.active_space_plot import ActiveSpacePlotter
from nps_active_space.viz.polyline_builders import TrackPolylineBuilder
from nps_active_space.viz.scene import ScenePlotter
from nps_active_space.viz.session import VizSession
from nps_active_space.viz.style import VizStyle
from nps_active_space.viz.tracks_plot import TracksPlotter
from nps_active_space.viz.widgets import PlotWidgets
from nps_active_space.viz import (
    Visualizer,
    annotation_z_profile,
    create_polyline_3d,
    densify_linestring,
    format_annotation_summary,
    is_surface_track,
    iter_plot_linestrings,
    parse_existing_file,
    parse_iso_date,
    parse_max_tracks,
    resolve_track_source_args,
    resolve_viz_plot_flags,
    sea_surface_z_profile,
    track_points_to_linestring,
    utm_orientation_axes_kwargs,
    vertex_z_from_coord,
)


class TestVizCliArgs:
    def test_parses_workflow_flags(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import sys

        from nps_active_space.viz.cli import main

        captured: list[tuple] = []

        def _capture_visualizer(*args, **kwargs):
            captured.append((args, kwargs))

        monkeypatch.setattr("nps_active_space.viz.cli.Visualizer", _capture_visualizer)
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "viz",
                "-e",
                "DENA_example",
                "-u",
                "DENA",
                "-s",
                "TRLA",
                "-y",
                "2024",
                "-A",
                "-a",
            ],
        )
        main()
        unit, site, year, environment, do_active = captured[0][0][:5]
        assert (unit, site, year, environment, do_active) == (
            "DENA",
            "TRLA",
            2024,
            "DENA_example",
            True,
        )

    def test_model_aam_implies_active_and_forwards_model(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import sys

        from nps_active_space.viz.cli import main

        captured: list[tuple] = []

        def _capture_visualizer(*args, **kwargs):
            captured.append((args, kwargs))

        monkeypatch.setattr("nps_active_space.viz.cli.Visualizer", _capture_visualizer)
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "viz",
                "-e",
                "DENA_example",
                "-u",
                "DENA",
                "-s",
                "TRLA",
                "-y",
                "2024",
                "--model",
                "aam",
            ],
        )
        main()
        assert captured[0][0][4] is True
        assert captured[0][1]["model"] is AcousticModel.AAM
        assert captured[0][1]["compare_models"] is False

    def test_compare_implies_active(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import sys

        from nps_active_space.viz.cli import main

        captured: list[tuple] = []

        def _capture_visualizer(*args, **kwargs):
            captured.append((args, kwargs))

        monkeypatch.setattr("nps_active_space.viz.cli.Visualizer", _capture_visualizer)
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "viz",
                "-e",
                "DENA_example",
                "-u",
                "DENA",
                "-s",
                "TRLA",
                "-y",
                "2024",
                "--compare",
            ],
        )
        main()
        assert captured[0][0][4] is True
        assert captured[0][1]["compare_models"] is True

    def test_requires_environment(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import sys

        from nps_active_space.viz.cli import main

        monkeypatch.setattr(
            sys,
            "argv",
            ["viz", "-u", "DENA", "-s", "TRLA", "-y", "2024"],
        )
        with pytest.raises(SystemExit):
            main()

    def test_requires_unit_site_year(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import sys

        from nps_active_space.viz.cli import main

        monkeypatch.setattr(sys, "argv", ["viz", "-e", "DENA_example"])
        with pytest.raises(SystemExit):
            main()


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


class TestResolveVizAnnotationFile:
    def test_resolves_basename_under_site_dir(self, tmp_path: Path) -> None:
        from nps_active_space.viz.cli import resolve_viz_annotation_file

        project_dir = tmp_path / "proj"
        site = project_dir / "DENATRLA"
        site.mkdir(parents=True)
        annot = site / "DENATRLA2025_saved_annotations.geojson"
        annot.write_text("{}")

        resolved = resolve_viz_annotation_file(
            str(project_dir), "DENA", "TRLA", "DENATRLA2025_saved_annotations.geojson"
        )
        assert Path(resolved) == annot.resolve()

    def test_absolute_path(self, tmp_path: Path) -> None:
        from nps_active_space.viz.cli import resolve_viz_annotation_file

        annot = tmp_path / "custom.geojson"
        annot.write_text("{}")
        resolved = resolve_viz_annotation_file(
            str(tmp_path / "proj"), "DENA", "TRLA", str(annot)
        )
        assert Path(resolved) == annot.resolve()

    def test_missing_raises(self, tmp_path: Path) -> None:
        from nps_active_space.viz.cli import resolve_viz_annotation_file

        with pytest.raises(argparse.ArgumentTypeError, match="also looked in"):
            resolve_viz_annotation_file(str(tmp_path), "DENA", "TRLA", "nope.geojson")


class TestResolveVizPlotFlags:
    def test_annotation_file_implies_annotations(self):
        flags = resolve_viz_plot_flags(annotation_file="/tmp/a.geojson")
        assert flags == (False, True, False, None)

    def test_transits_pkl_implies_transits(self):
        flags = resolve_viz_plot_flags(transits_pkl="/tmp/t.pkl")
        assert flags == (False, False, True, None)

    def test_all_enables_standard_layers(self):
        flags = resolve_viz_plot_flags(plot_all=True)
        assert flags == (True, True, True, None)

    def test_compare_implies_active_space(self):
        flags = resolve_viz_plot_flags(compare=True)
        assert flags[0] is True
        assert flags[1:] == (False, False, None)

    def test_model_implies_active_space(self):
        flags = resolve_viz_plot_flags(model=AcousticModel.AAM)
        assert flags[0] is True

    @pytest.mark.parametrize("source", list(TrackSource))
    def test_track_source_passed_through(self, source: TrackSource) -> None:
        flags = resolve_viz_plot_flags(track_source=source)
        assert flags == (False, False, False, source)


class TestResolveTrackSourceArgs:
    def test_dates_require_track_source(self) -> None:
        parser = argparse.ArgumentParser()
        args = argparse.Namespace(
            track_source=None,
            start_date="2024-05-24",
            end_date=None,
        )
        with pytest.raises(SystemExit):
            resolve_track_source_args(args, parser)

    def test_rejects_inverted_date_range(self) -> None:
        parser = argparse.ArgumentParser()
        args = argparse.Namespace(
            track_source=TrackSource.AIS,
            start_date="2024-05-25",
            end_date="2024-05-24",
        )
        with pytest.raises(SystemExit):
            resolve_track_source_args(args, parser)

    def test_passes_through_track_source_with_dates(self) -> None:
        parser = argparse.ArgumentParser()
        args = argparse.Namespace(
            track_source=TrackSource.ADSB,
            start_date="2024-05-24",
            end_date="2024-05-25",
        )
        assert resolve_track_source_args(args, parser) is TrackSource.ADSB


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


class TestDemElevationSampler:
    def test_sample_utm_many_vectorized(self) -> None:
        from rasterio.transform import from_bounds

        from nps_active_space.viz.elevation import DemElevationSampler

        band = np.arange(100, dtype=float).reshape(10, 10)
        dem_transform = from_bounds(0, 0, 100, 100, 10, 10)

        class _Dem:
            crs = "epsg:32608"
            nodata = None
            transform = dem_transform

        dem = _Dem()
        sampler = DemElevationSampler(dem, band, plot_crs="epsg:32608")

        xs = np.linspace(10, 90, 5)
        ys = np.linspace(10, 90, 5)
        elev = sampler.sample_utm_many(xs, ys)

        assert elev.shape == (5,)
        assert np.all(np.isfinite(elev))


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

    def test_utm_axes_use_east_north_up_labels(self):
        kwargs = utm_orientation_axes_kwargs()
        assert kwargs["xlabel"] == "E"
        assert kwargs["ylabel"] == "N"
        assert kwargs["zlabel"] == "Z"
        assert kwargs["viewport"] == (0.76, 0.03, 0.96, 0.21)


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

    def test_valid_linestring_is_returned_as_is(self):
        line = LineString([(0, 0), (1, 1), (2, 2)])
        lines = iter_plot_linestrings(line)
        assert lines == [line]


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

    @pytest.mark.parametrize("altitude", [2100, np.int64(2100), np.array(2100.0)])
    def test_scalar_layer_altitude_broadcasts(self, altitude):
        line = LineString([(0, 0), (1, 0), (2, 0)])
        poly = create_polyline_3d(line, z=altitude)
        assert poly.n_points == 3
        assert np.all(poly.points[:, 2] == pytest.approx(2100.0))


class TestTrackPointsToLinestring:
    def test_builds_line_from_points(self):
        import geopandas as gpd

        gdf = gpd.GeoDataFrame(
            geometry=[Point(0, 0), Point(1, 0), Point(2, 0)],
            crs="epsg:32608",
        )
        line = track_points_to_linestring(gdf)
        assert len(line.coords) == 3

    def test_single_point_track_returns_renderable_line(self):
        import geopandas as gpd

        gdf = gpd.GeoDataFrame(
            geometry=[Point(500_000.0, 6_000_000.0)],
            crs="epsg:32608",
        )
        line = track_points_to_linestring(gdf)
        assert len(line.coords) == 2

    def test_include_z_adds_altitude_to_coords(self):
        import geopandas as gpd

        gdf = gpd.GeoDataFrame(
            {"z": [1000.0, 1500.0, 2000.0]},
            geometry=[Point(0, 0), Point(1, 0), Point(2, 0)],
            crs="epsg:32608",
        )
        line = track_points_to_linestring(gdf, include_z=True)
        assert line.coords[1] == pytest.approx((1.0, 0.0, 1500.0))


def _attach_viz_plotters(vis: Visualizer) -> Visualizer:
    """Wire composed plotters for Visualizer instances built with __new__."""
    style = VizStyle()
    session = VizSession(
        unit=vis.unit,
        site=vis.site,
        year=vis.year,
        deployment=f"{vis.unit}{vis.site}{vis.year}",
        project_dir=vis.project_dir,
        crs=vis.crs,
        study_area=getattr(
            vis,
            "study_area",
            gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)], crs=vis.crs),
        ),
        plotter=vis.plotter,
        logger=vis.logger,
        style=style,
        fill_layers=getattr(vis, "fill_layers", False),
        max_tracks=getattr(vis, "max_tracks", 1000),
        to_wgs84=MagicMock(),
        dem_sampler=getattr(vis, "_dem_sampler", None),
        dem=getattr(vis, "dem", None),
    )
    vis._session = session
    vis._widgets = PlotWidgets(session)
    vis._scene = ScenePlotter(session)
    vis._polylines = TrackPolylineBuilder(session)
    vis._tracks = TracksPlotter(session, vis._widgets, vis._scene, vis._polylines)
    vis._active_space = ActiveSpacePlotter(session, vis._widgets)
    vis.study_area = session.study_area
    return vis


class TestPlotTracksSourceRouting:
    @staticmethod
    def _stub_visualizer() -> Visualizer:
        vis = Visualizer.__new__(Visualizer)
        vis.unit, vis.site, vis.year = "GLBA", "LSTL", 2024
        vis.project_dir = "/proj"
        vis.max_tracks = 500
        vis.crs = "epsg:32608"
        vis.study_area = gpd.GeoDataFrame(
            geometry=[box(500_000, 6_000_000, 510_000, 6_010_000)],
            crs=vis.crs,
        )
        vis.plotter = MagicMock()
        vis.logger = MagicMock()
        vis.sea_surface_offset_m = 5.0
        vis.sea_surface_densify_step_m = 100.0
        vis.vessel_track_color = "magenta"
        vis.flight_track_color = "yellow"
        vis = _attach_viz_plotters(vis)
        vis.dem = MagicMock()
        vis._dem_sampler = MagicMock()
        return vis

    @staticmethod
    def _loaded_tracks(*, flight: bool) -> LoadedTracks:
        id_col = "flight_id" if flight else "event_id"
        z = 5000.0 if flight else 0.0
        gdf = gpd.GeoDataFrame(
            {
                id_col: ["t1", "t1"],
                "TIME": [
                    pd.Timestamp("2024-05-24 10:00"),
                    pd.Timestamp("2024-05-24 10:05"),
                ],
                "altitude": [z, z],
            },
            geometry=gpd.points_from_xy([500_100.0, 500_200.0], [6_000_100.0, 6_000_200.0]),
            crs="epsg:32608",
        )
        tracks = Tracks(gdf, id_col=id_col, datetime_col="TIME", z_col="altitude")
        return LoadedTracks(tracks, None, None)

    @pytest.mark.parametrize(
        ("source", "expect_sea", "expect_flight", "use_flight_tracks"),
        [
            (TrackSource.AIS, True, False, False),
            (TrackSource.ADSB, False, True, True),
            (TrackSource.GPS, False, True, True),
        ],
    )
    def test_routes_polyline_builder_by_source(
        self,
        monkeypatch: pytest.MonkeyPatch,
        source: TrackSource,
        expect_sea: bool,
        expect_flight: bool,
        use_flight_tracks: bool,
    ) -> None:
        vis = self._stub_visualizer()
        captured_lines: list[LineString] = []
        sea = MagicMock(side_effect=lambda line: captured_lines.append(line) or MagicMock())
        flight = MagicMock(side_effect=lambda line: captured_lines.append(line) or MagicMock())
        monkeypatch.setattr(vis, "_sea_surface_polyline", sea)
        monkeypatch.setattr(vis, "_annotation_polyline", flight)
        monkeypatch.setattr(vis, "_add_track_line", lambda polyline, *, color: MagicMock())
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.get_deployment",
            lambda *args, **kwargs: MagicMock(lat=58.4, lon=-136.0),
        )
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.load_tracks",
            lambda *args, **kwargs: self._loaded_tracks(flight=use_flight_tracks),
        )

        vis.plot_tracks(source, "2024-05-24", "2024-05-24")

        assert sea.called is expect_sea
        assert flight.called is expect_flight
        coord_dim = len(captured_lines[0].coords[0])
        if use_flight_tracks:
            assert coord_dim == 3
        else:
            assert coord_dim == 2

    def test_forwards_load_tracks_kwargs(self, monkeypatch: pytest.MonkeyPatch) -> None:
        vis = self._stub_visualizer()
        calls: list[dict] = []
        monkeypatch.setattr(vis, "_add_track_line", lambda polyline, *, color: MagicMock())
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.get_deployment",
            lambda *args, **kwargs: MagicMock(lat=58.4, lon=-136.0),
        )
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.load_tracks",
            lambda *args, **kwargs: calls.append(kwargs) or self._loaded_tracks(flight=True),
        )

        vis.plot_tracks(TrackSource.ADSB, "2024-05-24", "2024-05-25")

        assert calls[0]["start_date"] == "2024-05-24"
        assert calls[0]["end_date"] == "2024-05-25"
        assert calls[0]["include_faa_paths"] is False
        assert calls[0]["study_area"] is vis.study_area

    def test_defaults_to_deployment_year_dates(self, monkeypatch: pytest.MonkeyPatch) -> None:
        vis = self._stub_visualizer()
        calls: list[dict] = []
        monkeypatch.setattr(vis, "_add_track_line", lambda polyline, *, color: MagicMock())
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.get_deployment",
            lambda *args, **kwargs: MagicMock(lat=58.4, lon=-136.0),
        )
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.load_tracks",
            lambda *args, **kwargs: calls.append(kwargs) or self._loaded_tracks(flight=True),
        )

        vis.plot_tracks(TrackSource.ADSB, None, None)

        assert calls[0]["start_date"] == "2024-01-01"
        assert calls[0]["end_date"] == "2024-12-31"

    def test_missing_config_reports_config_error(self, monkeypatch: pytest.MonkeyPatch) -> None:
        vis = self._stub_visualizer()
        messages: list[str] = []
        vis._session.status = lambda msg: messages.append(msg)
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.get_deployment",
            lambda *args, **kwargs: MagicMock(lat=58.4, lon=-136.0),
        )
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.load_tracks",
            lambda *args, **kwargs: (_ for _ in ()).throw(NoOptionError("adsb", "data")),
        )

        vis.plot_tracks(TrackSource.ADSB, "2024-05-24", "2024-05-24")

        assert any("missing config option" in msg for msg in messages)

    def test_parse_key_error_is_not_labeled_missing_config(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vis = self._stub_visualizer()
        messages: list[str] = []
        vis._session.status = lambda msg: messages.append(msg)
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.get_deployment",
            lambda *args, **kwargs: MagicMock(lat=58.4, lon=-136.0),
        )
        monkeypatch.setattr(
            "nps_active_space.viz.tracks_plot.load_tracks",
            lambda *args, **kwargs: (_ for _ in ()).throw(
                KeyError("MXAK AIS CSV lacks 'TIME' column")
            ),
        )

        vis.plot_tracks(TrackSource.AIS, "2024-05-24", "2024-05-24")

        assert messages[-1].startswith("No AIS tracks loaded:")
        assert "TIME" in messages[-1]
        assert not any("missing config option" in msg for msg in messages)


class TestAamActiveSpaceViz:
    @staticmethod
    def _stub_visualizer() -> Visualizer:
        vis = Visualizer.__new__(Visualizer)
        vis.unit, vis.site, vis.year = "DENA", "TRLA", 2024
        vis.project_dir = "/proj"
        vis.crs = "epsg:32608"
        vis.fill_layers = False
        vis.logger = MagicMock()
        vis._status = vis.logger.info
        vis.plotter = MagicMock()
        vis.nmsim_activespace_color = "orange"
        vis.aam_activespace_color = "cyan"
        vis.activespace_color = "orange"
        vis.vessel_track_color = "magenta"
        return _attach_viz_plotters(vis)

    def test_model_display_names(self) -> None:
        assert Visualizer._model_display_name(AcousticModel.AAM) == "AAM"
        assert Visualizer._model_display_name(AcousticModel.NMSIM) == "NMSim"

    def test_ais_track_color_is_not_aam_cyan(self) -> None:
        vis = self._stub_visualizer()
        assert vis.vessel_track_color != vis.aam_activespace_color

    def test_plot_activespace_loads_requested_model(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vis = self._stub_visualizer()
        loaded_models: list[AcousticModel] = []

        class _Active:
            layer_dirs = {500: "/layer"}
            activespaces = {500: gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)])}

        monkeypatch.setattr(
            vis._active_space, "resolve_gain", lambda model, gain: 12.0,
        )
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.load_layered_activespace",
            lambda *args, **kwargs: loaded_models.append(kwargs["model"]) or _Active(),
        )
        monkeypatch.setattr(vis._active_space, "plot_contoured", lambda *args, **kwargs: 0)
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.p.site_dir", lambda *args, **kwargs: "/site"
        )
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.p.model_activespaces_dir",
            lambda *args, **kwargs: "/as",
        )

        vis.plot_activespace(model=AcousticModel.AAM)

        assert loaded_models == [AcousticModel.AAM]
        assert vis._legend_models == [("AAM", "cyan")]

    def test_plot_compare_loads_both_models(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vis = self._stub_visualizer()
        loaded_models: list[AcousticModel] = []

        class _Active:
            layer_dirs = {500: "/layer"}
            activespaces = {500: gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)])}

        monkeypatch.setattr(
            vis._active_space, "resolve_gain", lambda model, gain: 10.0,
        )
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.load_layered_activespace",
            lambda *args, **kwargs: loaded_models.append(kwargs["model"]) or _Active(),
        )
        monkeypatch.setattr(vis._active_space, "plot_contoured", lambda *args, **kwargs: 1)

        vis.plot_compare_activespaces()

        assert loaded_models == [AcousticModel.NMSIM, AcousticModel.AAM]
        assert vis._legend_models == [("NMSim", "orange"), ("AAM", "cyan")]

    def test_plot_compare_ignores_explicit_gain(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vis = self._stub_visualizer()
        gain_args: list[float | None] = []

        class _Active:
            layer_dirs = {500: "/layer"}
            activespaces = {500: gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)])}

        def _resolve(model: AcousticModel, gain: float | None) -> float:
            gain_args.append(gain)
            return 10.0

        monkeypatch.setattr(vis._active_space, "resolve_gain", _resolve)
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.load_layered_activespace",
            lambda *args, **kwargs: _Active(),
        )
        monkeypatch.setattr(vis._active_space, "plot_contoured", lambda *args, **kwargs: 1)

        vis.plot_compare_activespaces(gain=5.0)

        assert gain_args == [None, None]

    def test_fill_layers_skips_empty_polygons(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vis = self._stub_visualizer()
        vis.fill_layers = True
        active = gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)])
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.active_to_polys", lambda _active: []
        )
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.active_to_linestrings",
            lambda _active: [LineString([(0, 0), (1, 1)])],
        )
        monkeypatch.setattr(vis, "_add_labeled_checkbox", lambda *args, **kwargs: MagicMock())
        poly_mesh_calls: list[object] = []
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.polygon_to_mesh",
            lambda *args, **kwargs: poly_mesh_calls.append(args) or MagicMock(),
        )
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.create_polyline_3d",
            lambda *args, **kwargs: MagicMock(),
        )

        vis.plot_active_layer(active, elevation=500.0, i=0)

        assert poly_mesh_calls == []

    def test_resolve_gain_uses_fitted_value(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vis = self._stub_visualizer()
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.resolve_3d_fit_gain",
            lambda *args, **kwargs: 7.5,
        )
        monkeypatch.setattr(
            "nps_active_space.viz.active_space_plot.p.fits_csv", lambda _proj: "/proj/fits.csv"
        )
        assert vis._resolve_activespace_gain(AcousticModel.NMSIM, None) == 7.5

    def test_color_legend_dedupes_entries(self) -> None:
        vis = self._stub_visualizer()
        vis._legend_models = [("NMSim", "orange"), ("NMSim", "orange"), ("AAM", "cyan")]
        vis._add_color_legend(compare_models=False)
        vis.plotter.add_legend.assert_called_once()
        entries = vis.plotter.add_legend.call_args[0][0]
        assert entries == [("NMSim", "orange"), ("AAM", "cyan")]


class TestScriptsVizShim:
    def test_scripts_viz_delegates_to_viz_cli(self):
        from nps_active_space.scripts import viz as scripts_viz
        from nps_active_space.viz.cli import main

        assert scripts_viz.main is main
