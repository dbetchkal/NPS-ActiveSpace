import datetime as dt

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from matplotlib.dates import date2num
from shapely.geometry import Point

from nps_active_space.ground_truthing.annotation_plot import (
    DEFAULT_SNAP_THRESHOLD_SECONDS,
    _cursor_status_text,
    _point_altitude_m,
    _show_track_altitude,
    _track_altitude_summary,
    closest_spline_at_cursor,
    seconds_to_matplotlib_days,
    snap_threshold_days,
)
from nps_active_space.utils.enums import TrackSource
from helpers import make_track_points


class TestCursorStatusText:
    def test_includes_altitude_when_finite(self):
        cursor_time = dt.datetime(2020, 1, 1, 12, 34, 56)
        assert _cursor_status_text(cursor_time, 842.4) == (
            "Cursor Time: 12:34:56\nAltitude: 842 m MSL"
        )

    def test_omits_altitude_when_missing(self):
        cursor_time = dt.datetime(2020, 1, 1, 12, 34, 56)
        assert _cursor_status_text(cursor_time, None) == "Cursor Time: 12:34:56"

    def test_omits_altitude_when_not_finite(self):
        cursor_time = dt.datetime(2020, 1, 1, 12, 34, 56)
        assert _cursor_status_text(cursor_time, np.nan) == "Cursor Time: 12:34:56"


class TestPointAltitudeM:
    def test_reads_z_from_3d_point(self):
        assert _point_altitude_m(Point(1.0, 2.0, 123.4)) == 123.4

    def test_returns_none_for_2d_point(self):
        assert _point_altitude_m(Point(1.0, 2.0)) is None


class TestShowTrackAltitude:
    def test_hidden_for_ais(self):
        assert _show_track_altitude(TrackSource.AIS) is False

    def test_shown_for_adsb(self):
        assert _show_track_altitude(TrackSource.ADSB) is True


class TestTrackAltitudeSummary:
    def test_summarizes_track_and_closest_altitudes(self):
        points = make_track_points(3)
        points["z"] = [100.0, 200.0, 150.0]
        points.geometry = gpd.points_from_xy(
            points.geometry.x,
            points.geometry.y,
            points["z"],
        )
        summary = _track_altitude_summary(points, Point(0.0, 0.0, 125.0))
        assert summary == (
            "Track altitude: 100–200 m MSL (mean 150)\n"
            "Closest point altitude: 125 m MSL\n"
        )

    def test_returns_empty_when_no_altitude(self):
        points = make_track_points(2)
        assert _track_altitude_summary(points, Point(0.0, 0.0)) == ""


class TestClosestSplineAtCursor:
    def test_snaps_to_track_end_when_cursor_is_in_padding(self):
        spline_time_num = date2num(
            pd.date_range("2020-01-01 12:00:00", periods=5, freq="s")
        )
        padding_seconds = 120.0
        cursor_x = spline_time_num[-1] + seconds_to_matplotlib_days(padding_seconds)
        closest_idx, delta_days = closest_spline_at_cursor(spline_time_num, cursor_x)
        assert closest_idx == len(spline_time_num) - 1
        assert delta_days == pytest.approx(seconds_to_matplotlib_days(padding_seconds))


class TestSnapThresholdDays:
    def test_defaults_to_one_second_when_spacing_unknown(self):
        assert snap_threshold_days(pd.NaT) == pytest.approx(
            seconds_to_matplotlib_days(DEFAULT_SNAP_THRESHOLD_SECONDS)
        )

    def test_converts_median_spacing(self):
        spacing_seconds = 2.0
        assert snap_threshold_days(pd.Timedelta(seconds=spacing_seconds)) == pytest.approx(
            seconds_to_matplotlib_days(spacing_seconds)
        )
