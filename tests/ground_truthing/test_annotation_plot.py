import datetime as dt

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from matplotlib.dates import date2num
from shapely.geometry import Point

from nps_active_space.ground_truthing.annotation_plot import (
    _cursor_status_text,
    _point_altitude_m,
    _track_altitude_summary,
    closest_spline_at_cursor,
    snap_threshold_days,
)
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
            "Closest approach: 125 m MSL\n"
        )

    def test_returns_empty_when_no_altitude(self):
        points = make_track_points(2)
        assert _track_altitude_summary(points, Point(0.0, 0.0)) == ""


class TestClosestSplineAtCursor:
    def test_snaps_to_track_end_when_cursor_is_in_padding(self):
        spline_time_num = date2num(
            pd.date_range("2020-01-01 12:00:00", periods=5, freq="s")
        )
        cursor_x = spline_time_num[-1] + (120.0 / 86400.0)
        closest_idx, delta_days = closest_spline_at_cursor(spline_time_num, cursor_x)
        assert closest_idx == len(spline_time_num) - 1
        assert delta_days == pytest.approx(120.0 / 86400.0)


class TestSnapThresholdDays:
    def test_defaults_to_one_second_when_spacing_unknown(self):
        assert snap_threshold_days(pd.NaT) == pytest.approx(1.0 / 86400.0)

    def test_converts_median_spacing(self):
        assert snap_threshold_days(pd.Timedelta(seconds=2)) == pytest.approx(2.0 / 86400.0)
