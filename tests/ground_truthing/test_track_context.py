import datetime as dt

import geopandas as gpd
import pandas as pd
import pytest
from matplotlib.dates import date2num, num2date

from nps_active_space.ground_truthing.track_context import (
    audible_ranges_from_annotations,
    limit_line_bounds,
)
from helpers import make_track_points


def _annotations(
    spline: gpd.GeoDataFrame,
    *,
    valid: list[bool],
    audible: list[bool],
    start_dt: list[pd.Timestamp],
    end_dt: list[pd.Timestamp],
) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "valid": valid,
            "audible": audible,
            "start_dt": start_dt,
            "end_dt": end_dt,
        },
        geometry=[spline.geometry.iat[i] for i in range(len(valid))],
        crs=spline.crs,
    )


class TestLimitLineBounds:
    def test_centered_window_within_spectrogram(self):
        closest = pd.Timestamp("2020-01-01 12:05:00")
        x_lims = date2num(pd.date_range("2020-01-01 12:00", periods=11, freq="min"))
        expected = (
            date2num(closest - dt.timedelta(seconds=60)),
            date2num(closest + dt.timedelta(seconds=60)),
        )
        assert limit_line_bounds(closest, x_lims) == expected

    def test_clamps_to_spectrogram_edges(self):
        closest = pd.Timestamp("2020-01-01 12:00:30")
        x_lims = date2num(pd.date_range("2020-01-01 12:00", periods=3, freq="min"))
        expected = (x_lims[0], date2num(closest + dt.timedelta(seconds=60)))
        assert limit_line_bounds(closest, x_lims) == expected

    def test_upper_clamped_at_end(self):
        closest = pd.Timestamp("2020-01-01 12:09:30")
        x_lims = date2num(pd.date_range("2020-01-01 12:00", periods=10, freq="min"))
        expected = (date2num(closest - dt.timedelta(seconds=60)), x_lims[-1])
        assert limit_line_bounds(closest, x_lims) == expected


class TestAudibleRangesFromAnnotations:
    def test_empty_annotations_returns_default_slider_range(self):
        spline = make_track_points(5)
        lower = date2num(spline.point_dt.iat[2] - dt.timedelta(seconds=30))
        upper = date2num(spline.point_dt.iat[2] + dt.timedelta(seconds=30))
        annotations = gpd.GeoDataFrame(
            columns=["valid", "audible", "start_dt", "end_dt"],
            geometry=[],
            crs="EPSG:4326",
        )
        expected = [[num2date(lower), num2date(upper)]]
        assert audible_ranges_from_annotations(annotations, spline, lower, upper) == expected

    def test_restores_audible_range_on_time_audible_axis(self):
        spline = make_track_points(6)
        annotations = _annotations(
            spline,
            valid=[True],
            audible=[True],
            start_dt=[spline.point_dt.iat[1]],
            end_dt=[spline.point_dt.iat[4]],
        )
        lower = date2num(spline.point_dt.iat[0])
        upper = date2num(spline.point_dt.iat[-1])
        expected = [[spline.time_audible.iat[1], spline.time_audible.iat[4]]]
        assert audible_ranges_from_annotations(annotations, spline, lower, upper) == expected

    def test_ignores_invalid_or_inaudible_rows(self):
        spline = make_track_points(5)
        annotations = _annotations(
            spline,
            valid=[False, True],
            audible=[True, False],
            start_dt=[spline.point_dt.iat[0], spline.point_dt.iat[1]],
            end_dt=[spline.point_dt.iat[2], spline.point_dt.iat[3]],
        )
        lower = date2num(spline.point_dt.iat[0])
        upper = date2num(spline.point_dt.iat[-1])
        assert audible_ranges_from_annotations(annotations, spline, lower, upper) == []

    def test_maps_annotations_onto_time_audible_axis(self):
        point_dt = pd.date_range("2020-01-01 12:00", periods=6, freq="min")
        time_audible = point_dt + pd.Timedelta(seconds=45)
        spline = make_track_points(6, time_audible=time_audible)
        annotations = _annotations(
            spline,
            valid=[True],
            audible=[True],
            start_dt=[point_dt[2]],
            end_dt=[point_dt[4]],
        )
        lower = date2num(point_dt[0])
        upper = date2num(point_dt[-1])
        expected = [[time_audible[2], time_audible[4]]]
        assert audible_ranges_from_annotations(annotations, spline, lower, upper) == expected

    def test_returns_none_when_start_outside_tolerance(self):
        spline = make_track_points(5)
        annotations = _annotations(
            spline,
            valid=[True],
            audible=[True],
            start_dt=[spline.point_dt.iat[0] - pd.Timedelta(seconds=5)],
            end_dt=[spline.point_dt.iat[2]],
        )
        lower = date2num(spline.point_dt.iat[0])
        upper = date2num(spline.point_dt.iat[-1])
        with pytest.warns(UserWarning, match="annotation start"):
            result = audible_ranges_from_annotations(annotations, spline, lower, upper)
        assert result is None

    def test_returns_none_when_end_outside_tolerance(self):
        spline = make_track_points(5)
        annotations = _annotations(
            spline,
            valid=[True],
            audible=[True],
            start_dt=[spline.point_dt.iat[1]],
            end_dt=[spline.point_dt.iat[-1] + pd.Timedelta(seconds=5)],
        )
        lower = date2num(spline.point_dt.iat[0])
        upper = date2num(spline.point_dt.iat[-1])
        with pytest.warns(UserWarning, match="annotation end"):
            result = audible_ranges_from_annotations(annotations, spline, lower, upper)
        assert result is None

    def test_accepts_subsecond_annotation_drift(self):
        spline = make_track_points(5)
        annotations = _annotations(
            spline,
            valid=[True],
            audible=[True],
            start_dt=[spline.point_dt.iat[1] + pd.Timedelta(milliseconds=500)],
            end_dt=[spline.point_dt.iat[3] - pd.Timedelta(milliseconds=500)],
        )
        lower = date2num(spline.point_dt.iat[0])
        upper = date2num(spline.point_dt.iat[-1])
        expected = [[spline.time_audible.iat[1], spline.time_audible.iat[3]]]
        assert audible_ranges_from_annotations(annotations, spline, lower, upper) == expected
