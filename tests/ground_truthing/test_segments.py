import datetime as dt

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from shapely.geometry import LineString

from nps_active_space.ground_truthing.segments import (
    ANNOTATION_MAX_VERTICES,
    build_annotation_segments,
    collapse_audible_ranges,
    downsample_linestring,
)
from helpers import SEGMENT_TABULAR_COLS, _tabular, make_track_points


class TestCollapseAudibleRanges:
    def test_merges_overlap(self):
        t0 = dt.datetime(2020, 1, 1, 12, 0, 0)
        ranges = [
            [t0, t0 + dt.timedelta(minutes=5)],
            [t0 + dt.timedelta(minutes=3), t0 + dt.timedelta(minutes=10)],
        ]
        expected_ranges = [[t0, t0 + dt.timedelta(minutes=10)]]
        assert collapse_audible_ranges(ranges) == expected_ranges

    def test_keeps_separate_non_overlapping(self):
        t0 = dt.datetime(2020, 1, 1, 12, 0, 0)
        ranges = [
            [t0, t0 + dt.timedelta(minutes=2)],
            [t0 + dt.timedelta(minutes=5), t0 + dt.timedelta(minutes=7)],
        ]
        assert collapse_audible_ranges(ranges) == ranges

    def test_single_range_unchanged(self):
        t0 = dt.datetime(2020, 1, 1, 12, 0, 0)
        ranges = [[t0, t0 + dt.timedelta(minutes=3)]]
        assert collapse_audible_ranges(ranges) == ranges


class TestBuildAnnotationSegments:
    def test_all_inaudible_when_no_ranges(self):
        points = make_track_points(3)
        result = build_annotation_segments("T1", points, audible_ranges=[], valid=True)
        actual = _tabular(result, SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1"],
                "start_dt": [points.point_dt.iat[0]],
                "end_dt": [points.point_dt.iat[-1]],
                "valid": [True],
                "audible": [False],
                "note": [None],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_invalid_track_single_segment(self):
        points = make_track_points(4)
        result = build_annotation_segments(
            "T1",
            points,
            audible_ranges=[[points.point_dt.iat[1], points.point_dt.iat[2]]],
            valid=False,
        )
        actual = _tabular(result, SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1"],
                "start_dt": [points.point_dt.iat[0]],
                "end_dt": [points.point_dt.iat[-1]],
                "valid": [False],
                "audible": [False],
                "note": [None],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_tail_inaudible_after_audible_window(self):
        points = make_track_points(6)
        t = points.point_dt
        audible_ranges = [[t.iat[1], t.iat[3]]]
        result = build_annotation_segments("T1", points, audible_ranges=audible_ranges)
        actual = _tabular(result.sort_values("start_dt"), SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1", "T1"],
                "start_dt": [t.iat[1], t.iat[3]],
                "end_dt": [t.iat[2], t.iat[4]],
                "valid": [True, True],
                "audible": [True, False],
                "note": [None, None],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)
        assert all(isinstance(geom, LineString) for geom in result.geometry)

    def test_two_audible_windows_with_gap(self):
        points = make_track_points(8)
        t = points.point_dt
        audible_ranges = [[t.iat[1], t.iat[3]], [t.iat[5], t.iat[7]]]
        result = build_annotation_segments("T1", points, audible_ranges=audible_ranges)
        actual = _tabular(result.sort_values("start_dt"), SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1", "T1", "T1"],
                "start_dt": [t.iat[1], t.iat[3], t.iat[5]],
                "end_dt": [t.iat[2], t.iat[4], t.iat[6]],
                "valid": [True, True, True],
                "audible": [True, False, True],
                "note": [None, None, None],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_overlapping_input_ranges_collapsed_before_split(self):
        points = make_track_points(6)
        t = points.point_dt
        audible_ranges = [[t.iat[1], t.iat[3]], [t.iat[2], t.iat[4]]]
        result = build_annotation_segments("T1", points, audible_ranges=audible_ranges)
        actual = _tabular(result, SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1"],
                "start_dt": [t.iat[1]],
                "end_dt": [t.iat[3]],
                "valid": [True],
                "audible": [True],
                "note": [None],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_splits_on_time_audible_not_point_dt(self):
        point_dt = pd.date_range("2020-01-01 12:00", periods=6, freq="min")
        time_audible = point_dt + pd.Timedelta(seconds=30)
        points = make_track_points(6, time_audible=time_audible)
        audible_ranges = [[time_audible[1], time_audible[3]]]
        result = build_annotation_segments("T1", points, audible_ranges=audible_ranges)
        actual = _tabular(result.sort_values("start_dt"), SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1", "T1"],
                "start_dt": [point_dt[1], point_dt[3]],
                "end_dt": [point_dt[2], point_dt[4]],
                "valid": [True, True],
                "audible": [True, False],
                "note": [None, None],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_note_propagates_to_segments(self):
        points = make_track_points(6)
        t = points.point_dt
        result = build_annotation_segments(
            "T1",
            points,
            audible_ranges=[[t.iat[1], t.iat[3]]],
            note="test note",
        )
        actual = _tabular(result.sort_values("start_dt"), SEGMENT_TABULAR_COLS)
        expected = pd.DataFrame(
            {
                "_id": ["T1", "T1"],
                "start_dt": [t.iat[1], t.iat[3]],
                "end_dt": [t.iat[2], t.iat[4]],
                "valid": [True, True],
                "audible": [True, False],
                "note": ["test note", "test note"],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_long_track_geometry_downsampled(self):
        points = make_track_points(500, freq="s")
        t = points.point_dt
        result = build_annotation_segments(
            "T1",
            points,
            audible_ranges=[[t.iat[0], t.iat[-1]]],
        )
        geom = result.iloc[0].geometry
        assert isinstance(geom, LineString)
        assert len(geom.coords) <= ANNOTATION_MAX_VERTICES
        assert geom.coords[0][0] == pytest.approx(0.0, abs=1e-6)
        assert geom.coords[0][1] == pytest.approx(0.0, abs=1e-6)
        assert geom.coords[-1][0] == pytest.approx(498.0, abs=1e-6)
        assert geom.coords[-1][1] == pytest.approx(498.0, abs=1e-6)


class TestDownsampleLinestring:
    def test_caps_vertices_and_preserves_endpoints(self):
        coords = [(float(i), float(i), float(i)) for i in range(1000)]
        line = LineString(coords)
        out = downsample_linestring(line, max_vertices=32)
        assert len(out.coords) == 32
        assert out.coords[0] == pytest.approx(coords[0], abs=1e-6)
        assert out.coords[-1] == pytest.approx(coords[-1], abs=1e-6)

    def test_short_line_only_rounds(self):
        line = LineString([(1.123456789, 2.987654321, 100.55), (3.0, 4.0, 200.44)])
        out = downsample_linestring(line)
        assert len(out.coords) == 2
        assert out.coords[0] == (1.123457, 2.987654, 100.5)
        assert out.coords[1] == (3.0, 4.0, 200.4)
