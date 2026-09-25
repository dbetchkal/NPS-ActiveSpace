import datetime as dt

import pandas as pd
from pandas.testing import assert_frame_equal
from shapely.geometry import LineString

from nps_active_space.ground_truthing.segments import (
    build_annotation_segments,
    collapse_audible_ranges,
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

    def test_empty_input(self):
        assert collapse_audible_ranges([]) == []

    def test_adjacent_ranges_stay_separate(self):
        t0 = dt.datetime(2020, 1, 1, 12, 0, 0)
        ranges = [
            [t0, t0 + dt.timedelta(minutes=3)],
            [t0 + dt.timedelta(minutes=3), t0 + dt.timedelta(minutes=6)],
        ]
        assert collapse_audible_ranges(ranges) == [[t0, t0 + dt.timedelta(minutes=6)]]


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

    def test_default_audible_ranges_not_shared_across_calls(self):
        """Calling without audible_ranges should never leak state between calls."""
        points = make_track_points(3)
        build_annotation_segments("T1", points)
        result = build_annotation_segments("T2", points)
        assert len(result) == 1
        assert result.audible.iat[0] is False

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
