"""Tests for the single-event peak-matching clock drift method.

These exercise the standalone helper functions directly with synthetic data,
since `ClockDriftFixer` itself requires a full site project (study area, FAA
databases, etc.) to instantiate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from nps_active_space.utils.clock_drift import (
    _find_acoustic_peaks,
    _format_missing_time_ranges,
    _predicted_flight_peaks,
    match_single_event,
)


def _make_sobel_signal(peak_times_and_heights: list[tuple[str, float]], baseline: float = 40.0) -> pd.Series:
    """Build a synthetic per-second Sobel signal with sharp peaks at given times."""
    index = pd.date_range("2025-06-23 00:00:00", "2025-06-24 00:00:00", freq="s", inclusive="left")
    values = np.full(len(index), baseline)
    signal = pd.Series(values, index=index)
    for time_str, height in peak_times_and_heights:
        signal.loc[pd.Timestamp(time_str)] = height
    return signal


def _make_pts_section(flights: list[tuple[str, str, float]]) -> pd.DataFrame:
    """Build a synthetic pts_section with one predicted-peak row per flight.

    Each flight is (track_id, time_audible, Lp_est).
    """
    return pd.DataFrame(
        {
            "track_id": [f[0] for f in flights],
            "time_audible": pd.to_datetime([f[1] for f in flights]),
            "Lp_est": [f[2] for f in flights],
        }
    )


class TestFindAcousticPeaks:
    def test_finds_prominent_isolated_peak(self):
        signal = _make_sobel_signal([("2025-06-23 12:00:00", 80.0)])
        peaks = _find_acoustic_peaks(signal, min_prominence_db=6.0, min_isolation_sec=30)
        assert len(peaks) == 1
        assert peaks.loc[0, "time"] == pd.Timestamp("2025-06-23 12:00:00")

    def test_no_peaks_below_prominence_threshold(self):
        # a small bump that doesn't clear the prominence threshold
        signal = _make_sobel_signal([("2025-06-23 12:00:00", 42.0)])
        peaks = _find_acoustic_peaks(signal, min_prominence_db=6.0, min_isolation_sec=30)
        assert peaks.empty

    def test_no_peaks_on_flat_signal(self):
        signal = _make_sobel_signal([])
        peaks = _find_acoustic_peaks(signal, min_prominence_db=6.0, min_isolation_sec=30)
        assert peaks.empty

    def test_sorted_by_descending_prominence(self):
        signal = _make_sobel_signal([
            ("2025-06-23 08:00:00", 70.0),
            ("2025-06-23 16:00:00", 90.0),
        ])
        peaks = _find_acoustic_peaks(signal, min_prominence_db=6.0, min_isolation_sec=30)
        assert len(peaks) == 2
        assert peaks.loc[0, "time"] == pd.Timestamp("2025-06-23 16:00:00")
        assert peaks["prominence"].is_monotonic_decreasing


class TestPredictedFlightPeaks:
    def test_one_row_per_flight_at_max_lp(self):
        pts_section = pd.DataFrame({
            "track_id": ["A", "A", "B"],
            "time_audible": pd.to_datetime([
                "2025-06-23 12:00:00", "2025-06-23 12:00:05", "2025-06-23 15:00:00",
            ]),
            "Lp_est": [30.0, 55.0, 40.0],
        })
        flight_peaks = _predicted_flight_peaks(pts_section)
        flight_peaks = flight_peaks.set_index("track_id")
        assert flight_peaks.loc["A", "time"] == pd.Timestamp("2025-06-23 12:00:05")
        assert flight_peaks.loc["A", "Lp_est"] == 55.0
        assert flight_peaks.loc["B", "time"] == pd.Timestamp("2025-06-23 15:00:00")

    def test_empty_input_returns_empty(self):
        pts_section = pd.DataFrame(columns=["track_id", "time_audible", "Lp_est"])
        flight_peaks = _predicted_flight_peaks(pts_section)
        assert flight_peaks.empty


class TestMatchSingleEvent:
    def test_matches_peak_to_nearby_flight_with_expected_drift(self):
        acoustic_time = pd.Timestamp("2025-06-23 12:00:03")
        predicted_time = pd.Timestamp("2025-06-23 12:00:00")
        signal = _make_sobel_signal([(str(acoustic_time), 80.0)])
        pts_section = _make_pts_section([("A", str(predicted_time), 50.0)])

        drift_sec, match_info = match_single_event(
            signal, pts_section, max_clock_drift=pd.Timedelta(minutes=5),
        )

        assert drift_sec == pytest.approx(3.0)
        assert match_info["track_id"] == "A"
        assert match_info["acoustic_time"] == acoustic_time
        assert match_info["predicted_time"] == predicted_time
        assert match_info["confidence"] == match_info["prominence"]
        assert match_info["confidence"] > 0

    def test_no_match_without_acoustic_peak(self):
        signal = _make_sobel_signal([])
        pts_section = _make_pts_section([("A", "2025-06-23 12:00:00", 50.0)])

        drift_sec, match_info = match_single_event(
            signal, pts_section, max_clock_drift=pd.Timedelta(minutes=5),
        )

        assert drift_sec is None
        assert match_info is None

    def test_no_match_without_nearby_flight(self):
        signal = _make_sobel_signal([("2025-06-23 12:00:00", 80.0)])
        # flight is way outside the max_clock_drift window
        pts_section = _make_pts_section([("A", "2025-06-23 18:00:00", 50.0)])

        drift_sec, match_info = match_single_event(
            signal, pts_section, max_clock_drift=pd.Timedelta(minutes=5),
        )

        assert drift_sec is None
        assert match_info is None

    def test_picks_loudest_candidate_when_multiple_flights_in_window(self):
        acoustic_time = pd.Timestamp("2025-06-23 12:00:00")
        signal = _make_sobel_signal([(str(acoustic_time), 80.0)])
        pts_section = _make_pts_section([
            ("quiet_flight", "2025-06-23 12:00:02", 20.0),
            ("loud_flight", "2025-06-23 12:00:01", 60.0),
        ])

        drift_sec, match_info = match_single_event(
            signal, pts_section, max_clock_drift=pd.Timedelta(minutes=5),
        )

        assert match_info["track_id"] == "loud_flight"
        assert drift_sec == pytest.approx(-1.0)


class TestFormatMissingTimeRanges:
    def test_single_contiguous_gap(self):
        missing = pd.date_range("2025-06-23 12:00:00", "2025-06-23 12:00:05", freq="s", inclusive="left")
        summary = _format_missing_time_ranges(missing)
        assert "2025-06-23 12:00:00" in summary
        assert "2025-06-23 12:00:05" in summary

    def test_empty_returns_none(self):
        assert _format_missing_time_ranges(pd.DatetimeIndex([])) == "none"
