"""Tests for clock drift utilities."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from nps_active_space.utils.clock_drift import (
    DEFAULT_DRIFT_SEC_PER_DAY,
    DEFAULT_POST_RESET_DRIFT_SEC,
    DEFAULT_SEASON_START_DRIFT_SEC,
    _format_missing_time_ranges,
    build_constant_drift_csv,
    correct_clock_drift,
    infer_correction_period,
)


EXAMPLE_DATA = Path(__file__).resolve().parents[2] / "example_data" / "site_projects" / "DENATRLA"


class TestBuildConstantDriftCsv:
    def test_single_segment_includes_start_noon_and_end(self):
        start_dt = pd.Timestamp("2026-06-05 00:00:00")
        end_dt = pd.Timestamp("2026-06-15 23:59:59.900")
        result = build_constant_drift_csv(
            start_dt=start_dt,
            end_dt=end_dt,
            maintenance_times=[],
            drift_sec_per_day=-10.0,
            season_start_drift_sec=-6.0,
        )

        assert result.iloc[0]["Time"] == start_dt
        assert result.iloc[0]["Seconds"] == pytest.approx(-6.0)
        assert pd.Timestamp(result.iloc[-1]["Time"]) == end_dt
        noon_rows = result[result["Time"].dt.hour == 12]
        assert len(noon_rows) == 11

    def test_maintenance_creates_discontinuity(self):
        maintenance = pd.Timestamp("2025-07-18 11:13:00")
        result = build_constant_drift_csv(
            start_dt=pd.Timestamp("2025-06-23 00:00:00"),
            end_dt=pd.Timestamp("2025-07-19 23:59:59.900"),
            maintenance_times=[maintenance],
            drift_sec_per_day=-9.0,
            season_start_drift_sec=-4.0,
            post_reset_drift_sec=-8.0,
        )

        maintenance_rows = result[
            (result["Time"] >= maintenance - pd.Timedelta(milliseconds=100))
            & (result["Time"] <= maintenance)
        ]
        assert len(maintenance_rows) == 2
        before_reset = maintenance_rows.iloc[0]
        after_reset = maintenance_rows.iloc[1]
        assert before_reset["Time"] == maintenance - pd.Timedelta(milliseconds=100)
        assert after_reset["Time"] == maintenance
        assert after_reset["Seconds"] == pytest.approx(-8.0)
        assert before_reset["Seconds"] < after_reset["Seconds"]

    def test_matches_historical_2026_shape_with_tuned_rate(self):
        reference = pd.read_csv(
            EXAMPLE_DATA / "DENATRLA2026_clock_drift_ADSB.csv",
            parse_dates=["Time"],
        )
        built = build_constant_drift_csv(
            start_dt=reference["Time"].iloc[0],
            end_dt=reference["Time"].iloc[-1],
            maintenance_times=[],
            drift_sec_per_day=-10.14,
            season_start_drift_sec=-5.99,
        )
        pd.testing.assert_frame_equal(
            built.reset_index(drop=True),
            reference.reset_index(drop=True),
            check_exact=False,
            rtol=1e-3,
            atol=0.05,
        )

class TestInferCorrectionPeriod:
    def test_uses_nvspl_bounds_and_encompasses_tracks(self):
        nvspl_dates = ["2025-06-23", "2025-09-18"]
        tracks = pd.DataFrame({
            "point_dt": pd.to_datetime([
                "2025-06-23 08:00:00",
                "2025-09-18 20:00:00",
            ]),
        })
        start_dt, end_dt = infer_correction_period(nvspl_dates, tracks)
        assert start_dt == pd.Timestamp("2025-06-23 00:00:00")
        assert end_dt == pd.Timestamp("2025-09-18 23:59:59.900")


class TestFormatMissingTimeRanges:
    def test_single_contiguous_gap(self):
        missing = pd.date_range("2025-06-23 12:00:00", "2025-06-23 12:00:05", freq="s", inclusive="left")
        summary = _format_missing_time_ranges(missing)
        assert "2025-06-23 12:00:00" in summary
        assert "2025-06-23 12:00:05" in summary

    def test_empty_returns_none(self):
        assert _format_missing_time_ranges(pd.DatetimeIndex([])) == "none"


def _make_tracks(times: list[str]) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"point_dt": pd.to_datetime(times)},
        geometry=[Point(0, 0) for _ in times],
    )


def _write_drift_file(path: Path) -> None:
    path.write_text(
        "Time,Seconds\n"
        "2020-01-01 00:00:00,10\n"
        "2020-01-01 01:00:00,20\n"
    )


class TestCorrectClockDrift:
    def test_inplace_false_returns_corrected_copy(self, tmp_path: Path):
        # drift is 10s at 00:00 and 20s at 01:00, so a point at 00:30 should
        # pick up the interpolated 15s.
        drift_file = tmp_path / "drift.csv"
        _write_drift_file(drift_file)

        tracks = _make_tracks(["2020-01-01 00:00:00", "2020-01-01 00:30:00"])
        original = tracks["point_dt"].tolist()

        result = correct_clock_drift(tracks, str(drift_file), inplace=False)

        expected = pd.to_datetime(["2020-01-01 00:00:10", "2020-01-01 00:30:15"])
        assert result["point_dt"].tolist() == list(expected)
        # inplace=False shouldn't touch the caller's frame
        assert tracks["point_dt"].tolist() == original

    def test_inplace_true_updates_tracks(self, tmp_path: Path):
        drift_file = tmp_path / "drift.csv"
        _write_drift_file(drift_file)

        tracks = _make_tracks(["2020-01-01 00:00:00", "2020-01-01 00:30:00"])
        result = correct_clock_drift(tracks, str(drift_file), inplace=True)

        expected = pd.to_datetime(["2020-01-01 00:00:10", "2020-01-01 00:30:15"])
        assert result is tracks
        assert tracks["point_dt"].tolist() == list(expected)

    def test_raises_when_drift_file_does_not_cover_tracks(self, tmp_path: Path):
        drift_file = tmp_path / "drift.csv"
        drift_file.write_text(
            "Time,Seconds\n"
            "2020-01-01 00:00:00,10\n"
            "2020-01-01 00:30:00,15\n"
        )
        tracks = _make_tracks(["2020-01-01 00:00:00", "2020-01-01 01:00:00"])

        with pytest.raises(Exception, match="Clock drift corrections must encompass"):
            correct_clock_drift(tracks, str(drift_file), inplace=False)
