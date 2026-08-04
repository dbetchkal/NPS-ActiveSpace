"""Tests for the run_clock_drift CLI script."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from nps_active_space.scripts.run_clock_drift import (
    build_parser,
    format_drift_status,
    parse_indices,
    parse_maintenance_times,
    print_drift_summary,
)


class TestRunClockDriftHelpers:
    def test_module_imports(self):
        import nps_active_space.scripts.run_clock_drift  # noqa: F401

    def test_parse_indices(self):
        assert parse_indices("0, 2, 5") == [0, 2, 5]

    def test_parse_maintenance_times(self):
        times = parse_maintenance_times("2025-06-01, 2025-07-15")
        assert len(times) == 2
        assert times[0] == pd.Timestamp("2025-06-01")

    def test_parse_maintenance_times_empty(self):
        assert parse_maintenance_times(None) == []

    @pytest.mark.parametrize(
        ("drift", "expected"),
        [
            (1.5, "valid"),
            (np.nan, "invalid"),
            (None, "invalid"),
        ],
    )
    def test_format_drift_status(self, drift, expected):
        assert format_drift_status(drift) == expected

    def test_print_drift_summary_includes_valid_and_invalid(self, capsys):
        times = pd.to_datetime(["2025-06-23 12:00:00", "2025-06-24 12:00:00"])
        drifts = np.array([2.5, np.nan])
        print_drift_summary(times, drifts)
        output = capsys.readouterr().out
        assert "index" in output
        assert "valid" in output
        assert "invalid" in output
        assert "2.50" in output


class TestRunClockDriftArgparse:
    def test_default_track_source_is_adsb(self):
        parser = build_parser()
        args = parser.parse_args(
            ["-e", "DENA_example", "-u", "DENA", "-s", "TRLA", "-y", "2025"]
        )
        assert args.track_source.value == "ADSB"

    def test_fit_requires_indices(self):
        from nps_active_space.scripts.run_clock_drift import main

        with pytest.raises(SystemExit):
            main(
                [
                    "-e", "DENA_example", "-u", "DENA", "-s", "TRLA", "-y", "2025",
                    "--fit",
                ]
            )
