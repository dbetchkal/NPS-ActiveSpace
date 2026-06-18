from unittest.mock import patch

import numpy as np
import pytest

from nps_active_space.utils.metrics import (
    _ais_dates_from_archive,
    _full_causal_days,
    get_obs_periods,
)


def _touch_csv(directory, name: str) -> None:
    (directory / name).write_text("MMSI,TIME,LAT,LON\n", encoding="utf-8")


class TestAisDatesFromArchive:
    def test_extracts_dates_from_mxak_filenames(self, tmp_path):
        _touch_csv(tmp_path, "MXAK-AIS-Anchorage-20240601.csv")
        _touch_csv(tmp_path, "MXAK-AIS-Juneau-20240603-1.csv")
        _touch_csv(tmp_path, "MXAK-AIS-Seward-20240602.csv")

        dates = _ais_dates_from_archive(str(tmp_path))

        assert dates.tolist() == [
            np.datetime64("2024-06-01"),
            np.datetime64("2024-06-02"),
            np.datetime64("2024-06-03"),
        ]

    def test_ignores_non_matching_csv_files(self, tmp_path):
        _touch_csv(tmp_path, "MXAK-AIS-Anchorage-20240601.csv")
        _touch_csv(tmp_path, "notes.csv")

        dates = _ais_dates_from_archive(str(tmp_path))

        assert dates.tolist() == [np.datetime64("2024-06-01")]


class TestFullCausalDays:
    def test_marks_interior_consecutive_days_as_full(self):
        dates = np.array(
            ["2024-06-01", "2024-06-02", "2024-06-03", "2024-06-05"],
            dtype="datetime64[D]",
        )

        full_day = _full_causal_days(dates)

        assert full_day.tolist() == [False, True, False, False]


class TestGetObsPeriodsAis:
    @patch("nps_active_space.utils.metrics.Srcid")
    @patch("nps_active_space.utils.metrics.iyore.Dataset")
    def test_ais_only_uses_full_interior_days(self, mock_dataset, mock_srcid, tmp_path):
        _touch_csv(tmp_path, "MXAK-AIS-Anchorage-20240601.csv")
        _touch_csv(tmp_path, "MXAK-AIS-Anchorage-20240602.csv")
        _touch_csv(tmp_path, "MXAK-AIS-Anchorage-20240603.csv")

        mock_dataset.return_value.srcid.return_value = []

        obs_periods = get_obs_periods(
            "DENA",
            "TRLA",
            "2024",
            "/fake/nvspl",
            ais_path=str(tmp_path),
        )

        assert obs_periods.tolist() == [["2024-06-02", "2024-06-02"]]
