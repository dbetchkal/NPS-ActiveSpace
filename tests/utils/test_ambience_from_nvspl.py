"""Tests for NVSPL-derived ambience."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from nps_active_space.utils.computation import ambience_from_nvspl, is_usable_spectral_ambience


def _make_nvspl_df(n_rows: int = 100, wind_speed: float = 0.0) -> pd.DataFrame:
    bands = {
        "12.5": 30.0,
        "20": 28.0,
        "1000": 25.0,
        "12500": 20.0,
        "20000": 18.0,
    }
    row = {
        "SiteID": "TEST",
        "dbA": 35.0,
        "dbC": 40.0,
        "Voltage": 13.0,
        "WindSpeed": wind_speed,
        "WindDir": 180.0,
        "TempIns": 20.0,
        "TempOut": 15.0,
        "Humidity": 50.0,
        "INVID": "",
        "INSID": "",
        "GChar1": "",
        "AdjustmentsApplied": "",
        "CalibrationAdjustment": 0,
        "GPSTimeAdjustment": 0,
        "GainAdjustment": 0,
        "Status": 0,
        **bands,
    }
    return pd.DataFrame([row] * n_rows)


class TestAmbienceFromNvspl:
    def test_missing_wind_skips_filter(self):
        ambience = ambience_from_nvspl(_make_nvspl_df(wind_speed=np.nan), quantile=90, broadband=False)
        assert ambience.loc["1000"] == pytest.approx(25.0, abs=0.01)
        assert is_usable_spectral_ambience(ambience)

    def test_all_high_wind_falls_back_to_all_rows(self):
        ambience = ambience_from_nvspl(_make_nvspl_df(wind_speed=10.0), quantile=90, broadband=False)
        assert ambience.loc["1000"] == pytest.approx(25.0, abs=0.01)

    def test_is_usable_spectral_ambience_rejects_all_nan(self):
        ambience = pd.Series({"20": np.nan, "1000": np.nan, "12500": np.nan})
        assert not is_usable_spectral_ambience(ambience)
