from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from nps_active_space.utils.models import Nvspl

REPO = Path(__file__).resolve().parents[2]
NVSPL_DIR = (
    REPO
    / "example_data"
    / "nvspl_archive"
    / "2024 GLBALSTL Lower South Tidal Inlet"
    / "01 DATA"
    / "NVSPL"
)
HOUR_FILE = NVSPL_DIR / "NVSPL_GLBALSTL_2024_05_24_00.txt"

NUMERIC_COLS = [
    "12.5",
    "1000",
    "20000",
    "dbA",
    "dbC",
    "WindSpeed",
    "WindDir",
]
DROPPED_METADATA_COLS = [
    "INVID",
    "INSID",
    "GChar2",
    "AdjustmentsApplied",
    "CalibrationAdjustment",
    "GPSTimeAdjustment",
]


class TestNvsplFileRead:
    def test_numeric_columns_are_float32(self):
        nvspl = Nvspl([str(HOUR_FILE)])
        for col in NUMERIC_COLS:
            assert nvspl[col].dtype == np.float32

    def test_drops_all_nan_metadata_columns(self):
        nvspl = Nvspl([str(HOUR_FILE)])
        for col in DROPPED_METADATA_COLS:
            assert col not in nvspl.columns
        for col in ["dbA", "WindSpeed", "12.5"]:
            assert col in nvspl.columns

    def test_preserves_nonempty_metadata_columns(self):
        nvspl = Nvspl([str(HOUR_FILE)])
        assert "GChar1" in nvspl.columns
        assert nvspl["GChar1"].notna().any()
        assert "GainAdjustment" in nvspl.columns
        assert "Status" in nvspl.columns

    def test_coerces_negative_infinity_strings(self, tmp_path: Path):
        source = HOUR_FILE.read_text()
        header, first_row = source.splitlines()[:2]
        mutated_row = first_row.replace("41.1", "-Infinity", 1)
        nvspl_path = tmp_path / "nvspl_negative_infinity.txt"
        nvspl_path.write_text(f"{header}\n{mutated_row}\n")

        nvspl = Nvspl([str(nvspl_path)])
        value = float(nvspl.loc[nvspl.index[0], "12.5"])
        assert value == -np.inf
        assert nvspl["12.5"].dtype == np.float32
