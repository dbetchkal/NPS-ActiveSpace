from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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

DROPPED_METADATA_COLS = [
    "INVID",
    "INSID",
    "GChar2",
    "AdjustmentsApplied",
    "CalibrationAdjustment",
    "GPSTimeAdjustment",
]


def _write_nvspl_copy(
    path: Path,
    source: Path = HOUR_FILE,
    invid: str | None = None,
) -> None:
    """Copy a bundled hour file, optionally setting INVID on the first data row."""
    with source.open(newline="") as f_in:
        reader = csv.reader(f_in)
        rows = list(reader)

    if invid is not None:
        invid_idx = rows[0].index("INVID")
        rows[1][invid_idx] = invid

    with path.open("w", newline="") as f_out:
        csv.writer(f_out).writerows(rows)


@pytest.fixture
def nvspl_hour_file() -> Nvspl:
    return Nvspl([str(HOUR_FILE)])


class TestNvsplFinalizeAfterConcat:
    def test_standard_numeric_columns_are_float32(self, nvspl_hour_file: Nvspl):
        expected = set(Nvspl._NUMERIC_COLS) & set(nvspl_hour_file.columns)
        assert expected, "fixture should include standard numeric NVSPL columns"
        for col in sorted(expected):
            assert nvspl_hour_file[col].dtype == np.float32, col

    def test_drops_all_nan_metadata_columns(self, nvspl_hour_file: Nvspl):
        for col in DROPPED_METADATA_COLS:
            assert col not in nvspl_hour_file.columns
        for col in ["dbA", "WindSpeed", "12.5"]:
            assert col in nvspl_hour_file.columns

    def test_drops_metadata_column_empty_across_multiple_files(self, tmp_path: Path):
        file_a = tmp_path / "nvspl_a.txt"
        file_b = tmp_path / "nvspl_b.txt"
        _write_nvspl_copy(file_a)
        _write_nvspl_copy(file_b)

        nvspl = Nvspl([str(file_a), str(file_b)])
        assert "INVID" not in nvspl.columns

    def test_preserves_metadata_column_populated_in_one_of_two_files(self, tmp_path: Path):
        file_a = tmp_path / "nvspl_a.txt"
        file_b = tmp_path / "nvspl_b.txt"
        _write_nvspl_copy(file_a)
        _write_nvspl_copy(file_b, invid="REC-001")

        nvspl = Nvspl([str(file_a), str(file_b)])
        assert "INVID" in nvspl.columns
        assert nvspl["INVID"].eq("REC-001").any()

    def test_preserves_nonempty_metadata_columns(self, nvspl_hour_file: Nvspl):
        assert "GChar1" in nvspl_hour_file.columns
        assert nvspl_hour_file["GChar1"].notna().any()
        assert "GChar3" in nvspl_hour_file.columns
        assert nvspl_hour_file["GChar3"].notna().any()
        assert "GainAdjustment" in nvspl_hour_file.columns
        assert nvspl_hour_file["GainAdjustment"].dtype == np.int64
        assert "Status" in nvspl_hour_file.columns
        assert nvspl_hour_file["Status"].dtype == np.int64

    def test_coerces_negative_infinity_strings(self, tmp_path: Path):
        rows = list(csv.reader(HOUR_FILE.open(newline="")))
        header = rows[0]
        h12p5_idx = header.index("H12p5")
        dba_idx = header.index("dbA")
        rows[1][h12p5_idx] = "-Infinity"
        expected_dba = float(rows[1][dba_idx])

        nvspl_path = tmp_path / "nvspl_negative_infinity.txt"
        with nvspl_path.open("w", newline="") as f_out:
            csv.writer(f_out).writerows(rows[:2])

        nvspl = Nvspl([str(nvspl_path)])
        assert float(nvspl.loc[nvspl.index[0], "12.5"]) == -np.inf
        assert nvspl["12.5"].dtype == np.float32
        assert nvspl.loc[nvspl.index[0], "dbA"] == expected_dba
