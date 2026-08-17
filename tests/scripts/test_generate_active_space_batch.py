import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from nps_active_space.scripts.generate_active_space_batch import (
    RESULT_COLUMNS,
    read_results_file,
    run_deployment,
)


STUB_SCRIPT = """
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("--results-out", required=True)
args = parser.parse_args()

with open(args.results_out, "w") as results_file:
    json.dump({
        "Number of valid annotated segments": 5,
        "Mean altitude": 3000,
        "KDE reduction (%)": "12.5%",
        "1/3rd Octave Gain (F1)": 12.5,
        "F1": 0.91,
    }, results_file)
"""


class TestReadResultsFile:
    def test_read_results_file(self, tmp_path: Path):
        results_path = tmp_path / "results.json"
        results_path.write_text(json.dumps({
            "Number of valid annotated segments": 12,
            "Mean altitude": 3000,
            "KDE reduction (%)": "15.0%",
            "1/3rd Octave Gain (F1)": 12.5,
            "F1": 0.87,
        }))

        series = read_results_file(results_path, "DENATRLA20253000m")

        expected = pd.Series({
            "Designator": "DENATRLA20253000m",
            "Number of valid annotated segments": 12,
            "Mean altitude": 3000,
            "KDE reduction (%)": "15.0%",
            "1/3rd Octave Gain (F1)": 12.5,
            "F1": 0.87,
        })
        pd.testing.assert_series_equal(series, expected)


class TestRunDeployment:
    def _stub_cmd(self, tmp_path: Path) -> list[str]:
        stub_path = tmp_path / "stub_generate_active_space.py"
        stub_path.write_text(STUB_SCRIPT)
        return [sys.executable, str(stub_path)]

    def test_run_deployment_reads_stub_script_results(self, tmp_path: Path):
        series = run_deployment("DENATRLA20253000m", self._stub_cmd(tmp_path))

        assert series is not None
        assert list(series.index) == RESULT_COLUMNS
        assert series["Designator"] == "DENATRLA20253000m"
        assert series["Number of valid annotated segments"] == 5
        assert series["Mean altitude"] == 3000
        assert series["KDE reduction (%)"] == "12.5%"
        assert series["1/3rd Octave Gain (F1)"] == 12.5
        assert series["F1"] == 0.91

    def test_run_deployment_returns_none_on_nonzero_exit(self, tmp_path: Path, monkeypatch):
        monkeypatch.setattr(
            "nps_active_space.scripts.generate_active_space_batch.subprocess.run",
            MagicMock(return_value=MagicMock(returncode=1)),
        )
        result = run_deployment("DENATRLA20253000m", self._stub_cmd(tmp_path))
        assert result is None

    def test_run_deployment_returns_none_when_results_missing(self, tmp_path: Path, monkeypatch):
        monkeypatch.setattr(
            "nps_active_space.scripts.generate_active_space_batch.subprocess.run",
            MagicMock(return_value=MagicMock(returncode=0)),
        )
        result = run_deployment("DENATRLA20253000m", self._stub_cmd(tmp_path))
        assert result is None

    def test_run_deployment_cleans_up_temp_results_file(self, tmp_path: Path):
        run_deployment("DENATRLA20253000m", self._stub_cmd(tmp_path))
        leftover = list(tmp_path.glob("*.json"))
        assert leftover == []


class TestBatchCsvIntegration:
    def test_run_deployment_results_write_expected_csv_row(self, tmp_path: Path):
        stub_path = tmp_path / "stub_generate_active_space.py"
        stub_path.write_text(STUB_SCRIPT)
        output_csv = tmp_path / "batch_output.csv"

        result_series = run_deployment(
            "DENATRLA20253000m",
            [sys.executable, str(stub_path)],
        )
        assert result_series is not None

        output_df = result_series.to_frame().T
        output_df.to_csv(output_csv, index=False)

        loaded = pd.read_csv(output_csv)
        row = loaded.iloc[0]
        assert row["Designator"] == "DENATRLA20253000m"
        assert row["Number of valid annotated segments"] == 5
        assert row["Mean altitude"] == 3000
        assert row["KDE reduction (%)"] == "12.5%"
        assert row["1/3rd Octave Gain (F1)"] == 12.5
        assert row["F1"] == 0.91
