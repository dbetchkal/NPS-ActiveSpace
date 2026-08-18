import json
import re
import shlex
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

import nps_active_space.scripts.generate_active_space_batch as batch_module
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

    def test_read_results_file_returns_none_on_missing_required_key(self, tmp_path: Path):
        results_path = tmp_path / "results.json"
        results_path.write_text(json.dumps({
            "Number of valid annotated segments": 12,
            "Mean altitude": 3000,
            "KDE reduction (%)": "15.0%",
            "1/3rd Octave Gain (F1)": 12.5,
        }))

        assert read_results_file(results_path, "DENATRLA20253000m") is None

    def test_read_results_file_returns_none_on_invalid_json(self, tmp_path: Path):
        results_path = tmp_path / "results.json"
        results_path.write_text("{not valid json")

        assert read_results_file(results_path, "DENATRLA20253000m") is None


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

    def test_run_deployment_returns_none_on_invalid_json_results(self, tmp_path: Path):
        stub_path = tmp_path / "invalid_json_stub.py"
        stub_path.write_text(
            """
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--results-out", required=True)
args = parser.parse_args()
with open(args.results_out, "w") as results_file:
    results_file.write("{not valid json")
"""
        )
        result = run_deployment("DENATRLA20253000m", [sys.executable, str(stub_path)])
        assert result is None

    def test_run_deployment_returns_none_on_missing_required_key(self, tmp_path: Path):
        stub_path = tmp_path / "partial_stub.py"
        stub_path.write_text(
            """
import argparse
import json
parser = argparse.ArgumentParser()
parser.add_argument("--results-out", required=True)
args = parser.parse_args()
with open(args.results_out, "w") as results_file:
    json.dump({"Number of valid annotated segments": 1}, results_file)
"""
        )
        result = run_deployment("DENATRLA20253000m", [sys.executable, str(stub_path)])
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


STUB_SCRIPT_WITH_ALTITUDE = """
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--altitude", type=int, required=True)
parser.add_argument("--results-out", required=True)
args, _ = parser.parse_known_args()

with open(args.results_out, "w") as results_file:
    json.dump({
        "Number of valid annotated segments": 10,
        "Mean altitude": args.altitude,
        "KDE reduction (%)": "5.0%",
        "1/3rd Octave Gain (F1)": 10.0,
        "F1": 0.85,
    }, results_file)
"""


class Test3dActiveSpaceWorkflow:
    """Verify the batch path used by generate_3d_active_space.py."""

    def test_multi_layer_commands_file_matches_3d_workflow(self, tmp_path: Path):
        package_root = tmp_path / "nps_active_space_pkg"
        scripts_dir = package_root / "scripts"
        scripts_dir.mkdir(parents=True)
        stub_path = scripts_dir / "generate_active_space.py"
        stub_path.write_text(STUB_SCRIPT_WITH_ALTITUDE)

        cmds_file = tmp_path / "DENATRLA2025_commands.txt"
        altitudes = [1800, 2100, 2400]
        cmds_file.write_text(
            "\n".join(
                f"DENATRLA2025_{alt}m\t-e GLBA_example -u DENA -s TRLA -y 2025 -l {alt}"
                for alt in altitudes
            ) + "\n"
        )

        output_df = pd.DataFrame()
        for line in cmds_file.read_text().splitlines():
            designator, options = re.split(r"\s+", line, maxsplit=1)
            cmd = [sys.executable, "-u", "-W", "ignore", str(stub_path)] + shlex.split(options)
            result_series = run_deployment(designator, cmd)
            assert result_series is not None
            output_df = pd.concat(
                [output_df, result_series.to_frame().T], ignore_index=True,
            )

        output_csv = tmp_path / "DENATRLA2025_commands_output.csv"
        output_df.to_csv(output_csv, index=False)

        loaded = pd.read_csv(output_csv)
        assert len(loaded) == len(altitudes)
        assert loaded["KDE reduction (%)"].notna().all()
        assert loaded["F1"].notna().all()
        assert set(loaded["Mean altitude"]) == set(altitudes)
        assert list(loaded["Designator"]) == [f"DENATRLA2025_{alt}m" for alt in altitudes]

    def test_generate_3d_active_space_invokes_batch_script(self):
        generate_3d_script = (
            Path(__file__).resolve().parents[2]
            / "nps_active_space"
            / "scripts"
            / "generate_3d_active_space.py"
        )
        source = generate_3d_script.read_text()
        assert "generate_active_space_batch.py" in source
        assert "fit_3d_active_space.py" in source
        assert '"-l", altitude' in source


class TestBatchMainOrchestration:
    def test_main_appends_csv_and_skips_existing_designator(self, tmp_path: Path, monkeypatch):
        cmds_file = tmp_path / "commands.txt"
        cmds_file.write_text(
            "RUN1\t-e test -u DENA -s TRLA -y 2025\n"
            "RUN2\t-e test -u DENA -s TRLA -y 2025\n"
        )
        output_csv = tmp_path / "output.csv"
        pd.DataFrame(
            {
                "Designator": ["RUN1"],
                "Number of valid annotated segments": [5],
                "Mean altitude": [3000],
                "KDE reduction (%)": ["12.5%"],
                "1/3rd Octave Gain (F1)": [12.5],
                "F1": [0.91],
            }
        ).to_csv(output_csv, index=False)

        run_calls: list[str] = []

        def fake_run_deployment(designator: str, cmd: list[str]) -> pd.Series:
            run_calls.append(designator)
            return pd.Series({
                "Designator": designator,
                "Number of valid annotated segments": 8,
                "Mean altitude": 3000,
                "KDE reduction (%)": "10.0%",
                "1/3rd Octave Gain (F1)": 11.0,
                "F1": 0.88,
            })

        monkeypatch.setattr(batch_module, "run_deployment", fake_run_deployment)
        monkeypatch.setattr(
            sys,
            "argv",
            ["generate_active_space_batch.py", str(cmds_file), "-o", str(output_csv)],
        )

        source = Path(batch_module.__file__).read_text()
        main_block = source[source.index("if __name__ == \"__main__\":"):]
        namespace = batch_module.__dict__.copy()
        namespace["__name__"] = "__main__"
        exec(compile(main_block, batch_module.__file__, "exec"), namespace)

        loaded = pd.read_csv(output_csv)
        assert run_calls == ["RUN2"]
        assert len(loaded) == 2
        assert set(loaded["Designator"]) == {"RUN1", "RUN2"}
        run2 = loaded.loc[loaded["Designator"] == "RUN2"].iloc[0]
        assert run2["Number of valid annotated segments"] == 8
        assert run2["F1"] == 0.88

