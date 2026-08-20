import json
import shlex
import sys
from pathlib import Path

import pandas as pd
import pytest

import nps_active_space.scripts.generate_active_space_batch as batch_module
from nps_active_space.scripts.generate_active_space_batch import (
    read_results_file,
    run_deployment,
)
from script_test_helpers import exec_script_main, stub_generate_active_space_cmd

VALID_RESULTS = {
    "Number of valid annotated segments": 12,
    "Mean altitude": 3000,
    "KDE reduction (%)": "15.0%",
    "1/3rd Octave Gain (F1)": 12.5,
    "F1": 0.87,
}
DESIGNATOR = "DENATRLA20253000m"


class TestReadResultsFile:
    @pytest.mark.parametrize(
        ("file_content", "expect_error_substring", "expected_values"),
        [
            (
                json.dumps(VALID_RESULTS),
                None,
                VALID_RESULTS,
            ),
            (
                json.dumps({
                    "Number of valid annotated segments": 12,
                    "Mean altitude": 3000,
                    "KDE reduction (%)": "15.0%",
                    "1/3rd Octave Gain (F1)": 12.5,
                }),
                "missing required keys",
                None,
            ),
            (
                "{not valid json",
                "Invalid JSON",
                None,
            ),
        ],
    )
    def test_read_results_file(
        self,
        tmp_path: Path,
        file_content: str,
        expect_error_substring: str | None,
        expected_values: dict[str, object] | None,
    ):
        results_path = tmp_path / "results.json"
        results_path.write_text(file_content)

        series, error_message = read_results_file(results_path, DESIGNATOR)

        if expect_error_substring is None:
            assert error_message is None
            expected = pd.Series({"Designator": DESIGNATOR, **expected_values})
            pd.testing.assert_series_equal(series, expected)
        else:
            assert series is None
            assert error_message is not None
            assert expect_error_substring in error_message


class TestRunDeployment:
    def test_run_deployment_reads_stub_subprocess_results(self, tmp_path: Path):
        altitude_m = 2100
        designator = f"DENATRLA2025_{altitude_m}m"
        options = f"-e GLBA_example -u DENA -s TRLA -y 2025 -l {altitude_m}"
        stub_cmd = stub_generate_active_space_cmd(tmp_path, altitude_m=altitude_m)
        cmd = [stub_cmd[0], "-u", "-W", "ignore", stub_cmd[1]] + shlex.split(options)

        result_series = run_deployment(designator, cmd)

        assert result_series is not None
        assert result_series["Designator"] == designator
        assert result_series["Mean altitude"] == altitude_m
        assert result_series["Number of valid annotated segments"] == 5


class TestBatchMainOrchestration:
    def test_main_skips_existing_appends_success_and_ignores_failed_run(
        self,
        tmp_path: Path,
        monkeypatch,
    ):
        cmds_file = tmp_path / "commands.txt"
        cmds_file.write_text(
            "RUN1\t-e test -u DENA -s TRLA -y 2025\n"
            "RUN2\t-e test -u DENA -s TRLA -y 2025\n"
            "RUN3\t-e test -u DENA -s TRLA -y 2025\n"
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

        def fake_run_deployment(designator: str, cmd: list[str]) -> pd.Series | None:
            run_calls.append(designator)
            if designator == "RUN3":
                return None
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

        exec_script_main(batch_module)

        loaded = pd.read_csv(output_csv)
        assert run_calls == ["RUN2", "RUN3"]
        assert len(loaded) == 2
        assert set(loaded["Designator"]) == {"RUN1", "RUN2"}
        run2 = loaded.loc[loaded["Designator"] == "RUN2"].iloc[0]
        assert run2["Number of valid annotated segments"] == 8
        assert run2["F1"] == 0.88
