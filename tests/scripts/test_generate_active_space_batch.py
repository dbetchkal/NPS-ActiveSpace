import json
import shlex
from pathlib import Path

import pandas as pd
import pytest

from nps_active_space.scripts.generate_active_space_batch import (
    read_results_file,
    run_deployment,
)
from script_test_helpers import stub_generate_active_space_cmd

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
            (json.dumps(VALID_RESULTS), None, VALID_RESULTS),
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
            ("{not valid json", "Invalid JSON", None),
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
        assert result_series["F1"] == 0.91
