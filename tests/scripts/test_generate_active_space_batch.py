import json
import shlex
from pathlib import Path

import pandas as pd
import pytest

import nps_active_space.utils.config as cfg
from nps_active_space.scripts.generate_active_space_batch import (
    batch_failure_hint,
    read_results_file,
    resolve_layer_output_dir,
    run_deployment,
    upsert_result_row,
)
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.utils.paths import layer_has_activespace_outputs
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

        result_series = run_deployment(designator, cmd, AcousticModel.NMSIM)

        assert result_series is not None
        assert result_series["Designator"] == designator
        assert result_series["Model"] == AcousticModel.NMSIM
        assert result_series["Mean altitude"] == altitude_m
        assert result_series["Number of valid annotated segments"] == 5


class TestLayerOutputSkip:
    def test_layer_has_activespace_outputs(self, tmp_path: Path) -> None:
        layer_dir = tmp_path / "DENATRLA2025_1500m"
        assert not layer_has_activespace_outputs(layer_dir)
        layer_dir.mkdir()
        assert not layer_has_activespace_outputs(layer_dir)
        (layer_dir / "DENATRLA2025_O_+000.geojson").write_text("{}")
        assert layer_has_activespace_outputs(layer_dir)

    def test_resolve_layer_output_dir_is_model_scoped(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        repo_root = Path(__file__).resolve().parents[2]
        monkeypatch.chdir(repo_root)
        cfg.initialize("GLBA_example", validate_config=False)

        nmsim_options = (
            "-e GLBA_example -u DENA -s TRLA -y 2025 -l 1500 --model nmsim"
        )
        aam_options = (
            "-e GLBA_example -u DENA -s TRLA -y 2025 -l 1500 --model aam"
        )

        nmsim_dir, nmsim_model = resolve_layer_output_dir(nmsim_options)
        aam_dir, aam_model = resolve_layer_output_dir(aam_options)

        assert nmsim_model is AcousticModel.NMSIM
        assert aam_model is AcousticModel.AAM
        assert nmsim_dir != aam_dir
        assert nmsim_dir.name == aam_dir.name == "DENATRLA2025_1500m"
        assert "Output_Data/nmsim/ACTIVESPACES" in nmsim_dir.as_posix()
        assert "Output_Data/aam/ACTIVESPACES" in aam_dir.as_posix()

    def test_batch_failure_hint_is_model_specific(self) -> None:
        site = "/data/DENATRLA"
        aam_hint = batch_failure_hint(site, AcousticModel.AAM)
        nmsim_hint = batch_failure_hint(site, AcousticModel.NMSIM)

        assert "Input_Data/aam" in aam_hint
        assert "TIG_TIS" not in aam_hint
        assert "03_TRAJECTORY" in nmsim_hint
        assert "nmsim/predictions" in nmsim_hint
        assert "TIG_TIS" not in nmsim_hint

    def test_upsert_result_row_replaces_same_model_only(self) -> None:
        df = pd.DataFrame([
            {"Designator": "DENATRLA2025_1500m", "Model": "nmsim", "F1": 0.5},
            {"Designator": "DENATRLA2025_1500m", "Model": "aam", "F1": 0.4},
        ])
        updated = upsert_result_row(
            df,
            pd.Series({
                "Designator": "DENATRLA2025_1500m",
                "Model": "aam",
                "F1": 0.9,
            }),
        )
        aam_rows = updated[updated["Model"] == "aam"]
        nmsim_rows = updated[updated["Model"] == "nmsim"]
        assert len(updated) == 2
        assert float(nmsim_rows.iloc[0]["F1"]) == 0.5
        assert float(aam_rows.iloc[0]["F1"]) == 0.9
