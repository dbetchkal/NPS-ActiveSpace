import json
from pathlib import Path

import pandas as pd

from nps_active_space.scripts.generate_active_space_batch import read_results_file


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
