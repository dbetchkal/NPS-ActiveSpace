"""Minimal stubs for batch subprocess tests."""

from __future__ import annotations

import sys
from pathlib import Path


def make_stub_generate_active_space_script(*, altitude_m: int | None = None) -> str:
    """Return stub script source that writes JSON results for ``run_deployment`` tests."""
    altitude_lines = ""
    mean_altitude = "3000"
    if altitude_m is not None:
        altitude_lines = 'parser.add_argument("-l", "--altitude", type=int, required=True)\n'
        mean_altitude = "args.altitude"

    return f"""
import argparse
import json

parser = argparse.ArgumentParser()
{altitude_lines}parser.add_argument("--results-out", required=True)
args, _ = parser.parse_known_args()

with open(args.results_out, "w") as results_file:
    json.dump({{
        "Number of valid annotated segments": 5,
        "Mean altitude": {mean_altitude},
        "KDE reduction (%)": "12.5%",
        "1/3rd Octave Gain (F1)": 12.5,
        "F1": 0.91,
    }}, results_file)
"""


def stub_generate_active_space_cmd(
    tmp_path: Path,
    *,
    altitude_m: int | None = None,
) -> list[str]:
    stub_path = tmp_path / "stub_generate_active_space.py"
    stub_path.write_text(make_stub_generate_active_space_script(altitude_m=altitude_m))
    return [sys.executable, str(stub_path)]
