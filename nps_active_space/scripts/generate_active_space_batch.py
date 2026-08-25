import json
from nps_active_space import ACTIVE_SPACE_DIR
import subprocess
import os
import sys
import pandas as pd
import re
from argparse import ArgumentParser
import shlex
import tempfile
from pathlib import Path

from nps_active_space.active_space.active_space_setup import BATCH_RESULT_KEYS

RESULT_COLUMNS = ["Designator", *BATCH_RESULT_KEYS]
REQUIRED_RESULT_KEYS = list(BATCH_RESULT_KEYS)


def read_results_file(results_path: Path, designator: str) -> tuple[pd.Series | None, str | None]:
    """Load structured run results written by generate_active_space.py."""
    try:
        with open(results_path) as results_file:
            data = json.load(results_file)
    except json.JSONDecodeError as exc:
        return None, f"Invalid JSON in results file {results_path}: {exc}"

    if not isinstance(data, dict):
        return None, f"Results file must contain a JSON object: {results_path}"

    missing_keys = [key for key in REQUIRED_RESULT_KEYS if key not in data]
    if missing_keys:
        return None, (
            f"Results file is missing required keys {missing_keys}: {results_path}"
        )

    series = pd.Series(data)
    series["Designator"] = designator
    return series.reindex(RESULT_COLUMNS), None


def run_deployment(designator: str, cmd: list[str]) -> pd.Series | None:
    """Run generate_active_space.py and read structured results from a JSON sidecar."""
    env = os.environ.copy()
    env["TQDM_NCOLS"] = "80"

    fd, results_name = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    os.unlink(results_name)
    results_path = Path(results_name)

    try:
        full_cmd = cmd + ["--results-out", str(results_path)]
        process = subprocess.run(full_cmd, env=env)
        if process.returncode != 0:
            return None
        if not results_path.exists():
            print(
                f"Run succeeded but no results file was written: {results_path}",
                flush=True,
            )
            return None
        result_series, error_message = read_results_file(results_path, designator)
        if error_message is not None:
            print(error_message, flush=True)
            return None
        return result_series
    finally:
        if results_path.exists():
            results_path.unlink()


if __name__ == "__main__":
    argparse = ArgumentParser()
    argparse.add_argument("input", help="Path to input .txt file containing commands to run. This file is a sequence of lines,"
                                        "where each line starts with a unique designator identifying the run, followed by some"
                                        "whitespace, then followed by the options for the generate_active_space.py script."
                                        "An example line is this: DENATRLA2025  -e DENA_streamline -u DENA -s TRLA -y 2025 --cleanup")
    argparse.add_argument("-o", "--output", help="Path to output .csv file")
    args = argparse.parse_args()

    if args.output is None:
        name, ext = os.path.splitext(args.input)
        args.output = name + "_output.csv"
        print(f"No output file provided, using: {args.output}")

    assert args.input.endswith(".txt")
    assert args.output.endswith(".csv")

    output_df = pd.DataFrame()
    if os.path.exists(args.output):
        output_df = pd.read_csv(args.output)

    with open(args.input) as file:
        lines = [line.rstrip('\n') for line in file]

    for line in lines:
        print("")
        print(line)
        designator, options = re.split(r'\s+', line, maxsplit=1)

        if not output_df.empty and designator in output_df["Designator"].values:
            print("Done already, skipping")
            continue

        script_path = Path(ACTIVE_SPACE_DIR) / "scripts" / "generate_active_space.py"
        cmd = [sys.executable, "-u", "-W", "ignore", str(script_path)] + shlex.split(options)
        result_series = run_deployment(designator, cmd)
        if result_series is None:
            print(
                f"Run failed for {designator}; skipping CSV update. "
                "See errors above (check site inputs and propagation model outputs).",
                flush=True,
            )
            continue

        output_df = pd.concat(
            [output_df, result_series.to_frame().T], ignore_index=True)
        output_df.to_csv(args.output, index=False)
