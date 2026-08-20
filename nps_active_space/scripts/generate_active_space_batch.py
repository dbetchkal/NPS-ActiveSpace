import json
import nps_active_space.utils.config as cfg
from nps_active_space import ACTIVE_SPACE_DIR
import subprocess
import os
import sys
import pandas as pd
import re
from argparse import ArgumentParser
import glob
import shlex
import shutil
import tempfile
from pathlib import Path

RESULT_COLUMNS = [
    "Designator",
    "Number of valid annotated segments",
    "Mean altitude",
    "KDE reduction (%)",
    "1/3rd Octave Gain (F1)",
    "F1",
]
REQUIRED_RESULT_KEYS = [column for column in RESULT_COLUMNS if column != "Designator"]


def read_results_file(results_path: Path, designator: str) -> tuple[pd.Series | None, str | None]:
    """Load structured run results written by generate_active_space.py.

    Returns
    -------
    series, error_message
        ``error_message`` is set when results cannot be loaded (for console display).
    """
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
    """
    Runs generate_active_space.py and reads structured results from a JSON output file.

    Parameters
    ----------
    designator: str
        Unique designator identifying this run. Included in the returned series.
    cmd: list
        Command to run, e.g. [sys.executable, "-u", "-W", "ignore", "...", "-e", "DENA_streamline", ...]

    Returns
    -------
    results: pd.Series or None
        If the run errored out, returns None.
        If the run succeeded, returns a pandas series with the results. Contains columns:
        "Designator", "Number of valid annotated segments", "Mean altitude", "KDE reduction (%)",
        "1/3rd Octave Gain (F1)", "F1"
    """
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


def copy_output_files(option_str, savedir, designator):
    """
    Copies output files from the site directory to another directory.
    This keeps things organized, and avoids these files being overwritten
    by future runs that use the same site directory.

    Parameters
    ----------
    option_str: str
        The string representing the command options for this run. E.g.,
        "-e DENA_streamline -u DENA -s TRLA -y 2025 --cleanup"
    savedir: str
        Path to a directory to save files to. Files will be copied to a subdirectory
        in savedir named with the designator.
    designator: str
        Unique string representing this run.
    """
    dst_dir = os.path.join(savedir, designator)
    os.makedirs(dst_dir, exist_ok=True)
    print(f"Copying output to {dst_dir}...")

    # Parse command options to get project dir, unit, and site.
    # We need these to locate the site directory where the files we want to copy are,
    # and to figure out which files are relevant to this run that we should copy
    argparse = ArgumentParser()
    argparse.add_argument("-e", "--environment")
    argparse.add_argument("-u", "--unit")
    argparse.add_argument("-s", "--site")
    argparse.add_argument("-y", "--year")
    argparse.add_argument('--annotation-file')
    # Use shlex to convert option str into a list for argparse.
    # shlex avoids splitting on spaces that are inside quotes
    args, _ = argparse.parse_known_args(shlex.split(option_str))

    cfg.initialize(environment=args.environment)
    site_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    deployment = f"{args.unit}{args.site}{args.year}"

    # Get filenames we wish to copy
    activespace_files = glob.glob(os.path.join(
        site_dir, "Output_Data", "ACTIVESPACES", f"{deployment}_O_*.geojson"))
    tested_pt_files = glob.glob(os.path.join(
        site_dir, "Output_Data", "TESTED_POINTS", f"{deployment}_O_*.pkl"))
    pr_plots = glob.glob(os.path.join(
        site_dir, f"PrecisionRecallPlot_{deployment}*.png"))
    # If a custom annotation file was used, copy that.
    # Otherwise, copy the default annotations file(s)
    if args.annotation_file is not None:
        annotation_files = [os.path.join(site_dir, args.annotation_file)]
    else:
        annotation_files = glob.glob(os.path.join(
            site_dir, f"{args.unit}{args.site}{args.year}*saved_annotations*.geojson"))

    # Copy files
    for src_path in activespace_files + tested_pt_files + pr_plots + annotation_files:
        basename = os.path.basename(src_path)
        dst_path = os.path.join(dst_dir, basename)
        shutil.copy2(src_path, dst_path)  # copy2 to preserve metadata


if __name__ == "__main__":
    argparse = ArgumentParser()
    argparse.add_argument("input", help="Path to input .txt file containing commands to run. This file is a sequence of lines,"
                                        "where each line starts with a unique designator identifying the run, followed by some"
                                        "whitespace, then followed by the options for the generate_active_space.py script."
                                        "An example line is this: DENATRLA2025  -e DENA_streamline -u DENA -s TRLA -y 2025 --cleanup")
    argparse.add_argument("-o", "--output", help="Path to output .csv file")
    # argparse.add_argument("-s", "--savedir", help="Parent directory to copy the output files to. Output files are the active spaces,"
    #                                               "the annotations used, and the precision-recall plot. A subdirectory named with"
    #                                               "the designator will contain the files.")
    args = argparse.parse_args()

    if args.output is None:
        name, ext = os.path.splitext(args.input)
        args.output = name + "_output.csv"
        print(f"No output file provided, using: {args.output}")

    assert args.input.endswith(".txt")
    assert args.output.endswith(".csv")

    # initialize output dataframe / load existing results
    output_df = pd.DataFrame()
    if os.path.exists(args.output):
        output_df = pd.read_csv(args.output)

    # read generate active space commands from input file
    with open(args.input) as file:
        lines = [line.rstrip('\n') for line in file]

    # run each generate active space command
    for line in lines:
        print("")
        print(line)
        designator, options = re.split(r'\s+', line, maxsplit=1)

        # check if this run has been done already, based on if the output file has a matching designator,
        # since designators should be unique to a run
        if not output_df.empty and designator in output_df["Designator"].values:
            print("Done already, skipping")
            continue

        # assemble and run the command
        script_path = Path(ACTIVE_SPACE_DIR) / "scripts" / "generate_active_space.py"
        cmd = [sys.executable, "-u", "-W", "ignore", str(script_path)] + shlex.split(options)
        result_series = run_deployment(designator, cmd)
        if result_series is None:
            print(
                f"Run failed for {designator}; skipping CSV update. "
                "See errors above (check site Input_Data/03_TRAJECTORY and Output_Data/TIG_TIS)."
            )
            continue

        output_df = pd.concat(
            [output_df, result_series.to_frame().T], ignore_index=True)
        output_df.to_csv(args.output, index=False)

        # if args.savedir is not None:
        #     copy_output_files(options, args.savedir, designator)
