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

from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel

RESULT_COLUMNS = [
    "Designator",
    "Model",
    "Number of valid annotated segments",
    "Mean altitude",
    "KDE reduction (%)",
    "1/3rd Octave Gain (F1)",
    "F1",
]
REQUIRED_RESULT_KEYS = [
    column for column in RESULT_COLUMNS if column not in {"Designator", "Model"}
]

_LAYER_OPTION_PARSER = ArgumentParser(add_help=False)
_LAYER_OPTION_PARSER.add_argument("-e", "--environment", required=True)
_LAYER_OPTION_PARSER.add_argument("-u", "--unit", required=True)
_LAYER_OPTION_PARSER.add_argument("-s", "--site", required=True)
_LAYER_OPTION_PARSER.add_argument("-y", "--year", type=int, required=True)
_LAYER_OPTION_PARSER.add_argument("-l", "--altitude", type=int, required=True)
_LAYER_OPTION_PARSER.add_argument(
    "--model",
    default=AcousticModel.NMSIM,
    type=AcousticModel,
    choices=list(AcousticModel),
)


def parse_layer_command_options(options: str):
    """Parse generate_active_space.py flags from one batch command line."""
    return _LAYER_OPTION_PARSER.parse_known_args(shlex.split(options))[0]


def resolve_layer_output_dir(options: str) -> tuple[Path, AcousticModel]:
    """Return the per-layer ACTIVESPACES directory for a batch command."""
    cmd_args = parse_layer_command_options(options)
    cfg.initialize(cmd_args.environment)
    site_dir = p.site_dir(cfg.read("project", "dir"), cmd_args.unit, cmd_args.site)
    model = AcousticModel.parse(cmd_args.model)
    layer_dir = Path(
        p.activespace_layer_dir(
            site_dir,
            cmd_args.unit,
            cmd_args.site,
            cmd_args.year,
            cmd_args.altitude,
            model,
        ),
    )
    return layer_dir, model


def layer_has_activespace_outputs(layer_dir: Path) -> bool:
    """True when the model-scoped layer folder already has omni geojson outputs."""
    if not layer_dir.is_dir():
        return False
    return any(layer_dir.glob("*_O_*.geojson"))


def batch_failure_hint(site_dir: str, model: AcousticModel) -> str:
    """Model-specific paths to inspect after a failed batch layer."""
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            return (
                f"check {site_dir}/Input_Data/aam, "
                f"{site_dir}/Output_Data/aam/runs, and active_space.log"
            )
        case AcousticModel.NMSIM:
            return (
                f"check {site_dir}/Input_Data/03_TRAJECTORY, "
                f"{site_dir}/Output_Data/nmsim/predictions, and "
                f"{site_dir}/Output_Data/nmsim/scratch"
            )


def upsert_result_row(output_df: pd.DataFrame, result_series: pd.Series) -> pd.DataFrame:
    """Replace any prior row for the same designator + model."""
    designator = result_series["Designator"]
    model = result_series["Model"]
    if not output_df.empty and "Designator" in output_df.columns:
        if "Model" in output_df.columns:
            mask = (output_df["Designator"] == designator) & (output_df["Model"] == model)
        else:
            mask = output_df["Designator"] == designator
        output_df = output_df[~mask]
    return pd.concat([output_df, result_series.to_frame().T], ignore_index=True)


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
    return series.reindex([c for c in RESULT_COLUMNS if c != "Model"]), None


def run_deployment(
    designator: str,
    cmd: list[str],
    model: AcousticModel,
) -> pd.Series | None:
    """
    Runs generate_active_space.py and reads structured results from a JSON output file.

    Parameters
    ----------
    designator: str
        Unique designator identifying this run. Included in the returned series.
    cmd: list
        Command to run, e.g. [sys.executable, "-u", "-W", "ignore", "...", "-e", "DENA_streamline", ...]
    model: AcousticModel
        Propagation model for this command (stored in the batch CSV).

    Returns
    -------
    results: pd.Series or None
        If the run errored out, returns None.
        If the run succeeded, returns a pandas series with the results. Contains columns:
        "Designator", "Model", "Number of valid annotated segments", "Mean altitude",
        "KDE reduction (%)", "1/3rd Octave Gain (F1)", "F1"
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
        result_series["Model"] = AcousticModel.parse(model)
        return result_series.reindex(RESULT_COLUMNS)
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

    output_csv = Path(args.output).resolve()

    # initialize output dataframe / load existing results
    output_df = pd.DataFrame()
    if output_csv.is_file():
        output_df = pd.read_csv(output_csv)

    # read generate active space commands from input file
    with open(args.input) as file:
        lines = [line.rstrip('\n') for line in file]

    # run each generate active space command
    for line in lines:
        print("")
        print(line)
        designator, options = re.split(r'\s+', line, maxsplit=1)

        layer_dir, model = resolve_layer_output_dir(options)
        if layer_has_activespace_outputs(layer_dir):
            print(
                f"Skipping {designator} ({model}): active-space outputs already in "
                f"{layer_dir} (delete that directory to force rerun)."
            )
            continue

        # assemble and run the command
        script_path = Path(ACTIVE_SPACE_DIR) / "scripts" / "generate_active_space.py"
        cmd = [sys.executable, "-u", "-W", "ignore", str(script_path)] + shlex.split(options)
        cmd_args = parse_layer_command_options(options)
        site_dir = p.site_dir(cfg.read("project", "dir"), cmd_args.unit, cmd_args.site)
        result_series = run_deployment(designator, cmd, model)
        if result_series is None:
            print(
                f"Run failed for {designator} ({model}); skipping CSV update. "
                f"See errors above ({batch_failure_hint(site_dir, model)})."
            )
            continue

        output_df = upsert_result_row(output_df, result_series)
        output_df.to_csv(output_csv, index=False)

        # if args.savedir is not None:
        #     copy_output_files(options, args.savedir, designator)
