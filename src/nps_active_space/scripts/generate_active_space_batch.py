import nps_active_space.utils.config as cfg
from nps_active_space import ACTIVE_SPACE_DIR
import subprocess
import os
import sys
import pandas as pd
import re
from argparse import ArgumentParser
import threading
import glob
import shlex
import shutil
from pathlib import Path
from typing import IO


def stream_and_capture(stream: IO[bytes], buffer: list[bytes], target: IO[bytes]) -> None:
    """
    Utility function for capturing stdout or stderr, and then forwarding it to the console.

    Parameters
    ----------
    stream
        Stream to capture. Should be process.stdout or process.stderr
    buffer: list
        A list to store the bytes in, for reading later.
    target
        Where to forward the stream to. Should be sys.stdout.buffer or sys.stderr.buffer
    """
    for chunk in iter(lambda: stream.read(1024), b''):
        buffer.append(chunk)
        target.write(chunk)
        target.flush()
    stream.close()


def run_deployment(designator: str, cmd: list[str]) -> pd.Series | None:
    """
    Runs generate active space, and parses the printed output to extract results we want.

    Parameters
    ----------
    designator: str
        Unique designator identifying this run. Included in the returned series.
    cmd: list
        Command line command to run as a list, e.g. ["python", "-u", "-W", "ignore", "nps_active_space/scripts/generate_active_space.py", "-e", "DENA_streamline", ...]

    Returns
    -------
    results: pd.Series or None
        If the run errored out, returns None.
        If the run succeeded, returns a pandas series with the results. Contains columns:
        "Designator", "Number of valid annotated segments", "Mean altitude", "KDE reduction (%)", "1/3rd Octave Gain (F1)", "F1"
    """

    # Configure tqdm - environment variables preceded by "TQDM_" get passed to all tqdm instances as default args.
    # Since we are capturing the tqdm stream and then re-printing, tqdm doesn't know anything about the width
    # of our display, and produces a very small progress bar as a result. So, make the progress bar a bit wider.
    env = os.environ.copy()
    env["TQDM_NCOLS"] = "80"  # 80 characters wide progress bar

    # Run the command
    process = subprocess.Popen(
        cmd,
        # capture printed output instead of printing to console
        stdout=subprocess.PIPE,
        # capture stderr output (e.g. tqdm) instead of printing to console
        stderr=subprocess.PIPE,
        bufsize=0,  # unbuffered
        env=env
    )

    # create lists to store byte chunks that are printed to the console
    stdout_chunks = []
    stderr_chunks = []

    # In order to print stdout and stderr to the console in the same way as they
    # would originally be printed, we need to capture them at the same time.
    #   - BTW, not capturing stderr and letting it print naturally will mess up the ordering
    #     of the stdout vs. stderr prints to the console (e.g. progress bars may come before
    #     the associated print statements)
    # To capture stdout and stderr at the same time, we need to iterate
    # through both their streams simultaneously. This requires doing two things
    # at the same time, so we use threading to allow this parallel operation.

    # Initialize the threads.
    stdout_thread = threading.Thread(
        target=stream_and_capture, args=(
            process.stdout, stdout_chunks, sys.stdout.buffer)
    )
    stderr_thread = threading.Thread(
        target=stream_and_capture, args=(
            process.stderr, stderr_chunks, sys.stderr.buffer)
    )
    # Run the threads.
    stdout_thread.start()
    stderr_thread.start()

    # Wait for the threads to finish.
    process.wait()
    stdout_thread.join()
    stderr_thread.join()

    # If we errored out, return nothing
    if process.returncode != 0:
        return None

    # Convert byte output to text output
    stdout_text = b''.join(stdout_chunks).decode(errors='replace')

    # Parse the printed output to get the results series
    return parse_output(stdout_text, designator)


def parse_output(s: str, designator: str) -> pd.Series:
    """
    Parses printed output to get the relevant information.

    Parameters
    ----------
    s: str
        String containing all printed output.
    designator: str
        Unique designator for this run. Included in the returned series.

    Returns
    -------
    results: pd.Series
        Pandas series with columns:
        "Designator", "Number of valid annotated segments", "Mean altitude", "KDE reduction (%)", "1/3rd Octave Gain (F1)", "F1"
    """
    # Use regular expressions to find the information we want.
    # Use capturing groups (parentheses in the regex) to extract just the part
    # we are interested in, and access these with the .group(index) function.

    altitude_match = re.search(r"Average altitude is: (\d+)m", s)
    avg_altitude = float(altitude_match.group(1)) \
        if altitude_match is not None else ""

    kde_match = re.search(
        r"Went from (\d+) to (\d+) points when normalizing point density", s)
    points_before_kde = int(kde_match.group(1))
    points_after_kde = int(kde_match.group(2))
    kde_reduction = f"{100 * (1 - (points_after_kde / points_before_kde))}%"

    n_annots = int(
        re.search(r"(\d+) valid annotated segments found", s).group(1))

    best_f1_match = re.search(
        r"The best performing omni source for F-1.0 is: O_(....) \(fbeta: (.+)\)", s)
    if best_f1_match is not None:
        gain = int(best_f1_match.group(1)) / 10
        f1 = float(best_f1_match.group(2))
    else:
        gain = None
        f1 = None

    return pd.Series({
        "Designator": designator,
        "Number of valid annotated segments": n_annots,
        "Mean altitude": avg_altitude,
        "KDE reduction (%)": kde_reduction,
        "1/3rd Octave Gain (F1)": gain,
        "F1": f1
    })


def copy_output_files(option_str: str, savedir: str, designator: str) -> None:
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

    config = cfg.load_config(args.environment)
    site_dir = f"{config.project.dir}/{args.unit}{args.site}"
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
        # if it ran with no errors, save the results
        if result_series is not None:
            output_df = pd.concat(
                [output_df, result_series.to_frame().T], ignore_index=True)
            output_df.to_csv(args.output, index=False)

        # if args.savedir is not None:
        #     copy_output_files(options, args.savedir, designator)
