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
sys.path.append(".")
import _DENA.resource.config as cfg
from _DENA import DENA_DIR


def stream_and_capture(stream, buffer, target):
    for chunk in iter(lambda: stream.read(1024), b''):
        buffer.append(chunk)
        target.write(chunk)
        target.flush()
    stream.close()


def run_deployment(designator, cmd):
    # configure tqdm
    env = os.environ.copy()
    env["TQDM_NCOLS"] = "80"

    process = subprocess.Popen(
        shlex.split(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0,  # unbuffered
        env=env
    )

    stdout_chunks = []
    stderr_chunks = []

    stdout_thread = threading.Thread(
        target=stream_and_capture, args=(process.stdout, stdout_chunks, sys.stdout.buffer)
    )
    stderr_thread = threading.Thread(
        target=stream_and_capture, args=(process.stderr, stderr_chunks, sys.stderr.buffer)
    )

    stdout_thread.start()
    stderr_thread.start()

    process.wait()
    stdout_thread.join()
    stderr_thread.join()

    if process.returncode != 0:
        return None

    stdout_text = b''.join(stdout_chunks).decode(errors='replace')
    return parse_output(stdout_text, designator)


def parse_output(s, designator):
    altitude_match = re.search(r"Average altitude is: (\d+)m", s)
    avg_altitude = float(altitude_match.group(1)) if altitude_match is not None else ""
    kde_match = re.search(r"Went from (\d+) to (\d+) points when normalizing point density", s)
    points_before_kde = int(kde_match.group(1))
    points_after_kde = int(kde_match.group(2))
    kde_reduction = f"{100 * (1 - (points_after_kde / points_before_kde))}%"
    n_annots = int(re.search(r"(\d+) valid annotated segments found", s).group(1))
    best_f1_match = re.search(r"The best performing omni source for F-1.0 is: O_(....) \(fbeta: (.+)\)", s)
    gain = int(best_f1_match.group(1)) / 10
    f1 = float(best_f1_match.group(2))

    return pd.Series({
        "Designator": designator,
        "Number of valid annotated segments": n_annots,
        "Mean altitude": avg_altitude,
        "KDE reduction (%)": kde_reduction,
        "1/3rd Octave Gain (F1)": gain,
        "F1": f1
    })


def copy_output_files(option_str, savedir, designator):
    dst_dir = os.path.join(savedir, designator)
    os.makedirs(dst_dir, exist_ok=True)
    print(f"Copying output to {dst_dir}...")

    # parse command options to get project dir, unit, and site
    argparse = ArgumentParser()
    argparse.add_argument("-e", "--environment")
    argparse.add_argument("-u", "--unit")
    argparse.add_argument("-s", "--site")
    argparse.add_argument("-y", "--year")
    argparse.add_argument('--annotation-file')
    args, _ = argparse.parse_known_args(shlex.split(option_str))
    cfg.initialize(f"{DENA_DIR}/config", environment=args.environment)
    site_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"

    deployment = f"{args.unit}{args.site}{args.year}"
    activespace_files = glob.glob(os.path.join(site_dir, f"{deployment}_O_*.geojson"))
    pr_plots = glob.glob(os.path.join(site_dir, f"PrecisionRecallPlot_{deployment}*.png"))
    if args.annotation_file is not None:
        annotation_files = [os.path.join(site_dir, args.annotation_file)]
    else:
        annotation_files = glob.glob(os.path.join(site_dir, f"{args.unit}{args.site}{args.year}*saved_annotations*.geojson"))
    
    for src_path in activespace_files + pr_plots + annotation_files:
        basename = os.path.basename(src_path)
        dst_path = os.path.join(dst_dir, basename)
        shutil.copy2(src_path, dst_path)  # copy2 to preserve metadata


if __name__ == "__main__":
    argparse = ArgumentParser()
    argparse.add_argument("input", help="Path to input .txt file containing commands to run")
    argparse.add_argument("-o", "--output", help="Path to output .csv file")
    argparse.add_argument("-s", "--savedir", help="Parent directory to copy the output files to. A subdirectory named with the designator will contain the files.")
    args = argparse.parse_args()

    if args.output is None:
        name, ext = os.path.splitext(args.input)
        args.output = name + "_output.csv"
        print(f"No output file provided, using: {args.output}")

    assert args.input.endswith(".txt")
    assert args.output.endswith(".csv")
    assert os.getcwd().endswith("NPS-ActiveSpace"), "Script needs to be run from root of NPS-ActiveSpace repository"

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

        cmd = Rf"python -u -W ignore '_DENA\scripts\generate_active_space.py' {options}"
        result_series = run_deployment(designator, cmd)
        if result_series is not None:
            output_df = pd.concat([output_df, result_series.to_frame().T], ignore_index=True)
            output_df.to_csv(args.output, index=False)
        
        if args.savedir is not None:
            copy_output_files(options, args.savedir, designator)
        