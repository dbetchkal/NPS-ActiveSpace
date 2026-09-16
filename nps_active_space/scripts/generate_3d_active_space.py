import argparse
import os
import subprocess
import sys
import time
from datetime import UTC, datetime

# Headless CLI plots (savefig). Avoid TkAgg on Windows hosts.
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from nps_active_space.utils.computation import (
    compute_ambience_from_nvspl_archive,
    load_spectral_ambience_pickle,
)
import nps_active_space.utils.config as cfg
from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.scripts.generate_3d_commands import (
    build_layer_command_parts,
    format_commands_file_line,
)

"""
This script creates a commands file for use with generate_active_space_batch.py, containing commands
for each active space layer. It then uses generate_active_space_batch.py to run the commands 
(this can be disabled with an argument).

Before making a commands file, it also precomputes ambience from NVSPL. This saves a lot of time
since ambience doesn't have to be recomputed for each active space layer.
"""

ALTITUDE_STEP = 300  # meters between 3D active space layers


def _format_elapsed_s(elapsed_s: float) -> str:
    hours, rem = divmod(elapsed_s, 3600)
    minutes, seconds = divmod(rem, 60)
    if hours:
        return f"{int(hours)}h {int(minutes):02d}m {seconds:05.1f}s"
    if minutes:
        return f"{int(minutes)}m {seconds:05.1f}s"
    return f"{seconds:.1f}s"


def _log_pipeline_timing(
    label: str,
    start_ts: datetime,
    start_wall: float,
) -> None:
    end_ts = datetime.now(UTC)
    elapsed_s = time.perf_counter() - start_wall
    print(
        f"{label} finished at {end_ts.strftime('%Y-%m-%dT%H:%M:%SZ')} "
        f"(started {start_ts.strftime('%Y-%m-%dT%H:%M:%SZ')}, "
        f"elapsed {_format_elapsed_s(elapsed_s)})",
        flush=True,
    )


if __name__ == "__main__":
    pipeline_start_wall = time.perf_counter()
    pipeline_start_ts = datetime.now(UTC)
    print(
        f"3D active-space pipeline started at "
        f"{pipeline_start_ts.strftime('%Y-%m-%dT%H:%M:%SZ')}",
        flush=True,
    )

    parser = argparse.ArgumentParser()

    # arguments used only by this script
    parser.add_argument("--min-altitude", type=int, required=True, 
                          help=f"Minimum layer altitude (meters) for 3D active space. Should be a multiple of {ALTITUDE_STEP}.")
    parser.add_argument("--max-altitude", type=int, required=True, 
                          help=f"Maximum layer altitude (meters) for 3D active space. Should be a multiple of {ALTITUDE_STEP}.")
    parser.add_argument("--only-prep", action="store_true",
                          help="Stop after creating the commands file; don't run generate_active_space_batch.py. " \
                               "This is useful if you want to combine several command files into a single one to run at once, e.g. overnight.")

    # generate_active_space.py arguments that this script needs to know about
    parser.add_argument('-e', '--environment', required=True,
                          help="The configuration environment to run the script in.")
    parser.add_argument('-u', '--unit', required=True,
                          help="Four letter unit code. E.g. DENA")
    parser.add_argument('-s', '--site', required=True,
                          help="Four letter site code. E.g. TRLA")
    parser.add_argument('-y', '--year', type=int, required=True,
                          help="Four digit year. E.g. 2018")
    parser.add_argument('-a', '--ambience', default='nvspl',
                          help="What type of ambience to use in NMSIM calculations. Choose from ['nvspl', 'mennitt', or a path to an ambience .pkl file]")
    parser.add_argument('--model', type=AcousticModel, choices=list(AcousticModel),
                          default=AcousticModel.NMSIM,
                          help="Propagation model for each active-space layer command.")
    
    # other arguments will just be forwarded to generate_active_space.py via the batch commands text file
    args, extra_args = parser.parse_known_args()
    assert args.min_altitude <= args.max_altitude
    assert args.min_altitude % ALTITUDE_STEP == 0, f"Min altitude not a multiple of {ALTITUDE_STEP}"
    assert args.max_altitude % ALTITUDE_STEP == 0, f"Max altitude not a multiple of {ALTITUDE_STEP}"

    cfg.initialize(args.environment)
    project_dir = cfg.read("project", "dir")
    layout = p.SiteModelPaths.from_project(
        project_dir, args.unit, args.site, args.year, args.model,
    )
    site_dir = layout.site_dir
    usy = layout.usy

    # determine altitudes and print to console, so user can quickly verify we're doing what they wanted
    # before we get into NVSPL processing
    altitudes = list(range(args.min_altitude, args.max_altitude + 1, ALTITUDE_STEP))
    print(f"Using altitudes (m): {altitudes}")

    # Precompute NVSPL ambience to save time, if applicable
    if args.ambience == "nvspl":
        ambience_dir = layout.ambience_dir
        ambience_pkl_path = os.path.join(ambience_dir, f"{usy}_ambience.pkl")
        ambience_plot_path = os.path.join(ambience_dir, f"{usy}_ambience_plot.png")

        cached_ambience = load_spectral_ambience_pickle(ambience_pkl_path)
        if cached_ambience is not None:
            print(f"Found existing NVSPL ambience, using it: {p.display_path(ambience_pkl_path)}")
        else:
            if os.path.exists(ambience_pkl_path):
                print(
                    f"Existing ambience pickle has no usable spectral bands; recomputing: "
                    f"{p.display_path(ambience_pkl_path)}"
                )
            else:
                print("Computing NVSPL ambience")
            archive = cfg.read('data', 'nvspl_archive')
            ambience_quantile = 90  # L90 = 90% exceedance = 10% quantile sound level
            ambience = compute_ambience_from_nvspl_archive(
                archive,
                args.unit,
                args.site,
                args.year,
                ambience_quantile,
                broadband=False,
            )

            # make a plot too
            ambience.plot()
            plt.title(f"Ambient Spectrum: {usy}")
            plt.ylabel("Band SPL (dB)")
            plt.xlabel("1/3 Octave Band Frequency")
            plt.tight_layout()

            os.makedirs(ambience_dir, exist_ok=True)
            ambience.to_pickle(ambience_pkl_path)
            plt.savefig(ambience_plot_path)
            print(f"Saved ambience to {p.display_path(ambience_pkl_path)}")
    
    # create commands file
    cmds_file = os.path.join(site_dir, f"{usy}_commands.txt")
    with open(cmds_file, "w") as f:
        for altitude in altitudes:
            ambience_arg = ambience_pkl_path if args.ambience == "nvspl" else args.ambience
            parts = build_layer_command_parts(
                args.environment,
                args.unit,
                args.site,
                args.year,
                altitude,
                ambience_arg,
                model=args.model,
                extra_args=extra_args,
            )
            line = format_commands_file_line(f"{usy}_{altitude}m", parts)
            f.write(f"{line}\n")

    if args.only_prep:
        _log_pipeline_timing("3D active-space prep (commands file only)", pipeline_start_ts, pipeline_start_wall)
        sys.exit(0)

    print("Running generate_active_space_batch.py on the commands file\n")
    batch_script = os.path.join(os.path.dirname(__file__), "generate_active_space_batch.py")
    batch_process = subprocess.run(
        [sys.executable, batch_script, cmds_file],
        check=False,
    )
    if batch_process.returncode != 0:
        print(
            f"generate_active_space_batch.py exited with code {batch_process.returncode}. "
            "Skipping fit_3d_active_space.py. Fix batch errors above and rerun, "
            "or run the batch script directly on the commands file.",
            flush=True,
        )
        _log_pipeline_timing("3D active-space pipeline (batch failed)", pipeline_start_ts, pipeline_start_wall)
        sys.exit(batch_process.returncode)

    print("\nRunning fit_3d_active_space.py to fit the active space\n")
    fit_script = os.path.join(os.path.dirname(__file__), "fit_3d_active_space.py")
    fit_process = subprocess.run(
        [
            sys.executable,
            fit_script,
            "-e", args.environment,
            "-u", args.unit,
            "-s", args.site,
            "-y", str(args.year),
            "--model", args.model,
        ],
        check=False,
    )
    _log_pipeline_timing("3D active-space pipeline", pipeline_start_ts, pipeline_start_wall)
    sys.exit(fit_process.returncode)