import argparse
import subprocess
import iyore
import os
import matplotlib.pyplot as plt
from nps_active_space.utils.models import Nvspl
from nps_active_space.utils.computation import ambience_from_nvspl
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


if __name__ == "__main__":
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
    site_dir = p.site_dir(project_dir, args.unit, args.site)
    usy = f"{args.unit}{args.site}{args.year}"

    # determine altitudes and print to console, so user can quickly verify we're doing what they wanted
    # before we get into NVSPL processing
    altitudes = list(range(args.min_altitude, args.max_altitude + 1, ALTITUDE_STEP))
    print(f"Using altitudes (m): {altitudes}")

    # Precompute NVSPL ambience to save time, if applicable
    if args.ambience == "nvspl":
        ambience_dir = os.path.join(site_dir, "Output_Data", "AMBIENCE")
        ambience_pkl_path = os.path.join(ambience_dir, f"{usy}_ambience.pkl")
        ambience_plot_path = os.path.join(ambience_dir, f"{usy}_ambience_plot.png")

        if os.path.exists(ambience_pkl_path):
            print(f"Found existing NVSPL ambience, using it: {ambience_pkl_path}")
        else:
            print("Computing NVSPL ambience")
            archive = iyore.Dataset(cfg.read('data', 'nvspl_archive'))
            nvspl_files = [e.path for e in archive.nvspl(unit=args.unit, site=args.site, year=str(args.year))]
            nvspl = Nvspl(nvspl_files)
            ambience_quantile = 90  # L90 = 90% exceedance = 10% quantile sound level
            ambience = ambience_from_nvspl(nvspl, ambience_quantile, broadband=False)

            # make a plot too
            ambience.plot()
            plt.title(f"Ambient Spectrum: {usy}")
            plt.ylabel("Band SPL (dB)")
            plt.xlabel("1/3 Octave Band Frequency")
            plt.tight_layout()

            os.makedirs(ambience_dir, exist_ok=True)
            ambience.to_pickle(ambience_pkl_path)
            plt.savefig(ambience_plot_path)
            print(f"Saved ambience to {ambience_pkl_path}")
    
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
    
    # if we're not only prepping a commands file, run generate_active_space_batch.py
    # and then fit the active space
    if not args.only_prep:
        print("Running generate_active_space_batch.py on the commands file\n")
        batch_script = os.path.join(os.path.dirname(__file__), "generate_active_space_batch.py")
        process = subprocess.Popen(
            ["python", batch_script, cmds_file]
        )
        process.wait()

        print("\nRunning fit_3d_active_space.py to fit the active space\n")
        fit_script = os.path.join(os.path.dirname(__file__), "fit_3d_active_space.py")
        process = subprocess.Popen(
            ["python", fit_script, "-e", args.environment, "-u", args.unit, "-s", args.site,
             "-y", str(args.year), "--model", args.model]
        )
        process.wait()