import os
import pandas as pd
from nps_active_space.utils.helpers import load_layered_activespace, load_annotations
import nps_active_space.utils.config as cfg
from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.active_space.active_space_setup import upsert_project_fit
from argparse import ArgumentParser

"""
This script fits the gain for a 3D active space, assuming that active space layers
have already been generated using generate_active_space_batch.py or generate_active_space.py
"""

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-u", "--unit", required=True, help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", required=True, help="Four letter site code. E.g. TRLA")
    parser.add_argument("-y", "--year", required=True, help="Four digit year. E.g. 2018")
    parser.add_argument("-e", "--environment", required=True,
                        help="The configuration environment to run the script in.")
    parser.add_argument("--model", type=AcousticModel, choices=list(AcousticModel),
                        default=AcousticModel.NMSIM,
                        help="Propagation model whose active-space layers to fit.")
    args = parser.parse_args()

    unit, site, year = args.unit, args.site, args.year
    cfg.initialize(args.environment)
    project_dir = cfg.read("project", "dir")
    usy = f"{unit}{site}{year}"
    print(usy)

    layered = load_layered_activespace(project_dir, unit, site, year, model=args.model)
    annots = load_annotations(project_dir, unit, site, year)
    plot_savepath = p.precision_recall_3d_plot(
        project_dir, unit, site, year, model=args.model, beta=1.0,
    )
    os.makedirs(os.path.dirname(plot_savepath), exist_ok=True)
    result = layered.fit(annots, plot_savepath=plot_savepath)

    csv_file = upsert_project_fit(project_dir, usy, args.model, result)
    print(f"Fit results saved to {p.display_path(csv_file)}")
