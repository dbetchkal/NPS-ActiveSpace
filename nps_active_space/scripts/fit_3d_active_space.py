import os
import pandas as pd
from nps_active_space.utils.helpers import load_layered_activespace, load_annotations
import nps_active_space.utils.config as cfg
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("unit", help="Four letter unit code. E.g. DENA")
    parser.add_argument("site", help="Four letter site code. E.g. TRLA")
    parser.add_argument("year", help="Four digit year. E.g. 2018")
    parser.add_argument("-e", "--environment", required=True,
                        help="The configuration environment to run the script in.")
    args = parser.parse_args()

    unit, site, year = args.unit, args.site, args.year
    cfg.initialize(args.environment)
    project_dir = cfg.read("project", "dir")
    print(unit+site+year)

    # fit active space model using annotations
    model = load_layered_activespace(project_dir, unit, site, year)
    annots = load_annotations(project_dir, unit, site, year)
    plot_savepath = Rf"{project_dir}/{unit}{site}/precision_recall_3d_f1_{unit}{site}{year}.png"
    result = model.fit(annots, plot_savepath=plot_savepath)
    
    # save to output csv
    csv_file = os.path.join(project_dir, "fits.csv")
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
        df = df[df["Designator"] != unit+site+year]
        df.loc[len(df)] = result
        df = df.sort_values(by="Designator")
    else:
        df = pd.DataFrame(result).T
    df.to_csv(csv_file, index=False)

    print(f"Fit results saved to {csv_file}")