import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import glob
from nps_active_space.utils.helpers import load_annotations
import nps_active_space.utils.config as cfg
from argparse import ArgumentParser
import re

"""
This module plots the distribution of annotated segment altitudes. It is intended to be used
before generating a 3D active space, to help the user determine which altitudes should be used.
"""

def plot_altitude_hist(project_dir, unit, site, year):
    """Plots three altitude histograms for a site using data from an annotation file.
    Histograms are - track point altitudes, mean track altitudes, median track altitudes.
    
    Parameters
    ----------
    project_dir: str
        Path to directory containing site directories (e.g. DENATRLA is a site directory)
    unit: str
    site: str
    year: str
    """
    print(f"Plotting altitudes for {unit}{site}{year}")

    # load annotations
    annotations = load_annotations(project_dir, unit, site, year)

    print("Extracting altitudes")
    # create a dataframe holding the z values
    rows = []
    for _, annot in annotations[annotations.audible].iterrows():
        for xyz in annot.geometry.coords:
            z = xyz[2]
            if z < 0 or z > 10000:
                continue
            rows.append({"_id": annot["_id"], "z": z})
    altitudes = pd.DataFrame(rows)
    mean_altitudes = altitudes.groupby("_id").mean()

    # calculate the altitude we used to fit the active space (mean of track means)
    z_active_space = mean_altitudes["z"].mean()

    # consider using the median of median track altitude
    median_altitudes = altitudes.groupby("_id").median()
    median_z = median_altitudes["z"].median()
    
    print("Plotting")

    def add_vlines(axis):
        axis.axvline(x=z_active_space, color="red", label=f"2D Active Space Altitude ({int(z_active_space)}m)\n= Mean of Mean Track Altitudes")
        axis.axvline(x=median_z, color="black", label=f"Median of Median Track Altitudes ({int(median_z)}m)")

    fig, ax = plt.subplots(nrows=3, figsize=(8, 12), sharex=True)

    ax[0].set_title(f"{unit}{site}{year} Track Point Altitudes")
    ax[0].hist(altitudes["z"], bins=50)
    add_vlines(ax[0])
    ax[0].set_xlabel("Altitude MSL (m)")
    ax[0].set_ylabel("Number of Track Points")
    ax[0].legend()

    ax[1].set_title(f"{unit}{site}{year} Mean Track Altitudes")
    ax[1].hist(mean_altitudes["z"], bins=20)
    add_vlines(ax[1])
    ax[1].set_xlabel("Altitude MSL (m)")
    ax[1].set_ylabel("Number of Tracks")
    ax[1].legend()

    ax[2].set_title(f"{unit}{site}{year} Median Track Altitudes")
    ax[2].hist(median_altitudes["z"], bins=20)
    add_vlines(ax[2])
    ax[2].set_xlabel("Altitude MSL (m)")
    ax[2].set_ylabel("Number of Tracks")
    ax[2].legend()

    # make ticks fall along the 300m intervals of a 3D active space
    # reenable tick labels for all axes, makes it easier to read in practice
    for axis in ax:
        span = altitudes["z"].max() - altitudes["z"].min()
        interval = 300 if span < 4500 else 600
        axis.xaxis.set_major_locator(MultipleLocator(interval))  # ticks w labels
        axis.xaxis.set_minor_locator(MultipleLocator(300))  # always mark 300m intervals, even if unlabeled
        axis.tick_params(labelbottom=True)

    fig.tight_layout()

    fig.savefig(f"{project_dir}/{unit}{site}/Altitude_Histogram_{unit}{site}{year}.png")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-e", "--environment", required=True,
                        help="The configuration environment to run the script in.")
    parser.add_argument("-u", "--unit", help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", help="Four letter site code. E.g. TRLA")
    parser.add_argument("-y", "--year",  help="Four digit year. E.g. 2018")
    parser.add_argument("-a", "--all", action="store_true",
                        help="If provided, plot altitude histograms for all sites in the project directory." \
                             "Should not be passed if -u -s -y flags are passed.")
    args = parser.parse_args()
    cfg.initialize(args.environment)
    project_dir = cfg.read("project", "dir")

    if not args.all:
        for arg in [args.unit, args.site, args.year]:
            assert arg is not None, "Must pass --unit, --site, and --year arguments"

        plot_altitude_hist(project_dir, args.unit, args.site, args.year)

    else:
        for arg in [args.unit, args.site, args.year]:
            assert arg is None, "Do not pass --unit, --site, or --year arguments if passing --all"

        print("Plotting altitudes for all sites\n")

        # use glob to find all annotation files
        annotation_files = glob.glob(f"{project_dir}/*/*saved_annotations*.geojson")
        for file in annotation_files:
            # get unit, site, year from path name
            site_name = os.path.basename(os.path.dirname(file))
            unit = site_name[:4]
            site = site_name[4:]
            year_match = re.search(rf"{unit}{site}(\d\d\d\d)", file)
            if year_match is None:
                print("Failed to find the year")
                continue
            year = int(year_match.group(1))

            plot_altitude_hist(project_dir, unit, site, year)
            print("")
