import logging
import geopandas as gpd
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import contextlib
import os
from shapely.geometry import Point
from nps_active_space.utils.helpers import load_annotations, load_studyarea, load_layered_activespace
from nps_active_space.utils.computation import normalize_point_density, NMSIM_bbox_utm

logger = logging.getLogger(__name__)

__all__ = [
    "get_csv_name",
    "fit_varying_n_tracks",
    "plot_stability"
]

def get_csv_name(project_dir, deployment):
    unit, site = deployment[:4], deployment[4:-4]
    return os.path.join(project_dir, unit+site, "Output_Data", "VALIDATION",
                        f"{deployment}_study_duration_stability.csv")


def fit_varying_n_tracks(project_dir, unit, site, year):
    """
    Fits a 3D active space with 1, 2, 3, ..., all tracks.
    Saves a CSV containing, for each number of tracks:
    - the best f1/gain,
    - the area of that fit and the altitude used for area calculation,
    - the number of track points,
    - the number of valid annotated segments,
    - rankings of the 1st, 2nd, 3rd best, etc. gains/f1s/areas

    Parameters
    ----------
    project_dir: str
        Path to the project directory, containing site folders.
    unit: str
        Four letter unit code. E.g. DENA
    site: str
        Site code. E.g. TRLA
    year: str
        Year as a string. E.g. "2025".
    """
    deployment = f"{unit}{site}{year}"
    logger.info(f"Processing deployment: {deployment}")

    savefile = get_csv_name(project_dir, unit+site+year)
    if os.path.exists(savefile):
        logger.info(f"Skipping, found existing results in {savefile}")
        return

    study_area = load_studyarea(project_dir, unit, site, year)
    utm_zone = NMSIM_bbox_utm(study_area)
    annots = load_annotations(project_dir, unit, site, year).to_crs(utm_zone)

    # Determine chronological order of tracks
    annots.sort_values(by="start_dt", inplace=True)
    track_order = annots.drop_duplicates(subset="_id", keep="first")["_id"]

    # Extract all valid points from their LineStrings.
    valid_points_lst = []
    for idx, row in tqdm(annots.iterrows(), total=annots.shape[0], desc='Extracting valid points', unit='valid track', colour='white'):
        valid_points_lst.extend([{'track_id': row._id, 'audible': row.audible, 'geometry': Point(coords)} for coords in row.geometry.coords])
    valid_points = gpd.GeoDataFrame(data=valid_points_lst, geometry='geometry', crs=annots.crs)

    # Load 3D active space for all gains
    model = load_layered_activespace(project_dir, unit, site, year, crs=utm_zone)
    model.preload_all_activespaces()  # load all gains into memory

    # Precompute area of the active space's middle layer, for each gain value
    areas = {}
    layers = list(model.layer_dirs.keys())
    mid_layer = layers[len(layers) // 2]
    logger.info(f"Using {mid_layer}m for areas")
    for gain in tqdm(model.all_activespaces, desc="Computing areas"):
        active = model.all_activespaces[gain][mid_layer]
        active_equal_area = active.to_crs("epsg:3338")  # albers AK equal area
        areas[gain] = active_equal_area.area.item() / 1e6  # square km

    # Precompute whether each point is within each active space. This avoids doing a ton of redundant computation later
    # For each active space gain, we add a boolean column to valid_points representing this
    inside = {}
    for gain in tqdm(model.all_activespaces, desc="Computing points inside/outside"):
        model.set_gain(gain)
        inside[f"in_AS_{gain}"] = model.predict(valid_points)
    valid_points = pd.concat([valid_points, pd.DataFrame(inside)], axis=1)
    
    # Keep adding tracks in chronological order, and see how the fits change
    result_rows = []
    for n in tqdm(range(1, 1+len(track_order)), desc="Fitting each # tracks"):
        track_ids = track_order.iloc[0:n]
        pts = valid_points[valid_points["track_id"].isin(track_ids)]
    
        # Reduce point density to median density, so very dense areas (e.g. airports) don't skew the fit
        # This takes a while relatively to the other things, but we have to redo it for each n value
        # Use a null context manager to suppress print statements for cleanliness
        with open(os.devnull, "w") as fnull:
            with contextlib.redirect_stdout(fnull):
                pts = normalize_point_density(pts, study_area, random_seed=679)
    
        # Compute F1 for each gain
        f1s = {}
        for gain in model.all_activespaces:
            in_AS = pts[f"in_AS_{gain}"].values
            audible = pts["audible"].values
            TP = np.all([in_AS, audible], axis=0).sum()
            FP = np.all([in_AS, ~audible], axis=0).sum()
            FN = np.all([~in_AS, audible], axis=0).sum()
            if TP == 0:
                f1 = 0.0
            else:
                precision = TP / (TP + FP)  # specificity... if a flight enters the active space, is it actually audible?
                recall = TP / (TP + FN)  # sensitivity... if a flight is audible, does it enter the active space?
                f1 = (2 * precision * recall) / (precision + recall)
            f1s[gain] = f1
        
        # Store F1s as a series indexed by gain. Sort to facilitate ranking (see below)
        f1s = pd.Series(f1s).sort_values(ascending=False)

        # For additional robustness, we would like to know gain/f1/area/etc. for the best fit, 2nd best fit, 3rd best, etc.,
        # and then later examine say the top 10 best gains to evaluate whether they have similar gains/areas/etc or chaotically different.

        # Create column prefixes ["1th_best", "2th_best", "3th_best", "4th_best", "5th_best", ...], never mind "1th" vs. "1st"
        rank_strs = ((np.arange(len(f1s)) + 1).astype(str) + "th_best").tolist()

        # Get values for columns ["1th_best_f1", "2th_best_f1", ...]
        ranked_f1s = {f"{s}_f1": f1 for (s, f1) in zip(rank_strs, f1s.values)}

        # Get values for columns ["1th_best_gain", "2th_best_gain", ...]
        gains = f1s.index
        ranked_gains = {f"{s}_gain": g for (s, g) in zip(rank_strs, gains)}

        # Get values for columns ["1th_best_area", "2th_best_area", ...]
        ranked_areas = {f"{s}_area": areas[g] for (s, g) in zip(rank_strs, gains)}
        
        # Combine all info into a row for the final dataframe
        result_rows.append({
            "n_tracks": n,
            "n_annotations": len(annots[annots["_id"].isin(track_ids)]),
            "n_points": len(pts),
            "gain": f1s.index[0],
            "area": areas[f1s.index[0]],
            "area_altitude_m": mid_layer,
            "f1": f1s.iloc[0]
        } | ranked_gains | ranked_f1s | ranked_areas)
    
    # Convert results to proper dataframe and save
    results = pd.DataFrame(result_rows)
    logger.info(f"Results:\n{results}")

    os.makedirs(os.path.dirname(savefile), exist_ok=True)
    results.to_csv(savefile, index=False)
    logger.info(f"Saved to {savefile}")



def plot_stability(project_dir, deployments, col="area", top_k=10, max_n_tracks=None):
    """
    Makes a plot showing whether gain/area/f1/etc has stabilized, for one or more deployments.

    If only one deployment is specified, the plot is saved in the VALIDATION output folder in that site directory.
    Otherwise, the plot is saved in the project directory.

    Parameters
    ----------
    project_dir: str
        Path to the project directory, containing site folders.
    deployments: List[str]
        A list of deployments to plot. E.g. ["DENATRLA2024", "DENATRLA2025", "DENAUWBT2023"]
    col: str, default "area"
        Which column in the study duration stability CSV (see fit_varying_n_tracks()) to plot.
        Common values are "area", "gain", "f1".
    top_k: int, default 10
        The range of the top k best fits will be plotted around the best fit.
    max_n_tracks: int, default None
        If specified, the rightmost x limit for plotting; maximum number of tracks to plot.
        This is useful when comparing deployments where one deployment has many more tracks than the others.
    """
    nrows = len(deployments)
    fig, ax = plt.subplots(nrows=nrows, ncols=1, figsize=(16, 2.5*nrows), sharex=True, sharey=True)
    if nrows == 1:
        ax = [ax]
    
    # determine title and y axis description
    if col == "area":
        title = "Area vs. Number of Tracks Used"
        y_desc = "Area (sq km)"
    elif col == "gain":
        title = "Gain vs. Number of Tracks Used"
        y_desc = "Gain (dB)"
    elif col == "f1":
        title = "F1 vs. Number of Tracks Used"
        y_desc = "F1 Socre"
    else:
        title = f"{col} vs. Number of Tracks Used"
        y_desc = col

    plt.suptitle(title)

    for i, deployment in enumerate(deployments):
        df = pd.read_csv(get_csv_name(project_dir, deployment))

        # look at top k best fits, find range of gain/f1/area
        if col in ["gain", "f1", "area"]:
            top_k_cols = (np.arange(top_k) + 1).astype(str) + f"th_best_{col}"
            top_k_vals = df[top_k_cols]
            min_vals = top_k_vals.min(axis=1)
            max_vals = top_k_vals.max(axis=1)
            ax[i].fill_between(df["n_tracks"], min_vals, max_vals, color="lightblue")

        # plot # tracks vs. value of interest
        ax[i].plot(df["n_tracks"], df[col], label=deployment)
        ax[i].xaxis.set_major_locator(MultipleLocator(50))
        ax[i].xaxis.set_minor_locator(MultipleLocator(10))

        # add y label, including which elevation was used for computing area, if applicable
        ylabel = f"{deployment}\n{y_desc}"
        if col == "area":
            ylabel += f" at {int(df.iloc[0]['area_altitude_m'])}m"
        ax[i].set_ylabel(ylabel)
        
        # apply rightmost xlim
        if max_n_tracks is not None and df["n_tracks"].max() > max_n_tracks:
            plt.xlim(right=max_n_tracks)
    
    # plt.gca().set_xscale("log")
    plt.xlabel("Number of Tracks Used")
    plt.tight_layout()

    # save result
    # if only one deployment, save it in that site's folder
    # otherwise, put it in the project directory
    if len(deployments) == 1:
        unit, site = deployments[0][:4], deployments[0][4:-4]
        filename = os.path.join(project_dir, unit+site, "Output_Data", "VALIDATION", f"{deployment}_{col}_stability.png")
    else:
        filename = os.path.join(project_dir, f"{col}_stability.png")
    plt.savefig(filename)
    logger.info(f"Saved plot to {filename}")