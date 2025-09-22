import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point, LineString
from nps_active_space.utils.computation import interpolate_spline, audible_time_delay, coords_to_utm
from nps_active_space.utils.models import Annotations, Adsb, Tracks
from _DENA.resource.helpers import get_deployment, load_studyarea
from tqdm import tqdm


def correct_clock_drift(tracks: Tracks, clock_drift_file: str, inplace: bool=True):
    """Fix the point_dt field of a set of tracks by correcting for clock drift.
    
    Parameters
    ----------
    tracks: Tracks
        A GeoDataFrame containing a "point_dt" field, with a row for each point.
        The times in this GeoDataFrame have not been corrected for clock drift.
    clock_drift_file: str
        Path to the clock drift file. This is a CSV file with columns "Time" and "Seconds",
        representing the drift at various points in time. Drift is defined as what you would
        need to add to the track point time to arrive at the corresponding NVSPL time.
        This function linearly interpolates between drifts in the CSV file to find the drift
        at any moment in time. Note that the range of the CSV file times must encompass
        the range of the track point times for interpolation to work.
    inplace: bool
        Whether to modify tracks["point_dt"] in place.
    
    Returns
    -------
    tracks: Tracks
        A GeoDataFrame of track points with times shifted to account for clock drift.
    """

    # read file into a pd.Series
    drifts = pd.read_csv(clock_drift_file, index_col="Time")["Seconds"]
    drifts.index = pd.to_datetime(drifts.index)
    # make sure the clock drifts encompass the track point times, so we can interpolate
    if tracks["point_dt"].min() < drifts.index.min() or tracks["point_dt"].max() > drifts.index.max():
        raise Exception(f"Clock drift corrections must encompass the whole track period ({tracks["point_dt"].min()} - {tracks["point_dt"].max()})")
    
    # add in entries corresponding to the track point times in between existing clock drift entries 
    # then interpolate to fill them in, then extract those interpolated values to get time adjustments
    # make sure there are no duplicate times in drifts_augmented, so that when we query drifts_augmented to get the adjustments,
    # each point_dt value will only correspond to one entry in drifts_augmented
    times_to_interpolate = pd.Series(data=np.nan, index=tracks["point_dt"].unique())
    times_to_interpolate = times_to_interpolate[~times_to_interpolate.index.isin(drifts.index)]
    drifts_augmented = pd.concat([drifts, times_to_interpolate])
    drifts_augmented.sort_index(inplace=True)
    drifts_augmented.interpolate(method="time", inplace=True)
    adjustments = drifts_augmented[tracks["point_dt"]]
    adjustments = pd.to_timedelta(adjustments, unit="s")
    correct_point_dt = tracks["point_dt"] + adjustments.values

    if inplace:
        tracks["point_dt"] = correct_point_dt
        return tracks
    else:
        track_copy = tracks.copy()
        track_copy["point_dt"] = correct_clock_drift
        return track_copy