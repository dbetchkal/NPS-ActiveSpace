#!/usr/bin/env python
# coding: utf-8

# # `NPS-ActiveSpace.analysis` module sketch
# purpose: **given an active space, vehicle tracks (as a subset of a transportation system), compute acoustic metrics attributive to the given subset**
# v0.0.1 alpha <br>
# TLogan@nps.gov, DBetchkal@nps.gov

# In[4]:


# GEOSPATIAL LIBRARIES
import gdal
from gdalconst import GA_ReadOnly
import geopandas as gpd
import pyproj
import rasterio
import rasterio.plot
from rasterio.windows import Window
from shapely.geometry import Point, LineString
from scipy.interpolate import UnivariateSpline

# OTHER LIBRARIES
import datetime as dt
from datetime import timedelta
import geopy as geopy
from geopy.distance import geodesic
import glob
import ipykernel
from itertools import islice
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.widgets import RangeSlider, Button, Slider, RadioButtons
import numpy as np
from operator import sub
import os
import pandas as pd
import pickle
import re
from scipy import interpolate
from shapely.geometry import LineString
from shapely.prepared import prep
import sqlalchemy
import subprocess
import sys
import time
from time import mktime
from tqdm import tqdm
import warnings

# SPECIALIZED LIBRARIES
# =======================================================================================================
repo_dir = r"C:\Users\talogan\PythonScripts\3_GITHUB_REPOSITORIES" # ADJUST TO YOUR LOCAL DIRECTORY
# =======================================================================================================
sys.path.append(r"C:\Users\talogan\PythonScripts\GITHUB_REPOS\NPS-ActiveSpace-main")
sys.path.append(r"C:\Users\talogan\PythonScripts\GITHUB_REPOS\NPS-ActiveSpace-main\_DENA\resource")
sys.path.append(r"C:\Users\talogan\PythonScripts\GITHUB_REPOS\NPS-ActiveSpace-main\nps_active_space\utils")

import nps_active_space
from helpers import query_adsb, query_tracks
from computation import interpolate_spline, coords_to_utm

#path to load tracks from api
RDS = r"\\inpdenaterm01\overflights"
sys.path.append(os.path.join(RDS, "scripts"))
from query_tracks import query_tracks, get_mask_wkt


def get_aspect(ax):
    
    """
    Return the aspect ratio of a given `matplotlib` axis object.
    """
    
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

    return disp_ratio / data_ratio


def load_activespace(u, s, y, gain, crs=None, PROJ_DIR=r"V:\NMSim\01 SITES"):
    
    '''The crucial input is an active space estimate, which will be taken as a 'given' logical condition in space.
       Provide a generalized loading process for this input
       with optional co-ordinate reference system conversion using `crs` parameter.
       
       Returns
       -------
       study_area : `gpd.GeoDataFrame` of user-defined study area polygon with designator metadata
    '''

    if gain < 0:
        sign = "-"
    else:
        sign = "+"
        
    if np.abs(gain) < 10:
        gain_string = "0" + str(np.abs(int(10*gain)))
    else:
        gain_string = str(np.abs(int(10*gain)))
    
    # find the estimate
    path = PROJ_DIR + os.sep + u+s + os.sep + u+s+str(y) + '_O_' + sign + gain_string + ".geojson"
    active_space = gpd.read_file(path) # read it in as a `gpd.GeoDataFrame`
    
    if crs is not None:
        active_space = active_space.to_crs(crs) # optional, user-defined coordinate reference system
    
    return active_space

def load_studyarea(u, s, y, crs=None, PROJ_DIR=r"V:\NMSim\01 SITES"):
    
    '''An important input is the Study Area polygon.
       Provide a generalized loading process 
       with optional co-ordinate reference system conversion using `crs` parameter.
       
       Returns
       -------
       study_area : `gpd.GeoDataFrame` of user-defined study area polygon with designator metadata
    '''

    # load in the study area as well
    study_area_path = glob.glob(PROJ_DIR + os.sep + u+s + os.sep + u+s + '*study*area*.shp')[0]
    study_area = gpd.read_file(study_area_path)
    
    if crs is not None:
        study_area = study_area.to_crs(crs)
        
    return study_area


def plot_by_id(col,id, df,active):
    
    """
    TO DO: write Docstring
    """
    
    condition = df[col] == id
    df_condition = df[condition]
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    df_condition.plot(ax=ax, color='red', zorder=1)
    active.boundary.plot(ax=ax, color="black", zorder=1)

    
def internal_interpolate_spline(points: 'Tracks', ds: int = 1) -> gpd.GeoDataFrame:
    
    """
    Interpolate points with a cubic spline between flight points, if possible.
    See https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html#spline-interpolation for docs

    Parameters
    ----------
    points : Tracks
        A Track gpd.GeoDataframe object containing known track points in a path. A minimum of 2 points is required.
    ds : int, default 1
        The second interval in which to calculate the spline for.
        E.g. ds = 1 is "calculate a spline point at every 1 second delta"

    Returns
    -------
    gpd.GeoDataFrame of all points in the interpolated spline.
    Columns: point_dt, geometry

    Raises
    ------
    AssertionError if there is fewer than 1 Track point.
    """
    # Calculate the order of polynomial to fit to the spline. The maximum is a cubic spline. If there are fewer than
    #  3 points, a cubic spline cannot be fit and lower order must be chosen.
    assert points.shape[0] > 1, "A minimum of 2 points is required for calculate a spline."
    k = min(points.shape[0] - 1, 3)

    points.sort_values(by='point_dt', ascending=True, inplace=True)
    starttime = points.point_dt.iat[0]
    endtime = points.point_dt.iat[-1]
    flight_times = (points.point_dt - starttime).dt.total_seconds().values  # Seconds after initial point

    coords = [points.geometry.x, points.geometry.y, points.z] if 'z' in points else [points.geometry.x, points.geometry.y]
    tck, u = interpolate.splprep(x=coords, u=flight_times, k=k)

    # Parametric interpolation on the time interval provided.
    duration = (endtime - starttime).total_seconds()
    tnew = np.arange(0, duration + ds, ds)
    spl_out = interpolate.splev(tnew, tck)
    track_spline = gpd.GeoDataFrame({'point_dt': [starttime + dt.timedelta(seconds=offset) for offset in tnew]},
                                    geometry=[Point(xyz) for xyz in zip(spl_out[0], spl_out[1])],
                                    crs=points.crs)
    return track_spline


def boundary_estimator(trajectory, cap_index, active_space_poly, id, date_time, t_avg, v_rms, key_order):
    
    '''
    Estimate the boundary point of a trajectory within an active space polygon based on certain parameters.
    
    This function estimates the boundary point of a given trajectory within an active space polygon. The boundary point is extrapolated by extending the trajectory line by a certain distance, considering the velocity of the trajectory and the specified time average.
    
    Parameters
    ----------
    trajectory : pandas.DataFrame
        A DataFrame containing the trajectory points with geometries.
    cap_index : int
        An index indicating whether the cap is an origin (0) or terminus (1).
    active_space_poly : geopandas.GeoDataFrame
        A GeoDataFrame representing the active space polygon.
    id : int
        An identifier for the flight.
    date_time : datetime
        The initial datetime of the trajectory.
    t_avg : int
        The time average used for extrapolation calculation (in seconds).
    v_rms : float
        The root mean square velocity of the trajectory.
    key_order : numeric array
        The array for identifying the entrence(0), exit(1) order
        
    Returns
    -------
    new_row : dict
        A dictionary containing flight information including the estimated boundary point geometry.
    boundary_gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing the new_row information for spatial operations.
    '''
    print("In estimator")
    
    aS = active_space_poly
    max_steps = 500
    
    # a cap is either an origin (i.e., 0)
    # in which case we 
    if(cap_index == 0):
        orientation = -1
        basis = trajectory.iloc[0:2]
    
    # ...or a terminus (i.e., 1)
    elif(cap_index == 1):
        orientation = 1
        basis = trajectory.iloc[-3:-1]


    traj_vector = np.array([basis.iloc[-1].geometry.x - basis.iloc[-2].geometry.x, 
                            basis.iloc[-1].geometry.y - basis.iloc[-2].geometry.y])
    
    
    
    
    print("\tExtrapolating to boundary... this may take a while...")
    

    # # calculate furthest end point
    # max_range_extender = orientation * v_rms * max_steps
    # end_point = Point(basis.iloc[-1].geometry + max_range_extender * (traj_vector / np.linalg.norm(traj_vector)))


    # # create a line from basis to end point
    # trajectory_line = LineString([basis.iloc[-1].geometry, end_point])

    # # find the points that intersect
    # intersection = aS.geometry.boundary.unary_union.intersection(trajectory_line)
    # # 2 represents basis point, 0 represents entrence point, 1 represents exit point


        # Create a single trajectory line using all the points in the trajectory
    current_trajectory_line = LineString(trajectory['geometry'].tolist())

    # Calculate the direction of the trajectory based on the last two points
    traj_vector = np.array([trajectory.iloc[-1].geometry.x - trajectory.iloc[-2].geometry.x, 
                            trajectory.iloc[-1].geometry.y - trajectory.iloc[-2].geometry.y])

    # Extrapolate the trajectory
    max_range_extender = v_rms * max_steps
    end_point = Point(trajectory.iloc[-1].geometry.x + max_range_extender * (traj_vector[0] / np.linalg.norm(traj_vector)),
                    trajectory.iloc[-1].geometry.y + max_range_extender * (traj_vector[1] / np.linalg.norm(traj_vector)))

    # Create the extrapolated trajectory line
    extrapolated_trajectory_line = LineString(current_trajectory_line.coords[:] + [end_point])

    # Check for intersections with the active space polygon's boundary (including holes)
    intersection = aS.geometry.boundary.unary_union.intersection(extrapolated_trajectory_line)

    
    if intersection.type == "MultiPoint":
        trajectory_lines = list()
        

        #append basis point to new frame
        # trajectory_lines.append({'flight_id': id,'DateTime': date_time, 'geometry' : basis.iloc[-1].geometry, 'direction': 2})


        # proportion = basis.iloc[-1].geometry.distance(intersection[0]) / basis.iloc[-1].geometry.distance(end_point)
        # t = steps * proportion
        # next_time = date_time + (orientation * (timedelta(seconds=t_avg) * t))
        # trajectory_lines.append({'flight_id': id,'DateTime': date_time,'geometry' : intersection[0], 'direction': key_order[0]})

        
        for i in range(0, len(intersection)):
            proportion = basis.iloc[-1].geometry.distance(intersection[i]) / basis.iloc[-1].geometry.distance(end_point)
            steps = max_steps * proportion
            # Get the exact time and point
            date_time = date_time + (orientation * (timedelta(seconds=t_avg ) * steps))
            new_row = {'flight_id': id, 'DateTime': date_time, 'geometry': intersection[i], 'direction': key_order[i%2] }
            trajectory_lines.append(new_row)
        
        return gpd.GeoDataFrame(trajectory_lines, crs=aS.crs)
    else:
        if intersection.is_empty: 
            print("No intersection found")
            return gpd.GeoDataFrame()
        else:

            # Calculate the proportion along the trajectory_line where the intersection occurs
            proportion = basis.iloc[-1].geometry.distance(intersection) / basis.iloc[-1].geometry.distance(end_point)
            steps = max_steps * proportion
            
            # Get the exact time and point
            next_time = date_time + (orientation * (timedelta(seconds=t_avg) * steps))
            new_row = {'flight_id': id, 'DateTime': next_time, 'geometry': intersection, 'direction': key_order[0]}
            return gpd.GeoDataFrame([new_row], crs=aS.crs)
    

def process_group(group, process_func):
    
    '''
    Apply a processing function to a group of data, returning the processed results as a DataFrame.

    This function is designed to process a group of data using the processing function interpolate spline. 
    It checks if the group contains enough data points to perform the processing and then applies the provided function to the group. 
    If thefunction returns a non-empty DataFrame, the processed results are returned. Otherwise, an empty DataFrame is returned.

    Parameters
    ----------
    group : pandas.DataFrame Group
        A group of data to be processed. This should be a subset of a larger DataFrame.
    process_func : callable function
        A processing function that takes the group of data and performs specific computations or transformations. 
        This function should accept the group DataFrame as the first argument and may have additional optional arguments as needed.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the processed results from applying the process_func to the input group of data. If the processing function does not yield valid results, an empty DataFrame is returned.
    '''
    
    if len(group) >= 2:
        spline_df = process_func(group, ds=5)
        if isinstance(spline_df, pd.DataFrame) and not spline_df.empty:
            return spline_df
    return pd.DataFrame()  # return empty df for groups that don't meet the condition


def find_points_by_flight_id(df, ls):
    
    '''
    Find all rows with a give flight_id in list(ls)
    '''
    
    frame = df[df['flight_id'].isin(ls)]
    return frame


def caculate_speed_distance(trajectory):
    
    '''
    Calculate time-average and root mean square velocity from a trajectory.

    This function calculates the time-average velocity (t_avg) and root mean square velocity (v_rms) from a given trajectory. 
    The trajectory should contain spatial and temporal information for each point.

    Parameters
    ----------
    trajectory : pandas.DataFrame
        A DataFrame containing trajectory points with geometry and timestamp information.

    Returns
    -------
    t_avg : float
        The calculated time-average velocity based on the provided trajectory.
    v_rms : float
        The calculated root mean square velocity based on the provided trajectory.
    '''
    
    trajectory["lon"] = trajectory.geometry.x
    trajectory["lat"] = trajectory.geometry.y
        
    # calculate range
    time_range = trajectory["DateTime"].min() - trajectory["DateTime"].max()

    # calc average time difference
    t_avg = time_range.total_seconds() / (len(trajectory) - 1)

    #Calculate the differences in latitude and longitude
    lat_diff = trajectory["lat"].diff()
    lon_diff = trajectory["lon"].diff()

    # calculate euclidian
    distance = np.sqrt(lat_diff**2 + lon_diff**2)

    # Calculate the total
    total_distance = distance.sum()

    # Calculate the time difference between the first and last timestamps
    time_difference = (trajectory["DateTime"].iloc[-1] - trajectory["DateTime"].iloc[0]).total_seconds()

    # Calculate the average velocity
    v_rms = total_distance / time_difference
    return t_avg, v_rms


def spatial_mitigation_metrics(intersections, active_space): 
    
    '''
    Calculate spatial mitigation metrics: median site distance and median inaudibility distance.

    This function calculates two spatial mitigation metrics based on the intersections of trajectories with an active space polygon. 
    The first metric is the median site distance, which measures the minimum distance from each intersection point to the boundary of the active space. 
    The second metric is the median inaudibility distance, which quantifies the distance from each intersection point to the exterior boundary of the active space polygon.

    Parameters
    ----------
    intersections : geopandas.GeoDataFrame
        A GeoDataFrame containing intersection points between trajectories and an active space polygon.
    active_space : geopandas.GeoDataFrame
        A GeoDataFrame representing the active space polygon.

    Returns
    -------
    median_site_distance : float
        The median value of the minimum distance from each intersection point to the active space boundary.
    median_inaudibility_distance : float
        The median value of the distance from each intersection point to the exterior boundary of the active space polygon.


    '''
    
    #calculate site distance
    site_distance = intersections.geometry.apply(lambda x: active.distance(x).min())

    inaudibility_distance = intersections.geometry.apply(lambda x: x.distance(active.unary_union.boundary))

    return site_distance.median(), inaudibility_distance.median()


def adjust_noise_free_intervals(noise_free_intervals, noise_intervals, begin_time, end_time):

    '''
    To properly compute metrics we need to binarize audibility across 
    the entire evaluation period from beginning to end. This function adjusts
    a record of noise-free intervals considering all possible boundary conditions.
    
    Parameters
    ----------
    noise_free_intervals : `np.array` of `np.int32` timestamps with shape (m,2) representing bounds of noise-free intervals
    noise_intervals : `np.array` of `np.int32` timestamps with shape (m,2) representing bounds of audible noise
    
    Returns
    -------
    noise_free_intervals : `np.array` of `np.int32` timestamps with shape (m,2) updated with correct boundary times
    '''

    nfi_starts = noise_free_intervals.T[0]
    nfi_ends = noise_free_intervals.T[1]

    # ------- Account for different beginning conditions -----------------------------------------

    # the record begins with noise...
    if(noise_intervals[0, 0] == 0):

        # ...the first noise free interval (and thus ALL intervals) need to start one second later
        nfi_starts = nfi_starts + 1


    # the record begins with quietude...
    else:

        # ...the first noise free interval stays the same, and equals the beginning timestamp
        # the rest are + 1
        nfi_starts = nfi_starts + 1
        nfi_starts[0] = begin_time


    # ------- Account for different ending conditions -----------------------------------------

        # the record ends with noise...
    if(noise_intervals[-1, 0] == 0):

        # ...the last noise free interval (and thus ALL intervals) need to end one second earlier
        nfi_ends = nfi_ends - 1


    # the record ends with quietude...
    else:

        # ...the last noise free interval stays the same, and equals the ending time stamp
        # the rest are - 1
        #save = nfi_ends[-1]
        nfi_ends = nfi_ends - 1
        nfi_ends[-1] = end_time

    # reset attribute to these updated, correct values
    new_noise_free_interval = np.array([nfi_starts, nfi_ends]).T
    
    return new_noise_free_interval


def contiguous_regions(condition):

    """
    Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    
    """

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero() 

    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)

    return idx


def summarize_duration(duration_list, dur_threshold=None):

    '''
    Given a set of interval durations [either overlapping or non-overlapping],
    generate a standardized battery of duration summaries:
    
        Mean, standard deviation, median, and median absolute deviation.
    
    Paramters
    ---------
    `duration_list` : (m,) `np.array` representing event durations
    `dur_threshold` : (float) defaults to None, if passed two additional outputs are generated. See Returns, below.
    
    Returns
    -------
    `out` : (6,) `np.array` of duration metrics:
        [0] the mean duration
        [1] the standard deviation of the durations
        [2] the median duration
        [3] the median absolute deviation of the durations
        [4] Optional: the percentage of intervals greater than the user-provided duration threshold
        [5] Optional: the number of intervals greater than the user-provided duration threshold
    '''

    # mean duration
    mean = np.mean(duration_list)

    # standard deviation duration
    stdev = np.std(duration_list)

    # median duration
    median = np.percentile(duration_list, 50)

    # median absolute deviation of duration
    mad = np.percentile(np.absolute(duration_list - median), 50)
    
    if(dur_threshold is not None):
        
        # how many events exceed the user-defined threshold?
        exceeding = duration_list[duration_list > dur_threshold]
        nExceed = len(exceeding)
        pExceed = nExceed/len(duration_list) # percentage on [0, 1]
        
        # combine the results
        out = np.array([mean, stdev, median, mad, nExceed, pExceed])

    else:
        # combine the results
        out = np.array([mean, stdev, median, mad, np.nan, np.nan])

    # it's convenient to return the results
    return out


# ## Step 1: Initialize global variables (for eventual `analysis` class)

# In[5]:


# A microphone deployment,
# is designated in three parts: UNIT, SITE, YEAR
# we also include a GAIN (in decibels)

# ===== A FEW TEST SCENARIOS ====================================
# from:  "V:\NMSim\2022 11 21 NPS-ActiveSpace gain results.xlsx"
# ===============================================================
# u, s, y, gain = "DENA", "RUTH", 2019, -9.5
# u, s, y, gain = "DENA", "KAHP", 2018, 22.5
# u, s, y, gain = "DENA", "TRLA", 2019, 15.0
u, s, y, gain = "DENA", "TRLA", 2021, 3.5
# u, s, y, gain = "DENA", "TOKO", 2021, 3.5 # our default test scenario
# u, s, y, gain = "DENA", "NTRL", 2019, 29.5

# Define the interval upon which to compute acoustic metrics ⊆ times where there are tracks
start_date, end_date = "2021-05-01", "2021-08-01" # notice this is a different year than fit!
draw_plots = True


# ## Step 2: Parse Data + Establish Geographic Intersection
# Geographic bugs remain as of v0.0.1

# In[6]:


print("========== NPS-ActiveSpace Analysis Module ==========")
print("           v0.0.1 alpha\n\n")
print("\tAnalyzing", u+s+str(y), "fit for a Gain of", gain, "dB:")
print("\t\tLoading geospatial inputs...")
start_full = time.time()
start = time.time()

# load the active space
# matching the above parameters
active = load_activespace(u, s, y, gain, crs="epsg:4326")

# load the study area
# matching the above parameters 
studyA = load_studyarea(u, s, y, crs="epsg:4326")

if(draw_plots):
    
    # an overview of the two requisite spatial inputs 
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    studyA.boundary.plot(ax=ax, color="blue", ls="--")
    active.boundary.plot(ax=ax, color="black", zorder=1)
    ax.set_title("Geospatial inputs for " + u+s+str(y) + " at Gain " + str(gain) + " dB", loc="left")
    simple_aspect = get_aspect(ax)
    plt.show()

end = time.time()
print(f"\t\t\tSite layers loaded in {(end - start):.3f} seconds.\n")
print("\t\tLoading and cleaning track data...")
start = time.time()

# initialize the `sqlalchemy` engine
engine = sqlalchemy.create_engine('postgresql://{username}:{password}@{host}:{port}/{name}'.format(username="overflights_admin", 
                                                                                                   password="0verfl!ghts",
                                                                                                   host="165.83.50.64",
                                                                                                   port="5432",
                                                                                                   name="overflights"))
# load tracks from the database over a certain daterange, using the buffered site
tracks = query_tracks(engine=engine,
                      start_date=start_date, end_date=end_date,
                      mask=studyA)

tracks = tracks.set_crs('WGS84') # define the spatial reference manually


# In[7]:



# ensure field names are uniform across all track sources
tracks.rename(columns={'ak_datetime': 'point_dt'}, inplace=True)
tracks.rename(columns={'geom' : 'geometry'}, inplace=True)

end = time.time()
print(f"\t\t\tRaw track data parsed in {(end - start):.3f} seconds.\n")
print("\t\tSplining track data...")
start = time.time()

# to spline, we first group by flight ID
grouped = tracks.groupby('flight_id')
result_dataframes = grouped.apply(process_group, process_func=internal_interpolate_spline)

# clean the data frame of empty rows
cleaned_df = result_dataframes.dropna()
cleaned_df = cleaned_df.reset_index()
cleaned_df = gpd.GeoDataFrame(cleaned_df) # apparently we need to coerce in some cases?

end = time.time()
print(f"\t\t\tTracks splined in {(end - start):.3f} seconds.\n")
print("\t\tComputing track intersections with active space...\n")
start = time.time()

#identify intersections
intersections = gpd.overlay(cleaned_df, active, how='intersection')
intersections = intersections.drop(columns=['level_1', 'altitude_m', 'mic_name'])
# intersections = intersections.rename(columns={'point_dt': 'DateTime'}) #              WHY?
intersections.dropna(axis=1, how='any')

if(draw_plots):
    
    fig, ax = plt.subplots(1,1, figsize=(6,6))
    studyA.boundary.plot(ax=ax, color="blue", ls="--")
    intersections.plot(ax=ax, color="red", markersize=1, lw=0.5, zorder=-2)
    active.boundary.plot(ax=ax, color="black", zorder=1)
    ax.set_title("Pre-processed track intersections", loc="left")
    plt.show()
    
# glean the UTM zone using the active space centroid
# `geopandas` developers will helpfully warn us 
# about geographic co-ordinate systems, here
# but it's OK - we suppress the warning...
warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS.")
x_coord,y_coord = active["geometry"].centroid[0].xy 
utm_zone = coords_to_utm(lat=y_coord[0], lon=x_coord[0])

# to accurately compute distance, we need to be in an equal-area projection
intersections_ea = intersections.to_crs(utm_zone)
tracks_ea = tracks.to_crs(utm_zone)
grouped_ea = intersections_ea.groupby('flight_id')
active_ea = active.to_crs(utm_zone)


# In[8]:


intersections_ea.rename(columns={'point_dt': 'DateTime'}, inplace=True)


# In[9]:


origin = gpd.GeoDataFrame([])
terminus = gpd.GeoDataFrame([])
intersections_ea['direction'] = 2

count = 0
for j, frame in grouped_ea:
    count += 1
    if frame.empty or len(frame['geometry']) < 2: #                CAN YOU EXPLAIN THESE CONDITIONS?
        print("Empty")
        continue
    else:
        # every interval has two caps (or boundaries)

        start_cap = frame['geometry'].iloc[0].within(active_ea.geometry.unary_union)
        cap_end = frame['geometry'].iloc[-1].within(active_ea.geometry.unary_union)
        print(start_cap, cap_end )

        t_avg, v_rms = caculate_speed_distance(frame)
        print('\n\t\t\tFlight ID: ', j)
        
        if start_cap:
            start = time.time()
            print("\t\t\tOrigin at", frame["DateTime"].iloc[0], "is inside the active space. Extrapolation required!\n")
            print("\t\t\t\tWe extrapolate the first two points backwards in time:\n")

            fwd = boundary_estimator(frame, cap_index=0, active_space_poly=active_ea, id=j,
                                          date_time=frame["DateTime"].iloc[0], t_avg=t_avg, v_rms=v_rms, key_order=[0,1])
            origin = origin.append(fwd, ignore_index=True)
            
            # intersections_ea = intersections_ea.append(origin, ignore_index=True)

            
            end = time.time()
            print(f"\t\t\tCalculating origin took {(end - start):.3f} seconds.")
        else:
            print("\t\t\tOrigin at", frame["DateTime"].iloc[0], "is outside the active space. No extrapolation needed.")

        if cap_end:
            start = time.time()
            print("\t\t\tTerminus at", frame["DateTime"].iloc[-1], "is inside the active space. Extrapolation required!\n")
            print("\t\t\t\tWe extrapolate the last two points forward in time:\n")
            
            bkwd = boundary_estimator(frame, cap_index=1, active_space_poly=active_ea, id=j,
                                          date_time=frame["DateTime"].iloc[len(frame) - 1], t_avg=t_avg, v_rms=v_rms, key_order = [1, 0])
            
            terminus = terminus.append(bkwd, ignore_index=True)
            
            # intersections_ea = intersections_ea.append(terminus, ignore_index=True)
            
            end = time.time()
            print(f"\t\t\tCalculating origin took {(end - start):.3f} seconds.")

        else:
            print("\t\t\tTerminus at",frame.iloc[-1]["DateTime"], "is outside the active space. No extrapolation needed.")
            # where are the entries and exit points along the active space?
        fig, ax = plt.subplots(1,1, figsize=(8,8))
        active_ea.boundary.plot(ax=ax, color="black", zorder=-1)
        # sA.boundary.plot(ax=ax, color="blue", ls="--", zorder=-2)
        frame.plot(ax=ax, markersize=10, color='blue', marker="o", zorder=-1)
        origin[origin['flight_id']  == j].plot(ax=ax, markersize=40, c=origin['direction'], marker="x", zorder=1, label="entry points")
        terminus[terminus['flight_id']  == j].plot(ax=ax, markersize=40, c=terminus['direction'], marker="x", zorder=1, label="exit points")
        ax.set_title("Post-processed track intersections", loc="left")
        ax.legend()
        # ax.set_aspect(simple_aspect)
        plt.show()
        plt.close()


# In[10]:


origin


# In[ ]:





# In[11]:


# create a temporary column
intersections_ea['direction'] = 2

#merge all data frames
intersections_ea = intersections_ea.append(origin, ignore_index=True)
intersections_ea = intersections_ea.append(terminus, ignore_index=True)


intersections_ea


# In[12]:


# extract just the entry and exit points
intersections_ea = intersections_ea.to_crs(active.crs)
# along with their timing
entries = []
entry_times = []
exits = []
exit_times = []
for idx, intersect in intersections_ea.iterrows():
    #entry times
    if intersect.direction == 0:
        entries.append(intersect.geometry)
        entry_times.append(intersect.DateTime)
    #exit times
    elif intersect.direction == 1:
        exits.append(intersect.geometry)
        exit_times.append(intersect.DateTime)
    else:
        continue
intersections_ea = intersections_ea.to_crs(utm_zone)


# In[13]:


print(len(exits), len(entries))


# In[14]:


grouped = intersections_ea.groupby('flight_id')
for id, group in grouped:
    local_points = []
    group.sort_values(by='DateTime',inplace=True)
    print(group.DateTime, group.direction)


# In[15]:


grouped = intersections_ea.groupby('flight_id')


lines = []

for id, group in grouped:
    local_points = []
    group.sort_values(by='DateTime',inplace=True)
    print(group)
    for g in group.geometry:
        coords = np.array(g)
        if g.type == 'Point':
            lat_lon = coords.reshape(-1, 2)
        elif g.type == 'Multipoint':
            lat_lon = coords
        local_points.extend([tuple(row) for row in lat_lon])
    
    if len(local_points) >= 2:
        lines.append({
            'flight_id': id, 
            'geometry': LineString(local_points), 
            'DateTime': group.DateTime.values
        })

# Finalize into geopandas dataframe
# and transform back into geographic co-ordinates (WGS-84)
final_intersections_ea = gpd.GeoDataFrame(lines, crs=utm_zone)
final_intersections = final_intersections_ea.to_crs(active.crs)


    
# organize the information as two `gpd.GeoSeries` objects
enter = gpd.GeoSeries(entries, index=entry_times)
exit = gpd.GeoSeries(exits, index=exit_times)
final_intersections = gpd.overlay(final_intersections, active, how='intersection')
final_intersections_ea.to_crs(utm_zone)


if(draw_plots):
    
    # where are the entries and exit points along the active space?
    fig, ax = plt.subplots(1,1, figsize=(8,8))
    active.boundary.plot(ax=ax, color="black", zorder=-1)
    studyA.boundary.plot(ax=ax, color="blue", ls="--", zorder=-2)
    final_intersections.plot(ax=ax, color="gray", markersize=1, lw=0.5, zorder=0)
    enter.plot(ax=ax, markersize=10, color="orange", marker="o", zorder=-1, label="entry points")
    exit.plot(ax=ax, markersize=10, color="magenta", marker="x", zorder=-1, label="exit points")
    ax.set_title("Post-processed track intersections", loc="left")
    ax.legend()
    ax.set_aspect(simple_aspect)
    plt.show()

end = time.time()
print(f"\n\t\t\tIntersections have been determined! Took {(end - start):.3f} seconds.\n")

print("\n\n\n\nEVENTUALLY THE ACOUSTIC METRIC COMPUTATIONS WILL HAPPEN HERE\n\n\n\n")

end_full = time.time()
print(f"\n\tAnalysis complete. The entire process took {(end_full - start_full):.3f} seconds.\n")


# In[16]:


intersections_ea[intersections_ea['flight_id'] == 'N185AR_202106111515']


# In[17]:


final_intersections


# In[18]:


print(final_intersections.iloc[0].geometry)


# In[19]:


def calculate_duration(datetime_array):
    duration = pd.to_timedelta(max(datetime_array) - min(datetime_array))
    return duration.total_seconds()

# Calculate duration for each flight
final_intersections['duration_seconds'] = final_intersections['DateTime'].apply(calculate_duration)


# In[20]:


def compute_metrics(column, quantile, df):
    print(f"Quantiles:\n{df[column].quantile(quantile)}\n")
    print(f"Mean Duration: {df[column].mean()}")
    print(f"Mean Absolute Deviation: {df[column].mad()}")
    print(f"Standard Deviation: {df[column].std()}")


# In[21]:


print("Noise Event Duration")
compute_metrics('duration_seconds', [0.25, 0.5, 0.75], final_intersections)


# In[22]:


chunk = final_intersections[final_intersections['flight_id'] == 'N185AR_202105251000']
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.set_ylim([63.65, 63.80])
ax.set_xlim([-149.2, -148.825])
rand_colors = np.random.random(size=(len(chunk), 3))
active.boundary.plot(color="k", ax=ax, label="active space")
chunk.plot(color=rand_colors, ax=ax, markersize=50, label="input trajectory")
enter.plot(ax=ax, markersize=10, color="orange", marker="o", zorder=2, label="entry points")
exit.plot(ax=ax, markersize=50, color="magenta", marker="x", zorder=1, label="exit points")
ax.legend(bbox_to_anchor=(1.55, 1.0))
plt.show()


# In[ ]:





# ## Step 3: Compute (Attributive) Acoustic Metrics 
# 
# We use the term *attributive* because we are computing these metrics given tracks potentially representing a subset of the overall transportation system. We expect our results to differ from cumulative analyses - such as results derived from an acoustic record - in proportion to the subset's throughput when compared to the transportation system throughput as a whole.
# 
# ### We must be capable of performing analyses at *three temporal scales*: <font color="salmon">overall</font>, <font color="orange">by day</font>, and <font color="green">by hour</font>.
# 
# ----
# 
# ## Metric types:
# 
# #### Spatial mitigation
# 
# >median distance from site <br> 
# >median distance from inaudibility <br>
# 
# #### Duration
# precursor `duration_list`, return:
# >noise event duration list <br>
# >quantile, mad, mean, std Event duration
# 
# 
# #### Noise-free Intervals
# precursor `NFI_list`, return:
# >noise-free interval duration list <br>
# >quantile, mad, mean, std NFI duration <br>
# >longest NFI in a given summary period <br>
# >longest NFI in a given summary period exceeds $t_{thresh}$
# 
# 
# #### Audibility
# >overall TA (%) <br>
# >quantile, mad Daily TA (%) <br>
# >mean, std Daily TA (%) <br>
# >quantile, mad Hourly TA (%) <br>
# >mean, std Hourly TA (%)
# 
# #### Event Rates
# >overall Event count <br>
# >quantile, mad Daily Events <br>
# >mean, std Daily Events <br>
# >quantile, mad Hourly Events <br>
# >mean, std Hourly Events
# 
#  
# #### Specific Park Standards
# >ATMP-related standards? <br>
# >ZION standards <br>
# >DENA PA standard <br>
# >DENA Events standard

# In[ ]:


import geopandas as gpd
from shapely.geometry import MultiLineString

gdf = final_intersections

for index, row in gdf.iterrows():
    geometry = row['geometry']
    date_times = row['DateTime']

    if isinstance(geometry, MultiLineString):
        total_noise_free_duration = 0
        for i in range(len(geometry) - 1):
            end_point_datetime = date_times[len(geometry[i].coords) - 1]
            
            start_point_datetime_next_line = date_times[len(geometry[i].coords)]
            
            noise_free_interval = start_point_datetime_next_line - end_point_datetime
            
            noise_free_interval_seconds = noise_free_interval / np.timedelta64(1, 's')
            total_noise_free_duration += noise_free_interval_seconds
        
        total_duration = (date_times[-1] - date_times[0]) / np.timedelta64(1, 's')
        
        percentage_noise_free = (total_noise_free_duration / total_duration) * 100
        
        print(f"Trajectory {index} is {percentage_noise_free:.2f}% noise-free.")


# In[ ]:


print("\t\tSummarizing results.\n\t\tGenerating boolean audibility array, this may take a while...")
start = time.time()

# convert the boundaries of the evaluation period into true `dt.datetime` objects
start_dt = dt.datetime.strptime(start_date, "%Y-%m-%d")
end_dt = dt.datetime.strptime(end_date, "%Y-%m-%d")

# an array containing every second of the record as datetime
all_seconds = np.array([start_dt + dt.timedelta(seconds=i) for 
                        i in range(int((end_dt - start_dt).total_seconds()))])

# organize the intervals into an (n, 2) array
noise_intervals = np.array([[b, e] for b,e in zip(enter.index.to_pydatetime(), 
                                                  exit.index.to_pydatetime())])

# a boolean array the same length as the record
within_times = np.zeros(all_seconds.shape)

# ...and create a `pd.Series` that unites time with audibility status
aud_binary = pd.Series(within_times, index=all_seconds)

end = time.time()
print(f"\t\t\tGenerating the boolean array took {(end - start):.3f} seconds.\n")
print("\t\tTransforming audibility results...")
start = time.time()

for noise_event in tqdm(noise_intervals, unit="intervals"):

    try:

        # convert the enter/exit times into integer indices
        ind_begin = np.argwhere(all_seconds == noise_event[0])[0][0]
        ind_end = np.argwhere(all_seconds == noise_event[1])[0][0]

        # set values between the two indices equal to 'True'
        within_times[ind_begin:ind_end] = 1.0

    except IndexError:
        pass # basically the event had no duration

# finally, we reduce:

# first, reduce by finding the times where vehicles were not in the active space (i.e., inaudibility)
# and do a bit of rearranging to arrive at (geometric) NFI values
nfi_bounds = contiguous_regions(aud_binary.values == 0)
NFI_list = (nfi_bounds.T[1] - nfi_bounds.T[0])

# then, reduce by finding times where vehicles were in the active space (i.e., audibility)
noise_bounds = contiguous_regions(aud_binary.values == 1)
TA_list = (noise_bounds.T[1] - noise_bounds.T[0])

# we'll also summarize the durations of each unique aircraft event
# which is usually not possible with the microphone, alone
simplify_times = lambda t: t.total_seconds()
simple = np.vectorize(simplify_times)
duration_list = simple(noise_intervals.T[1] - noise_intervals.T[0]).astype('int')

end = time.time()
print(f"\t\t\tTransformed in {(end - start):.3f} seconds. TO DO NEXT... TRUE SUMMARIES!\n\n\n\n\n")

# TA_list, duration_list, NFI_list # a list of noise event and noise-free interval durations


# #### A template sketch for one of several summary metric procedural functions
# will eventually be used on both `duration_list` and `NFI_list` at different time scales

# In[ ]:


t_thresh = 8500.0
# mean, standard deviation, median, median absolute deviation, number of intervals exceeding threshold, percentage exceeding
M,std,m,mad,nE,pE = summarize_duration(NFI_list, dur_threshold=t_thresh) 

print("The mean interval lasted {0} ± {1}.\nThe typical (median) interval lasted {2} ± {3}.\nA total of {4} intervals ({5:.1f}%) exceeded the {6} threshold.".format(dt.timedelta(seconds=int(M)),
                                                                                                                                                                     dt.timedelta(seconds=int(std)),
                                                                                                                                                                     dt.timedelta(seconds=int(m)),
                                                                                                                                                                     dt.timedelta(seconds=int(mad)),
                                                                                                                                                                     int(nE),
                                                                                                                                                                     100*pE,
                                                                                                                                                                     dt.timedelta(seconds=int(t_thresh))))


# #### A template sketch for finer temporal scales using `pd.Grouper`
# For <font color="orange">**daily** </font> analyses we group using `freq="1 D"`</font> and for <font color="green">**hourly** </font> analyses using `freq="1 H"`</font> 

# In[ ]:


# to implement finer time resolution 
for meta, t_block in aud_binary.groupby(pd.Grouper(freq="1 D", axis=0)):
    
    print("This analysis block begins at", meta)

    # we perform 
    nfi_bounds = contiguous_regions(t_block.values == 0)
    NFI_list = (nfi_bounds.T[1] - nfi_bounds.T[0])

    noise_bounds = contiguous_regions(t_block.values == 1)
    duration_list = (noise_bounds.T[1] - noise_bounds.T[0])
    
    print("\tNFIs")
    print("\t", NFI_list)
    
    print("\tDurations")
    print("\t", duration_list)

