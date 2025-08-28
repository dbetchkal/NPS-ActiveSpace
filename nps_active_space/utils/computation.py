import datetime as dt
import math
from typing import Iterable, List, Optional, Tuple, TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from osgeo import gdal
from scipy import interpolate
from scipy.spatial import KDTree
from shapely.geometry import Point
from KDEpy import FFTKDE

if TYPE_CHECKING:
    from nps_active_space.utils.models import Microphone, Nvspl, Tracks


__all__ = [
    'audibility_to_interval',
    'ambience_from_nvspl',
    'ambience_from_raster',
    'audible_time_delay',
    'build_src_point_mesh',
    'calculate_duration_summary',
    'climb_angle',
    'compute_fbeta',
    'contiguous_regions',
    'coords_to_utm',
    'create_overlapping_mesh',
    'expected_relative_Lp',
    'interpolate_spline',
    'NMSIM_bbox_utm',
    'normalize_point_density',
    'project_raster'
]


def NMSIM_bbox_utm(study_area: gpd.GeoDataFrame) -> str:
    """
    NMSIM references an entire project to the westernmost extent of the elevation (or landcover) file.
    Given that, return the UTM Zone the project will eventually use. NMSIM uses NAD83 as its geographic
    coordinate system, so the study area will be projected into NAD83 before calculating the UTM zone.

    Parameters
    ----------
    study_area : gpd.GeoDataFrame
        A study area (Polygon) to find the UTM zone of the westernmost extent for.

    Returns
    -------
    UTM zone projection name (e.g.  'epsg:26905' for UTM 5N) that aligns with the westernmost extent of a study area.
    """
    if study_area.crs.to_epsg() != 4269:
        study_area = study_area.to_crs(epsg='4269')
    study_area_bbox = study_area.geometry.iloc[0].bounds  # (minx, miny, maxx, maxy)
    lat = study_area_bbox[3]  # maxy
    lon = study_area_bbox[0]  # minx

    return coords_to_utm(lat, lon)


def coords_to_utm(lat: float, lon: float) -> str:
    """
    Takes the latitude and longitude of a point and outputs the EPSG code corresponding to the UTM zone of the point.

    Parameters
    ----------
    lat : float
        Latitude of a point in decimal degrees in a geographic coordinate system.
    lon : float
        Longitude of a point in decimal degrees in a geographic coordinate system.

    Returns
    -------
    utm_proj : str
        UTM zone projection name (e.g.  'epsg:26905' for UTM 5N)

    Notes
    -----
    Remember: x=longitude, y=latitude
    """
    # 6 degrees per zone; add 180 because zone 1 starts at 180 W.
    utm_zone = int((lon + 180) // 6 + 1)

    # 269 = northern hemisphere, 327 = southern hemisphere
    utm_proj = 'epsg:269{:02d}'.format(utm_zone) if lat > 0 else 'epsg:327{:02d}'.format(utm_zone)
    return utm_proj


def climb_angle(v: Iterable) -> np.ndarray:
    """
    Compute the 'climb angle' of a vector.
    A = 𝑛•𝑏=|𝑛||𝑏|𝑠𝑖𝑛(𝜃)

    Parameters
    ----------
    v : array-like
        Vector to compute the climb angle for.

    Returns
    -------
    degrees : ndarray of floats
        Corresponding climb angle value in degrees.
    """
    n = np.array([0, 0, 1])  # A unit normal vector perpendicular to the xy plane
    degrees = np.degrees(np.arcsin(np.dot(n, v) / np.linalg.norm(n) * np.linalg.norm(v)))
    return degrees


def interpolate_spline(points: 'Tracks', ds: int = 1, s: float = None) -> gpd.GeoDataFrame:
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
    s : float, optional, default None
        The amount of smoothing to apply. Larger `s` means more smoothing while smaller values of `s`
        indicate less smoothing.

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
    tck, u = interpolate.splprep(x=coords, u=flight_times, k=k, s=s)

    # Parametric interpolation on the time interval provided.
    duration = (endtime - starttime).total_seconds()
    tnew = np.arange(0, duration + ds, ds)
    spl_out = interpolate.splev(tnew, tck)
    track_spline = gpd.GeoDataFrame({'point_dt': [starttime + dt.timedelta(seconds=offset) for offset in tnew]},
                                    geometry=[Point(xyz) for xyz in zip(spl_out[0], spl_out[1], spl_out[2])],
                                    crs=points.crs)
    return track_spline


def audible_time_delay(points: gpd.GeoDataFrame, time_col: str, target: Point,
                       m1: float = 343., drop_cols: bool = False) -> gpd.GeoDataFrame:
    """
    Given a set of points and a target location, calculate when a sound made at each point could be heard at
    the target.

    **IMPORTANT**: The points GeoDataFrame and the target Point should be in the same crs for accurate calculations.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        A gpd.GeoDataFrame of sound location points.
    time_col : str
        Name of the column in the points gpd.GeoDataFrame with time of sound occurrence at each point.
    target : Point
        The target point.
    m1 : float, default 343 m/s
        The speed of sound to use for calculations. Make sure this value uses the same units as the crs of
        the points GeoDataFrame and the target Point.
    drop_cols : bool, default False
        If True, drop the intermediate columns used to determine time of audibility.

    Returns
    -------
    The points GeoDataFrame with added columns:
    Standard: time_audible
    Optional: distance_to_target, audible_delay_sec
    """
    points['distance_to_target'] = points.geometry.apply(lambda geom: target.distance(geom))
    points['audible_delay_sec'] = points['distance_to_target'] / m1
    points['time_audible'] = points.apply(lambda row: row[time_col] + dt.timedelta(seconds=row.audible_delay_sec), axis=1)

    if drop_cols:
        points.drop(['distance_to_target', 'audible_delay_sec'], inplace=True)

    return points


def expected_relative_Lp(points: gpd.GeoDataFrame, target: Point):
    """Get expected Lp values for a set of points and a target observer location, using a crude acoustic propagation model.
    Note that these are only relative Lp values, because an arbitrary speaker power is used.
    
    **IMPORTANT**: The points GeoDataFrame and the target Point should be in the same crs for accurate calculations.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        A gpd.GeoDataFrame of sound location points.
    target : Point
        The target point.
    
    Returns
    -------
    The points GeoDataFrame with added columns: Lp_est
    """
    distances = points.geometry.apply(lambda geom: target.distance(geom))
    Lw = 150
    atm_abs = -0.002
    A_geometric = 10*np.log10(1/(4*np.pi*np.power(distances,2)))
    A_atmosphere = np.array(distances)*atm_abs
    points["Lp_est"] = Lw + A_geometric + A_atmosphere
    return points


def build_src_point_mesh(area: gpd.GeoDataFrame, density: int = 48, altitude: Optional[int] = None) -> List[Point]:
    """
    Given a polygon and a density, create a square mesh of evenly spaced points throughout the polygon.

    Parameters
    ----------
    area : gpd.GeoDataFrame
        A GeoDataFrame of the area to create the square point mesh over.
    density : int
        The number of points along each mesh axis. The mesh will contain density x density points.
    altitude : int, default None
        A standard altitude to apply to every point in the mesh.

    Returns
    -------
    mesh points : List[Point]
        A list of shapely Points in the mesh.
    """
    # Start out with a grid of N = density x density points. Polygon bounds:  (minx, miny, maxx, maxy)
    x = np.linspace(area.total_bounds[0], area.total_bounds[2], density)
    y = np.linspace(area.total_bounds[1], area.total_bounds[3], density)
    x_ind, y_ind = np.meshgrid(x, y)

    # Create an array of mesh points. np.ravel linearly indexes an array into a row.
    mesh_points = np.array([np.ravel(x_ind), np.ravel(y_ind)]).T

    # Convert coordinate tuples into shapely points.
    mesh_points = [Point(point[0], point[1]) if not altitude
                   else Point(point[0], point[1], altitude) for point in mesh_points]

    return mesh_points


def create_overlapping_mesh(area: gpd.GeoDataFrame, spacing: int = 1,
                            mesh_size: int = 25) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Create a mesh of polygons as close to size mesh_size x mesh_size as possible over a specific area.

    Parameters
    ----------
    area : gpd.GeoDataFrame
        The area to cover with the mesh. CRS should be a geographic coordinate system that uses D.d.
    spacing : int, default 1 km
        Distance apart receiver points should be in kilometers
    mesh_size : int, default 25 km
        The target size in kilometers of a mesh square (mesh_size x mesh_size)

    Returns
    -------
    An overlapping mesh of squares that cover the requested area.
    A GeoDataFrame of the center points used to create the mesh squares.
    """
    equal_area_crs = coords_to_utm(area.centroid.iat[0].y, area.centroid.iat[0].x)
    area_m = area.to_crs(equal_area_crs)

    minx, miny, maxx, maxy = area_m.total_bounds
    x = np.linspace(minx, maxx, math.ceil((maxx-minx)/(spacing*1000)))
    y = np.linspace(miny, maxy, math.ceil((maxy-miny)/(spacing*1000)))
    x_ind, y_ind = np.meshgrid(x, y)

    # np.ravel linearly indexes an array into a row.
    mesh_points = [Point(point[0], point[1]) for point in np.array([np.ravel(x_ind), np.ravel(y_ind)]).T]
    mesh_points = gpd.GeoDataFrame({'geometry': mesh_points}, geometry='geometry', crs=equal_area_crs)

    # Only keep points that fall within the study area.
    mesh_points = gpd.sjoin(mesh_points, area_m, op='within')[['geometry']]

    # Create mesh around points.
    mesh = mesh_points.buffer(mesh_size*1000, cap_style=3)

    mesh.reset_index(drop=True, inplace=True)
    mesh_points.reset_index(drop=True, inplace=True)

    return mesh.to_crs(area.crs), mesh_points.to_crs(area.crs)


def project_raster(input_raster: str, output_raster: str, crs: str):
    """
    Project a raster to a new crs

    Parameters
    ----------
    input_raster : str
        Absolute path to the raster to project.
    output_raster : str
        Absolute path to where the projected raster should be written.
    crs : crs
        The CRS to project the input raster to. Of the format: 'epsg:XXXX...'
    """
    gdal.Warp(output_raster, input_raster, dstSRS=crs)


def ambience_from_raster(ambience_src: str, mic: 'Microphone') -> float:
    """
    Select the ambience level from a broadband raster at a specific microphone location.

    Parameters
    ----------
    ambience_src : str
        The absolute file path to a raster of broadband ambience.
    mic : Microphone
        A Microphone object whose location to select the broadband ambience for.

    Returns
    -------
    Lx : float
        The ambience level at the microphone location.
    """
    with rasterio.open(ambience_src) as raster:
        projected_mic = mic.to_crs(raster.crs)
        band1 = raster.read(1)
        Lx = band1[raster.index(projected_mic.x, projected_mic.y)]

    return Lx


def ambience_from_nvspl(ambience_src: 'Nvspl', quantile: int = 50,
                        low_hz: float = "12.5", high_hz: float = "20000", broadband: bool = False):
    """

    Parameters
    ----------
    ambience_src : Nvspl
        An NVSPL object to calculate ambience from.
    quantile : int, default 50
        This quantile of the data will be used to calculate the ambience.
    low_hz : float, default 12.5
        Lowest 1/3 octave band to include.
    high_hz : float, default 20000
        Highest 1/3 octave band to include.
    broadband : bool, default False
        If True, quantiles will be calculated from the dBA column instead of the 1/3rd octave band columns.

    Returns
    -------
    Lx
    """
    # filter out high wind periods
    low_wind = ambience_src.loc[ambience_src["WindSpeed"] <= 5.0, :]

    if broadband:
        Lx = low_wind.loc[:, 'dbA'].quantile(1 - (quantile / 100))
    else:
        Lx = low_wind.loc[:, low_hz:high_hz].quantile(1 - (quantile / 100))

    return Lx


def compute_fbeta(valid_points: gpd.GeoDataFrame, active_space: gpd.GeoDataFrame,
                  beta: float = 1.0) -> Tuple[float, float, float, int]:
    """
    Given a set of annotated points and an active space geometry, compute accuracy metrics such as F1 score, precision,
    and recall.
        TP = True Positives
        FP = False Positives
        FN = False Negatives

    Parameters
    ------
    valid_points : gpd.GeoDataFrame
        Annotated points. Must include geometry and an 'audible' column.
    active_space : gpd.GeoDataFrame
        Polygon or Multipolygon of a computed active space.
    beta : float, default 1.0
        Beta value to use when calculating fbeta

    Returns
    -------
    fbeta : float
        fbeta score (more here: https://en.wikipedia.org/wiki/F-score)
    precision : float
        Defined TP/(TP+FP), measure of how well a positive test corresponds with an actual audible flight
    recall : float
        Defined TP/(TP+FN), measure of how well an audible flight is marked as audible by the given active space
    n_tot: int
        number of points annotated.
    """
    # Before computing anything, make sure projections match:
    if valid_points.crs != active_space.crs:
        valid_points.to_crs(active_space.crs, inplace=True)

    # iterate through all valid points and check if they are in the active space... this takes a while.
    in_AS_gdf = gpd.clip(valid_points, active_space)

    # make an `in_activespace` column and set to true for points inside mask
    valid_points['in_AS'] = False
    valid_points.loc[in_AS_gdf.index, 'in_AS'] = True

    in_AS = valid_points.in_AS.values  # convert both of these columns to boolean arrays for easier
    audible = valid_points.audible.values

    # compute true positives, etc.
    TP = np.all([in_AS, audible], axis=0).sum()
    FP = np.all([in_AS, ~audible], axis=0).sum()
    FN = np.all([~in_AS, audible], axis=0).sum()
    n_tot = len(valid_points)

    precision = TP / (TP + FP)  # specificity... if a flight enters the active space, is it actually audible?
    recall = TP / (TP + FN)  # sensitivity... if a flight is audible, does it enter the active space?
    fbeta = (1 + np.power(beta, 2)) * ((precision * recall) / ((np.power(beta, 2) * precision) + recall))

    return fbeta, precision, recall, n_tot

def contiguous_regions(condition):

    """
    Finds contiguous True regions of an input boolean array. 
    
    Parameters
    ----------
    condition : `np.ndarray` of dtype "bool" 
                 or `np.ndarray` and conditional logic statement to produce the boolean array 
                 e.g., arr    or    arr > 5.5

    Returns
    -------
    idx : 2-D numpy.ndarray
        A 2-D int array where the first column is the start index of each contiguous True region, 
        and the second column is the end index of each contiguous True region.

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

def audibility_to_interval(aud, invert=False):

    '''
    Given an audibility time series in 1-D binary format (e.g., detection/non-detection sequence)
    separate it into two, 2-D arrays of {begin, end} interval pairs. The first represents the temporal 
    bounds of each audible noise event, the second, each noise-free interval.
    
    Noise intervals are closed via observation, but to close noise-free intervals 
    we must account for the beginning and end of the observation record. Use of closed intervals 
    ensures that no index is considered to be both noise and not-noise.

    Parameters
    ----------
    aud : 1-D array-like
        A 1-D boolean array representing audibility states of detection (e.g., True or 1) and non-detection (False or 0).
    invert : bool
        Detection and non-detection are mapped to 0 and 1, respectively. If True, invert the mapping. Default: False.

    Returns
    -------
    noise_intervals: 2-D numpy.ndarray
        A 2-D int array of closed intervals bounding audible noise events.  
        The first value in the pair is the start index, the second value is the end index.
    noise_free_intervals: 2-D numpy.ndarray
        A 2-D int array bounding closed noise-free intervals.  
        The first value in the pair is the start index, the second value is the end index.
    '''
    aud = aud.astype('bool')
    if invert == True:
        aud = np.invert(aud) # invert detection mappings
    
    # compute naiive intervals
    noise_intervals = contiguous_regions(aud == True)
    noise_free_intervals_naiive = contiguous_regions(aud == False)
    
    nfi_starts = noise_free_intervals_naiive.T[0]
    nfi_ends = noise_free_intervals_naiive.T[1]

    # the record begins with noise...
    if(noise_intervals[0, 0] == 0):
        # ...the first noise free interval (and thus ALL intervals) 
        #    need to start one second later
        nfi_starts = nfi_starts + 1
    
    # the record begins with quietude...
    else:
        # ...the first noise free interval stays the same, and equals zero
        # the rest are + 1
        nfi_starts = nfi_starts + 1
        nfi_starts[0] = 0

    # the record ends with noise...
    if(noise_intervals[-1, 0] == 0):
        # ...the last noise free interval (and thus ALL intervals) need to end one second earlier
        nfi_ends = nfi_ends - 1
    
    # the record ends with quietude...
    else:
        # ...the last noise free interval stays the same, and equals zero
        # the rest are - 1
        save = nfi_ends[-1]
        print(save)
        nfi_ends = nfi_ends - 1
        nfi_ends[-1] = save

    # recompose NFIs using updated, correct values
    noise_free_intervals = np.array([nfi_starts, nfi_ends]).T
    
    return noise_intervals, noise_free_intervals

def calculate_duration_summary(noise_intervals):

    '''
    Compute durations from interval-based noise event data. 
    Summarize the central tendency and variability using parametric and non-parametric estimators.

    Parameters
    ----------
    noise_intervals: 2-D numpy.ndarray
        A 2-D int array where the first column is the start index of each contiguous True region, 
        and the second column is the end index of each contiguous True region.

    Returns
    -------
    duration_summary : tuple
        A tuple of duration-based acoustic metrics:
            idx [0] a list of each event's duration
            idx [1] the mean duration
            idx [2] the standard deviation of the durations
            idx [3] the median duration
            idx [4] the median absolute deviation of the durations
    '''

    # the durations, themselves, are found by differencing (end - begin)
    duration_list = noise_intervals.T[1] - noise_intervals.T[0]

    # mean duration
    mean = np.mean(duration_list)

    # standard deviation duration
    stdev = np.std(duration_list)

    # median duration
    median = np.percentile(duration_list, 50)

    # median absolute deviation of duration
    mad = np.percentile(np.absolute(duration_list - median), 50)

    # combine the results into a single array
    duration_summary = (duration_list, mean, stdev, median, mad)

    return duration_summary

def select_optimal(unit: str, site: str, year: int,
                   valid_points, active_space_polygons: list, beta_=1.0,
                   verbose=True, plot=True, plot_savepath=None):
    """
    From a ground-truthed causal dataset and a set of active space polygons, 
    select the optimal geospatial prediction of observed audibility.

    Parameters
    ------
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA', '009'
    year : int
        Deployment year. YYYY
    valid_points : gpd.GeoDataFrame
        Annotated points. Must include geometry and an 'audible' column.
    active_space_polygons : list of gpd.GeoDataFrame objects
        A list of 2-tuples. The first value in each tuple is a gain string (ex. "O_+055" or "O_-125").
        The second value is a single-row `gpd.GeoDataFrame` of Polygons/Multipolygons representing an active space prediction.
        The entire list usually represents the geometries computed over a range of gains.
    verbose : bool, default True
        If True, print out F-score, precision, and recall for each gain.
    plot : bool, default True
        If True, plot the precision recall curve.
    plot_savepath : str
        Filename for saving the precision recall plot. If provided, the precision recall curve will be saved instead of shown. Requires plot=True to work.
    beta_ : float, default 1.0
        Beta value to use when calculating F-Beta

    Returns
    -------
    best_omni : float
        The optimal gain, corresponding to a maximum F-Beta
    max_fbeta : float
        The maximum F-Beta score at the optimal gain
    best_precision: float
        The precision component of the maximum F-Beta score at the optimal gain
    best_recall: float
        The recall component of the maximum F-Beta score at the optimal gain
    detection_results: pd.DataFrame
        A table of the F-Beta, precision, and recall values indexed by gain. The `.name` attribute contains deployment and F-Beta information.
    """

    # map symbolic sign onto a numeric equivalent
    sign = {"+":1, "-":-1}

    # intialize a table to write each test result
    detection_results = pd.DataFrame([])
    detection_results.name = f"{unit}{site}{year} {beta_}"
    
    # calculate the F-Beta score to determine the active space with optimal detection properties
    precisions = []
    recalls = []
    f_betas = []
    max_fbeta = 0
    best_omni = None
    best_recall = None
    best_precision = None
    for omni, res in active_space_polygons:

        # it is convenient to store the gain in a simple, numeric representation
        numeric_gain = sign[omni[-4:-3]]*int(omni[-3:])/10

        # compute the f-measure for the geometry at a given Beta
        fbeta, precision, recall, n_tot = compute_fbeta(valid_points, res, beta_)

        # store the overall result and its subcomponents in a `pd.DataFrame` object
        detection_results.loc[omni, f"F-{beta_}"] = fbeta
        detection_results.loc[omni, "Recall"] = recall
        detection_results.loc[omni, "Precision"] = precision
        detection_results.loc[omni, "gain"] = numeric_gain

        if verbose:
            print(f"omni: {omni} --> F-{beta_}: {fbeta:0.3f} precision: {precision:0.3f} recall: {recall:0.3f}")
        
        if fbeta > max_fbeta:
            max_fbeta = fbeta
            best_omni = omni
            best_recall = recall
            best_precision = precision

        detection_results.sort_values("gain", ascending=True, inplace=True) # it makes sense to plot (and return) the object sorted by gain
    
    if plot:
        # create Precision-Recall Plot.
        fig, ax = plt.subplots()
        ax.plot(detection_results["Recall"], detection_results["Precision"], ls="-", lw=0.2, marker="o", ms=2, color="k")
        ax.plot(best_recall, best_precision, ls="", marker="o", ms=5, color="None", markeredgecolor="limegreen")
        ax.set_title(f'Precision-Recall Curve, {beta_:.1f}', loc="left")
        ax.set_ylabel('Precision')
        ax.set_xlabel('Recall')

        if plot_savepath is None:
            plt.show()
        else:
            plt.savefig(plot_savepath)

    print(f"The best performing omni source for F-{beta_} is: {best_omni} (fbeta: {max_fbeta:0.3f})")
    return best_omni, max_fbeta, best_precision, best_recall, detection_results


def normalize_point_density(points: gpd.GeoDataFrame, study_area: gpd.GeoDataFrame, bandwidth: float=100.0, cellsize: float=100.0, random_seed=None, visualize=False):
    """
    Drop points from a dataframe that are in very dense areas, to normalize the point density.
    This uses kernel density estimation to get a density for each point, specifically FFTKDE from KDEpy,
    which is a very fast algorithm. This speed is critical for processing as many points as we do.
    
    Parameters
    ----------
    points: gpd.GeoDataFrame
        Points to process
    study_area: gpd.GeoDataFrame
        Study area polygon. Used to determine the UTM zone, and any points outside will be removed.
    bandwidth: float, default 100.0
        Standard deviation of the gaussian kernel, in meters
    cellsize: float, default: 100.0
        Resolution of the grid used to evaluate point density. Points are assigned the density at the nearest grid point.
    visualize: bool, default False
        If true, show a plot of the relative point densities before and after.

    Returns
    -------
    points: gpd.GeoDataFrame
        Copy of the input points with some points removed.
    """
    # convert to UTM because we are doing distance-based processing
    orig_crs = points.crs  # remember this so we can project back later
    crs = NMSIM_bbox_utm(study_area)
    points = points.to_crs(crs)
    
    # clip to study area to reduce the amount of grid cells needed, since the KDEpy grid needs to include all points
    # also those points can be safely excluded from consideration when fitting the active space
    points = points.clip(study_area.to_crs(points.crs))

    positions = np.stack([points.geometry.x, points.geometry.y], axis=1)

    # construct grid - FFTKDE requires evaluating densities on an evenly spaced grid
    # the grid needs to include all the points with some padding
    pad = 3 * bandwidth  # enough to be outside of the range of the gaussian kernel
    xmin, ymin = positions.min(axis=0)
    xmax, ymax = positions.max(axis=0)
    # use add cellsize to np.arange() end value since end value is exclusive, but we want it to be inclusive
    x_coords = np.arange(xmin-pad, xmax+pad+cellsize, cellsize)
    y_coords = np.arange(ymin-pad, ymax+pad+cellsize, cellsize)
    mgrid = np.meshgrid(x_coords, y_coords)
    # put x and y coords for each point in the grid in an array of shape (n_grid_points, 2)
    # note that we use the transpose of the meshgrid so that the order is as KDEpy expects,
    # where x starts at the first value while y sweeps through, then x goes to the next value and y sweeps through
    grid = np.stack([mgrid[0].T.ravel(), mgrid[1].T.ravel()], axis=1)

    # run kernel density estimation
    grid_densities = FFTKDE(bw=bandwidth).fit(positions).evaluate(grid)

    # assign points the density of their nearest grid point
    # use a KDTree for faster nearest-neighbor querying
    tree = KDTree(grid)
    distance, idx = tree.query(positions)
    point_densities = grid_densities[idx]

    # drop points to reach a certain target density
    # empirically and intuitively, the median density is a good target
    # so for points with higher densities than the median, keep each point
    # with probability median_density/point_density
    n_before = len(points)
    median = np.median(point_densities)
    p_keep = median / np.maximum(point_densities, median)  # = 1 for points with density <= median
    rng = np.random.default_rng(seed=random_seed)
    keep = rng.random(len(points)) <= p_keep
    points = points[keep]
    print(f"Went from {n_before} to {len(points)} points when normalizing point density")

    if visualize:
        # run it again to visualize the difference
        new_pos = np.stack([points.geometry.x, points.geometry.y], axis=1)
        new_grid_densities = FFTKDE(bw=bandwidth).fit(new_pos).evaluate(grid)

        # reshape densities to work with matplotlib
        z0 = grid_densities.reshape(x_coords.size, y_coords.size).T
        z1 = new_grid_densities.reshape(x_coords.size, y_coords.size).T

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))
        ax[0].set_title("Relative Density Before")
        ax[0].contourf(mgrid[0], mgrid[1], z0, levels=100)
        ax[0].axis("equal")
        ax[1].set_title("Relative Density After")
        ax[1].contourf(mgrid[0], mgrid[1], z1, levels=100)
        ax[1].axis("equal")
        fig.tight_layout()
        plt.show()

    # convert points back to the original crs
    points = points.to_crs(orig_crs)

    return points