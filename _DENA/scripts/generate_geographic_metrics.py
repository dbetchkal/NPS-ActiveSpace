import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import stats
from scipy.spatial.distance import directed_hausdorff, cdist
from scipy.signal import find_peaks
from shapely.geometry import Point, LineString
from tqdm import tqdm
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.patches as patches

__all__ = [
    'clip_events_to_time_period',
    'tracks2events',
    'get_all_stats',
    'calculate_spatial_stats',
    'plot_events',
    'circular_sliding_avg',
    'find_circular_peaks',
    'endpoints_around_active',
    'identify_stereotypical_tracks'
]

## ========================================== STATISTICS/METRICS ======================================== ##

def clip_events_to_time_period(df, start_col, end_col, start_date, end_date, months=list(range(1,13))):
    df = df.copy()
    during_time_period = (df[end_col] > start_date) & (df[start_col] < end_date)
    during_correct_months = df[start_col].dt.month.isin(months) | df[end_col].dt.month.isin(months)
    df = df[during_time_period & during_correct_months]
    # shorten events that partly exist outside of the time period
    df[start_col] = np.maximum(df[start_col].values, start_date)
    df[end_col] = np.minimum(df[end_col].values, end_date)
    return df


def tracks2events(tracks, start_date, end_date, min_dur=30):
    '''
    Performs audibility binarization on outputs of the `NPS-ActiveSpace.audible_transits` module. 
    
    This function collapses a set of audible transits spanning a certain time period into 
    unique noise events, which may contain multiple overlapping transits.
    
    The two output event dataframes capture alternating periods of time: 
        (0)  Noise Events (e.g., 'Noisy intervals')
        (1)  Noise-free intervals (NFI)
    
    Parameters
    ----------
    tracks : GeoDataFrame
        A GeoDataFrame containing the fully interpolated, cleaned, clipped, and extrapolated tracks
        resembling those produced by `NPS-ActiveSpace.audible_transits`. Only requires columns 'entry_time' and 'exit_time'
    start_date : string
        The initial date of tracks to include, formatted as 'yyyy-mm-dd'. Note that midnight at the beginning of this day should fall within the monitoring period (not before).
    end_date : string
        The last date of tracks to include, formatted as 'yyyy-mm-dd'. Note that midnight at the beginning of this day should fall within the monitoring period (not after).
    min_dur : float, default 30
        The minimum event duration to include, in seconds
         
        
    Returns
    -------
    event_df : pd.DataFrame
        A DataFrame containing each noise event. Looks like:
            start_time | end_time | duration
    NFI_df : pd.DataFrame
        A DataFrame containing each noise-free interval. Looks like:
            start_time | end_time | duration
    TA : float
        Total fraction of time audible. Ranges from 0 to 1.
    ''' 

    print("Combining audible transits into a binary event time series.")
    
    buffer = np.timedelta64(30, 's') # Max number of seconds between tracks in order to combine
    
    start_date = np.datetime64(start_date)  # conversion to datetime64
    end_date = np.datetime64(end_date)      # conversion to datetime64

    tracks = clip_events_to_time_period(tracks, "entry_time", "exit_time", start_date, end_date)

    if tracks.empty:
        event_df = pd.DataFrame(columns=["start_time", "end_time", "duration"])
        NFI_df = pd.DataFrame(columns=["start_time", "end_time", "duration"])
        return event_df, NFI_df, 0.0
    
    tracks.sort_values(by=['entry_time'], inplace=True)
    entry_times = np.asarray(tracks.entry_time) # datetime format
    exit_times = np.asarray(tracks.exit_time)   # datetime format
    elapsed_times = exit_times - entry_times    # timedelta64 format
    
    # Make copies of the entry and exit times arrays for calculation, 
    # we may want to use the originals later
    entry_times_cp = entry_times.copy()
    exit_times_cp = exit_times.copy()
    
    # Loop through each entry time and look forward to find the corresponding event end time 
    # (not necessarily the same as the exit time, if there are overlapping flights)
    i = 0      # increment for entry times
    while i < len(entry_times_cp)-1:
        j = 0       # increment for exit times

        # set event end time to the first track's exit time (overwrite if needed in while loop)
        event_end = exit_times_cp[i]  

        # check if the current track exits after the next track starts; 
        # if so, adjust the event end time and increment
        while event_end >= (entry_times_cp[i+j+1] - buffer):
            event_end = max(exit_times_cp[i+j+1], event_end)  # account for case where first track starts first but also ends later
            j += 1
            if (i+j+1) >= len(entry_times_cp):
                break  # NEED to leave while loop if we are about to exceed the size our entry_times array -- !can't look forward beyond the last track!
        # entry time stays the same, set new exit time
        exit_times_cp[i] = event_end

        # Get rid of all the tracks that were blobbed together with the first one (e.g., i). 
        # the 'j' increment tells us how many tracks we combined
        exit_times_cp = np.delete(exit_times_cp, slice(i+1, i+j+1))
        entry_times_cp = np.delete(entry_times_cp, slice(i+1, i+j+1))
        i += 1

    # Filter out short events and rename these arrays; they now represent the start and end times of each combined event
    event_start_times = entry_times_cp[exit_times_cp - entry_times_cp >= np.timedelta64(min_dur, 's')]
    event_end_times = exit_times_cp[exit_times_cp - entry_times_cp >= np.timedelta64(min_dur, 's')]
    
    # Account for event-less time at beginning and end of timeframe in question
    # This is needed for the inaudible_begins and inaudible_ends indices to match up properly so we can subtract to compute durations.
    inaudible_begins = np.insert(event_end_times, 0, start_date)
    inaudible_ends = np.append(event_start_times, end_date)
    # Filter out zero-duration NFIs caused by the timeframe starting or ending with an event
    zero_duration = inaudible_begins == inaudible_ends
    inaudible_begins = inaudible_begins[~zero_duration]
    inaudible_ends = inaudible_ends[~zero_duration]

    # Get durations for events and NFIs
    audible_times = event_end_times - event_start_times        # List of event durations in timedelta64 format
    inaudible_times = inaudible_ends - inaudible_begins            # List of noise-free interval durations in timedelta64 format

    # Calculate time audible
    total_audible = sum(audible_times)/np.timedelta64(1,'s')   # Add up all event durations to get total audible time (convert to float in seconds)
    total_time = (end_date - start_date)/np.timedelta64(1,'s')     # Calculate total time of timeframe in seconds 
    TA = 100 * total_audible / total_time                # (%) of total time with an audible event (Time Audible)

    # We organize this information into two dataframes -> Noise events, Noise-free intervals
    duration_list = audible_times/np.timedelta64(1,'s')  # list of event duration lengths in seconds
    NFI_list = inaudible_times/np.timedelta64(1,'s')     # list of noise-free interval lengths in seconds
    event_df = pd.DataFrame(data={'start_time': event_start_times, 'end_time': event_end_times, 'duration': duration_list})
    NFI_df = pd.DataFrame(data={'start_time': inaudible_begins, 'end_time': inaudible_ends, 'duration': NFI_list})

    return event_df, NFI_df, TA


def _split_events(df, freq):
    """Split events that span hour/day/etc boundaries into consecutive events.
    
    Parameters
    ----------
    df: pd.DataFrame
        A DataFrame containing each event. Looks like:
            start_time | end_time | duration
    freq: str
        A string representing a time frequency in pandas, e.g. 'h', 'd', etc.
    
    Returns
    -------
    df: pd.DataFrame
        A new DataFrame containing split events.
    """

    split_rows = []
    for _, row in df.iterrows():
        current_start = row['start_time']
        while current_start < row['end_time']:
            # Find the end of the current frequency bucket
            next_boundary = current_start.floor(freq) + pd.to_timedelta(1, unit=freq)
            current_end = min(row['end_time'], next_boundary)
            # Add segment
            split_duration = (current_end - current_start).total_seconds()
            split_rows.append({
                'start_time': current_start,
                'end_time': current_end,
                'duration': split_duration
            })
            current_start = current_end

    return pd.DataFrame(split_rows)


def _time_binned_df(event_df, start_date, end_date, months, freq):
    # determine appropriate human-readable time period prefix
    prefix_map = {"h": "hourly", "d": "daily", "w": "weekly", "m": "monthly"}
    prefix = f"{prefix_map[freq]}_" if freq in prefix_map else ""

    # prepare a dataframe to hold values for each time period, indexed by the time at the start of that period
    # make sure that the entire time range of interest is represented, so that periods without data are accounted for
    date_index = pd.date_range(start_date, end_date, freq=freq, inclusive="left")
    date_index = date_index[date_index.month.isin(months)]
    periods_df = pd.DataFrame(index=date_index, columns=[prefix+"time_audible", prefix+"event_count"])

    # group by each time period and calculate the number of events STARTING in each period
    # (which is different than the number of events overlapping that period)
    for period, group in event_df.groupby(event_df["start_time"].dt.floor(freq)):
        periods_df.loc[period, prefix+"event_count"] = len(group)  # event count is just length of dataframe group

    # group by each time period and calculate time audible
    # first, split events overlapping multiple periods, so that the relevant chunks of each event
    # are assigned to the correct time period
    split = _split_events(event_df, freq)
    period_duration = pd.to_timedelta(1, freq).total_seconds()
    for period, group in split.groupby(split["start_time"].dt.floor(freq)):
        periods_df.loc[period, prefix+"time_audible"] = 100 * group["duration"].sum() / period_duration

    periods_df.fillna(0, inplace=True)  # Fill in any periods that had no audible events with count=0 and % time audible=0

    return periods_df


# nice trick for aggregating quantiles cleanly, inspired from here:
# https://skeptric.com/pandas-aggregate-quantile/
class Quantile:
    def __init__(self, q):
        self.q = q
        self.__name__ = f"{int(100*q)}% quantile"
    
    def __call__(self, series: pd.Series):
        return series.quantile(self.q)


def get_all_stats(event_df, NFI_df, start_date, end_date, months=list(range(1,13)), quantiles=.5):
    start_date = np.datetime64(start_date)  # Convert to datetime64
    end_date = np.datetime64(end_date)      # Convert to datetime64

    # Input validation. Both 'quantiles' and 'months' paramters must be converted to lists
    quantiles = [quantiles] if type(quantiles)!=type([]) else quantiles
    months = [months] if type(months)!=type([]) else months

    # Make sure months are between 1 and 12
    for month in months:
        if (month < 1) | (month > 12):
            print("Warning: Invalid months. Must be a list of integers from 1-12. Ignoring months parameter...")
            months=list(range(1,13))

    event_df = clip_events_to_time_period(event_df, "start_time", "end_time", start_date, end_date, months)
    NFI_df = clip_events_to_time_period(NFI_df, "start_time", "end_time", start_date, end_date, months)

    # prepare the values we want statistics for
    values = {
        "event_duration": event_df["duration"],
        "NFI_duration": NFI_df["duration"]
    }
    # include time audible and event count, binned by hour and by day
    for freq in ['d', 'h']:
        binned_df = _time_binned_df(event_df, start_date, end_date, months, freq)
        for col in binned_df.columns:
            values[col] = binned_df[col]
    
    # prepare the statistics we want
    agg_stats = ['min', 'max', 'mean', 'std', stats.median_abs_deviation, 'sem'] + [Quantile(q) for q in quantiles]

    # compute statistics
    statistics = {}
    for col, series in values.items():
        statistics[col] = series.agg(agg_stats)
        # convert duration fields to timedeltas, more human readable
        if col in ["event_duration", "NFI_duration"]:
            statistics[col] = pd.to_timedelta(statistics[col], unit="s")
    
    return pd.DataFrame(statistics)


def calculate_spatial_stats(tracks, active):
    '''
    Calculates spatial statistics on audible transits through a given active space. 
    On an event-wise basis, the displacement distance that could render the track inaudible:
                                                                            - mean
                                                                            - maximum
    In aggregate:
                                                                             - min
                                                                             - max
                                                                             - mean
                                                                             - median

                                                    
    Parameters
    ----------
    tracks : gpd.GeoDataFrame
        GeoDataFrame containing the audible transit flight tracks. Need to have geometry column.
    active : gpd.GeoDataFrame
        GeoDataFrame containing the active space geometries.

    Returns
    -------
    distance_from_audibility_stats : gpd.GeoDataFrame
        GeoDataFrame of overall stats containing the min, max, mean, and median for:
            mean distance from inaudibility and max distance from inaudibility

    '''
    boundary = np.asarray(active.iloc[0].exterior.coords)
    
    # Calculate distance between each line point and active space
    distances = tracks.geometry.apply([lambda line: cdist(np.asarray(line.coords)[:,:2], boundary).min(axis=1)])
    
    # Calculate mean of these distances for each track
    tracks['mean_distance_from_inaudibility'] = distances.apply([lambda dists: dists.mean()])
    
    # Calculate max of these distances for each track
    tracks['max_distance_from_inaudibility'] = distances.apply([lambda dists: dists.max()])
    
    # Create statistics dataframe consisting of min, max, mean, and median for each
    distance_from_inaudibility_stats = tracks.agg({'max_distance_from_inaudibility':['min', 'max', 'mean', 'median'], 'mean_distance_from_inaudibility':['min', 'max', 'mean', 'median']})
    
    return distance_from_inaudibility_stats


## ========================================== VISUALIZATION =============================================== ##

def _split_interval_by_hour(start, end):
    """Yield (segment_start, segment_end) tuples split at each hour boundary."""
    while start < end:
        next_hour = (start + timedelta(hours=1)).replace(minute=0, second=0, microsecond=0)
        segment_end = min(end, next_hour)
        yield (start, segment_end)
        start = segment_end
        

def plot_events(start_times, end_times, title="Events", labels=False):
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    earliest_hour = start_times.min().floor("h")
    latest_hour = end_times.max().ceil("h")

    for start, end in zip(start_times, end_times):
        for s, e in _split_interval_by_hour(start, end):
            x_start = (s.minute + s.second / 60)
            duration = (e - s).total_seconds() / 60
            y = s.floor('h').timestamp()

            # Plot rectangle
            ax.add_patch(patches.Rectangle((x_start, y), duration, 0.8 * 3600,
                                      linewidth=1, edgecolor='black', facecolor='skyblue', alpha=0.5))
            # Plot Text
            if labels:
                ax.text(x_start, y,
                        f"{start.strftime('%H:%M')}-{end.strftime('%H:%M')}",
                        ha='left', va='top', fontsize=8)

    # Configure axes
    hours = pd.date_range(earliest_hour, latest_hour, freq="h")
    ax.set_xlim(0, 60)
    ax.set_ylim(latest_hour.timestamp(), earliest_hour.timestamp())  # inverted y axis so time flows down
    ax.set_yticks(hours.astype(int) // 1e9)
    yticklabels = hours.map(lambda t: t.strftime("%m-%d  %H:00") if t.hour == 0 else t.strftime("%H:00")).tolist()
    yticklabels[0] = hours[0].strftime("%m-%d  %H:00")  # make sure the first hour is labeled with the date
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('Minutes within the hour')
    ax.set_ylabel('Hour')
    ax.set_title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

## ========================================== STEREOTYPICAL TRACKS ======================================== ##

def circular_sliding_avg(vector, window_len):
    '''
    Smooth circular data by calculating sliding averages on a circular array. 
    Used to create an integrated map of entry/exit points around the boundary of an active space.

    Parameters
    ----------
    vector : numpy arrary
        Array that contains the values to be calculated on
    window_len: int (odd, positive)
        Length of sliding average window. Creates window, centered around each sample, of np.ones(window_len). In other words, equally weighted.

    Returns
    -------
    smoothed : numpy array
        Array of same length as input vector, contains respective sliding average values.
    
    '''
    # prepend end half-window of vector and append start half-window of vector
    vector_wrapped = np.concatenate((vector[-(window_len//2):], vector, vector[:window_len//2]))

    # Convolution using mode 'valid' will only output frames where an entire window can fit
    smoothed = np.convolve(vector_wrapped, np.ones(window_len), mode='valid')/window_len
    
    return smoothed

def find_circular_peaks(column, distance_delta, peak_distance):
    '''
    A specific adaptation of SciPy.Signal's 'find_peaks' algorithm for use on an active space polygon.
    Critical to implementation is the ability to wrap the input array. In this case, we use a minimum-distance
    criterion to separate peaks on a polygon (e.g., values on a closed LineString).
    
    This is important if you wish to only use peaks separated by at least a certain distance for a circular array, 
    such as values on a closed LineString.

    Parameters
    ----------
    column : gpd.GeoSeries
        Column from a GeoDataFrame to find circular peaks from. Each value corresponds with a segment of the active space boundary.
    distance_delta : int (m)
        The length of the segments that the active space boundary has been broken up into.
    peak_distance : int (m)
        Minimum distance between peaks, in meters.

    Returns
    -------
    sorted_peak_indices : numpy array of ints
        Array containing the indices of the peaks, sorted by descending peak height.
    
    '''
    samples_between_peaks = peak_distance//distance_delta 
    entries_vector = column.to_list()
    
    # Add bits of the end to the beginning and vice versa to give the impression of a circular array
    entries_wrapped = np.concatenate((entries_vector[-samples_between_peaks:], entries_vector, entries_vector[:samples_between_peaks]))

    # Use SciPy's find_peaks to identify peaks that are above the 80% quantile and meet the distance between peaks requirement
    peaks, properties = find_peaks(entries_wrapped, height=column.quantile(.8), distance=samples_between_peaks)

    # Sort indices by peak_height
    sort_indices = properties["peak_heights"].argsort()
    sorted_peaks = peaks[sort_indices[::-1]] # reverse sort

    # Adjust indices to undo the wrapping of the indices; 
    # get rid of peaks in the prepended AND appended parts of the wrapped vector.
    sorted_peak_indices = sorted_peaks[(sorted_peaks >= samples_between_peaks) & (sorted_peaks < len(entries_vector)+samples_between_peaks)] - samples_between_peaks
    return sorted_peak_indices

def endpoints_around_active(active, tracks, distance_delta, peak_distance, endpoint_type):
    '''
    Function to calculate the position of the most frequent entry and exit points through an active space.
    (e.g., transportation fluxes along the polygon). A sliding average is applied to spatially smooth
    Uses a sliding average on a segmented active space to generate smoothed spatial data, 
    while also detecting the peaks of this data.

    Parameters
    ----------
    active : gpd.GeoDataFrame
        GeoDataFrame containing the active space geometries.
    tracks : gpd.GeoDataFrame
        GeoDataFrame containing the audible transit flight tracks. Need to have entry_position and exit_position columns.
    distance_delta : int (m)
        The length of the segments that the active space boundary has been broken up into.
    peak_distance : int (m)
        Minimum distance between peaks, in meters.
    endpoint_type : str
        Either 'entry' or 'exit', depending on what we're searching for

    Returns
    -------
    active_gdf : gpd.GeoDataFrame
        GeoDataFrame of segmentized active space boundary. Segmented into n distance_delta length pieces, 
        where n is boundary.length/distance_delta.
        Looks like:
            midpoints  |  segments  |  entries  | exits
    peak_indices : numpy array of ints
        Array containing the indices of the peaks, sorted by descending peak height. 
    '''   
    line = active.exterior.iloc[0]
    # generate the equidistant points
    distances = np.arange(0, line.length, distance_delta)
    
    pts = [line.interpolate(distance) for distance in distances] # perform the interpolation
    if len(pts)%2 == 0 : pts = pts + [Point(line.coords[-1])]
    line_pts = pts[::2]
    midpoints = pts[1::2]
    
    lines=[]
    for i in range(len(line_pts)-1):
        lines.append(LineString(line_pts[i:i+2]))
    
    active_gdf = gpd.GeoDataFrame()
    active_gdf['midpoints'] = gpd.GeoSeries(midpoints, crs=active.crs)
    active_gdf['segments'] = gpd.GeoSeries(lines, crs=active.crs)
    active_gdf.set_geometry('segments', inplace=True)
    
    entries = active_gdf.segments.apply([lambda line: tracks.entry_position.intersects(line.buffer(1)).sum()])
    exits = active_gdf.segments.apply([lambda line: tracks.exit_position.intersects(line.buffer(1)).sum()])
    active_gdf['entries'] = circular_sliding_avg(entries.values.T[0], 11)
    active_gdf['exits'] = circular_sliding_avg(exits.values.T[0], 11)
    active_gdf.set_geometry('midpoints', inplace=True)

    if endpoint_type == 'entry': 
        peak_indices = find_circular_peaks(active_gdf.entries, distance_delta, peak_distance)
    elif endpoint_type == 'exit':
        peak_indices = find_circular_peaks(active_gdf.exits, distance_delta, peak_distance)
    else:
        print("error: endpoint_type must be 'entry' or 'exit'")
        return 0
                                            
    return active_gdf, peak_indices

def identify_stereotypical_tracks(active, tracks, distance_delta=100, peak_distance=1000):
    '''
    Function to identify stereotypical tracks across an active space. 
    Picks out specific paths that closely resemble the most tracks.

    Parameters
    ----------
    active : gpd.GeoDataFrame
        GeoDataFrame containing the active space geometries.
    tracks : gpd.GeoDataFrame
        GeoDataFrame containing the audible transit flight tracks. Need to have entry_position and exit_position columns.
    distance_delta : int (m)
        Defaults to 100 meters, the length of the segments that the active space boundary has been broken up into.
    peak_distance : int (m)
        Defaults to 1000 meters, the minimum distance between peaks.

    Returns
    -------
    active_gdf : gpd.GeoDataFrame
        GeoDataFrame of segmentized active space boundary. Segmented into n distance_delta length pieces, where n is boundary.length/distance_delta.
        Looks like:
            midpoints  |  segments  |  entries  | exits
    peak_indices : numpy array of ints
        Array containing the indices of the peaks, sorted by descending peak height.
    '''
    print("Identifying stereotypical tracks...")
    print(" - finding common entries and exits")
    
    # Identify common entry points and generate segmented active space boundary
    active_gdf, peak_entry_indices = endpoints_around_active(active, tracks, distance_delta, peak_distance, endpoint_type='entry')
    
    common_endpoints_rows = []
    for peak_entry_index in peak_entry_indices:
        
        # For each commmon entry point, isolate tracks that start at near the entry point (within peak_distance)
        entry_point = active_gdf.midpoints.iloc[peak_entry_index]
        tracks_entryA = tracks[tracks.entry_position.distance(entry_point) < peak_distance]
        
        # Identify common exit points for each common entry point, generate segmented active space
        active_gdf_entryA, peak_exit_indices = endpoints_around_active(active, tracks_entryA, distance_delta, peak_distance, endpoint_type='exit')
        
        for peak_exit_index in peak_exit_indices:

            # For each entry-exit pair, create a row that contains the number of tracks that start and end near the respective entry/exit (within peak_distance)
            exit_point = active_gdf_entryA.midpoints.iloc[peak_exit_index]
            transit_length = entry_point.distance(exit_point)
            corresponding_tracks = tracks_entryA[(tracks_entryA.exit_position.distance(exit_point) < peak_distance) & (tracks_entryA.entry_position.distance(entry_point) < peak_distance)]
            common_endpoints_rows.append({"entry_point":entry_point,"exit_point":exit_point, "transit_length":transit_length, "number_of_corresponding_tracks":len(corresponding_tracks)})
            
    common_endpoints_gdf = gpd.GeoDataFrame(common_endpoints_rows, geometry='entry_point', crs = tracks.crs)

    
    # Pick out the top 10 entry/exit pairs
    if len(common_endpoints_gdf) > 10:
        #common_endpoints_gdf = common_endpoints_gdf[common_endpoints_gdf.number_of_corresponding_tracks >
        common_endpoints_gdf = common_endpoints_gdf[common_endpoints_gdf.transit_length > 2000]
    common_endpoints_gdf = common_endpoints_gdf.sort_values(by='number_of_corresponding_tracks', ascending=False)

    print(f" - of {len(common_endpoints_gdf)}, identifying possible stereotypical tracks")
    
    all_tracks = tracks.copy()
    stereotype_locs = []   # list to store gdf indices of stereotypical tracks
    bundle_count = []      # list to store number of tracks that each potential stereotypical track represents        
    
    for idx, common_endpoints in tqdm(common_endpoints_gdf.iterrows(), unit=' possible stereotypes'):

        # Identify possible stereotypical tracks, should start/end near the common entry/exit point
        buffer = peak_distance # buffer around each track for identifying other similar tracks  
        subset_pct = 100
        while subset_pct > 5:
            # This while loop helps optimize performance. Make buffer smaller until only 5% or fewer tracks make the cut
            near_entry = all_tracks.entry_position.intersects(common_endpoints.entry_point.buffer(buffer))
            near_exit = all_tracks.exit_position.intersects(common_endpoints.exit_point.buffer(buffer))
            possible_tracks = all_tracks[near_entry & near_exit]
            subset_pct = 100 * len(possible_tracks)/len(all_tracks)
            buffer *= 0.95

        # Calculate number of tracks that could be bundled up with each possible stereotypical track
        possible_tracks['bundle_count'] = possible_tracks.interp_geometry.apply([lambda track: all_tracks.within(track.buffer(1000)).sum()])
        # The maximum bundle count for each entry/exit pair will be the representative track for a given path
        print(len(possible_tracks), len(all_tracks))
        stereotype_locs.append(possible_tracks[possible_tracks.bundle_count == possible_tracks.bundle_count.max()].index[0])
        bundle_count.append(possible_tracks.bundle_count.max())

    print(" - narrowing down results")
    
    stereotypical_tracks = all_tracks.loc[stereotype_locs]
    stereotypical_tracks['represents_pct'] = np.array(bundle_count)/len(all_tracks) * 100
    
    return stereotypical_tracks
