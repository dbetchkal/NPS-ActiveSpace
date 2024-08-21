import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import stats
## ========================================== STATISTICS/METRCIS ======================================== ##

def tracks2events(tracks, start_date, end_date, min_dur=30):
    '''
    Collapses a group of tracks into individual noise 'events', which may contain multiple overlapping tracks (in time). 
    The resulting event dataframe has updated start and end times for each event, as well as the total time.
    Also outputs the inverse, AKA the noise-free intervals (NFI). Uses the start and end date to bookend the NFI dataframe
    
    Parameters
    ----------
    tracks : GeoDataFrame
        A GeoDataFrame containing processed tracks. These should be the fully interpolated, cleaned, clipped, and extrapolated tracks 
        for best results. Only requires columns 'entry_time' and 'exit_time'
    start_date : string
        The start date to begin making events out of tracks, formatted as 'yyyy-mm-dd'
     end_date : string
        The end date to stop making events out of tracks, formatted as 'yyyy-mm-dd'       
        
    Returns
    -------
    event_df : GeoDataFrame
        A GeoDataFrame containing each noise event. Looks like:
            start_time | end_time | duration
    NFI_df : GeoDataFrame
        A GeoDataFrame containing each noise-free interval. Looks like:
            start_time | end_time | duration
    '''    
    print("Combining tracks into binary event time series.")
    
    buffer = np.timedelta64(30, 's') # max number of seconds between tracks in order to combine
    
    start_date = np.datetime64(start_date)  # Convert to datetime64
    end_date = np.datetime64(end_date)      # Convert to datetime64
    
    tracks.sort_values(by=['entry_time'], inplace=True)
    entry_times = np.asarray(tracks.entry_time) # datetime format
    exit_times = np.asarray(tracks.exit_time)   # datetime format
    elapsed_times = exit_times - entry_times    # timedelta64 format
    
    # Make copies of the entry and exit times arrays for calculation, may want to use the originals later
    entry_times_cp = entry_times.copy()
    exit_times_cp = exit_times.copy()
    
    # Loop through each entry time and look forward to find the corresponding event end time 
    # (not necessarily the same as the exit time, if there are overlapping flights)
    i = 0      # increment for entry times
    while i < len(entry_times_cp)-1:
        j = 0       # increment for exit times
        event_end = exit_times_cp[i]  # Set event end time to the first track's exit time (overwrite if needed in while loop)
        # Check if the current track exits after the next track starts; if so, adjust the event end time and increment
        while event_end >= (entry_times_cp[i+j+1] - buffer):
            event_end = max(exit_times_cp[i+j+1], event_end)  # Account for case where first track starts first but also ends later
            j += 1
            if (i+j+1) >= len(entry_times_cp):
                break  # NEED to leave while loop if we are about to exceed the size our entry_times array -- !can't look forward beyond the last track!
        # entry time stays the same, set new exit time
        exit_times_cp[i] = event_end
        # Get rid of all the tracks that were blobbed together with the first one (i). The 'j' increment tells us how many tracks we combined
        exit_times_cp = np.delete(exit_times_cp, slice(i+1, i+j+1))
        entry_times_cp = np.delete(entry_times_cp, slice(i+1, i+j+1))
        i += 1

    # Filter out short events and rename these arrays; they now represent the start and end times of each combined event
    event_start_times = entry_times_cp[exit_times_cp - entry_times_cp >= np.timedelta64(min_dur, 's')]
    event_end_times = exit_times_cp[exit_times_cp - entry_times_cp >= np.timedelta64(min_dur, 's')]
    
    # Create a list of tuples of the event intervals (audible intervals)
    audible_intervals = [(a, b) for a, b in zip(event_start_times, event_end_times)]
    audible_times = event_end_times - event_start_times        # List of event durations in timedelta64 format
    total_audible = sum(audible_times)/np.timedelta64(1,'s')   # Add up all event durations to get total audible time (convert to float in seconds)
    
    # Account for event-less time at beginning and end of timeframe in question
    # Note that if timeframe starts/ends with an event, no inaudible time will actually be added
    inaudible_begins = np.insert(event_end_times, 0, start_date)
    inaudible_ends = np.append(event_start_times, end_date)
    # Create a list of tuples of the inaudible intervals (time in between audible intervals)
    inaudible_intervals = [(a, b) for a, b in zip(inaudible_begins, inaudible_ends)]
    inaudible_times = inaudible_ends - inaudible_begins            # List of noise-free interval durations in timedelta64 format
    total_inaudible = sum(inaudible_times)/np.timedelta64(1,'s')   # Add up all noise-free durations to get inaudible time (convert to float in seconds)
    
    total_time = (end_date - start_date)/np.timedelta64(1,'s')     # Calculate total time of timeframe in seconds 
    
    # GENERAL STATS
    duration_list = audible_times/np.timedelta64(1,'s')  # list of event duration lengths in seconds
    NFI_list = inaudible_times/np.timedelta64(1,'s')     # list of noise-free interval lengths in seconds
    TA = 100 * total_audible / total_time                # (%) of total time with an audible event (Time Audible)
    event_df = pd.DataFrame(data={'start_time': event_start_times, 'end_time': event_end_times, 'duration': duration_list})
    NFI_df = pd.DataFrame(data={'start_time': inaudible_begins, 'end_time': inaudible_ends, 'duration': NFI_list})
    return event_df, NFI_df
    
def compute_audibility_stats(event_dataframe, start_date, end_date, months=list(range(1,13)), quantiles=.5):
    '''
    Computes audibility using the event dataframe from tracks2events(). 
    The output includes the following, both hourly and daily in %:
                                                                 - min
                                                                 - max
                                                                 - mean
                                                                 - standard deviation
                                                                 - median absolute deviation
                                                                 - desired quantiles                                                  
    
    Parameters
    ----------
    event_dataframe : GeoDataFrame
        A GeoDataFrame containing each noise event. Looks like:
            start_time | end_time | duration
    start_date : string
        The start date to begin calculating audibility, formatted as 'yyyy-mm-dd'. Partitions by event start date
    end_date : string
        The end date to stop calculating audibility, formatted as 'yyyy-mm-dd'. Partitions by event start date
    months : int or list of ints (between 1 and 12)
        Default is the full year, an optional input to specify the months of interest as a list of integers, 1-12. 
        This is helpful for highly seasonal flight patterns, such as Denali's summer vs winter splits.
    quantiles : float or list of floats (between 0 and 1)
        Default is .5 (the median), specifies which quantiles to output. E.g., [.1, .5., .9] will output 10th, 50th, and 90th quantiles
    
        
    Returns
    -------
    daily_audibility_stats : GeoDataFrame
        A GeoDataFrame containing the daily audibility % min, max, mean, standard deviation, median absolute deviation, and quantiles
    hourly_audibility_stats : GeoDataFrame
        A GeoDataFrame containing the hourly audibility % min, max, mean, standard deviation, median absolute deviation, and quantiles
    '''        
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
    print("Computing daily and hourly audibility metrics (% audible time)")
            
    # Extract events within the start and end date
    event_df = event_dataframe[(event_dataframe.start_time >= start_date) & (event_dataframe.start_time <= end_date)]
    # Filter just the months specified in the 'months' parameter
    event_df = event_df[[date.month in months for date in event_df.start_time]]

    # Create a datetimeindex that spans all days within the start and end date; add column time_audible
    daily_audibility_df = pd.DataFrame(index=pd.date_range(start_date, end_date, freq='d'), columns=['time_audible'])
    daily_audibility_df = daily_audibility_df[[date.month in months for date in daily_audibility_df.index]]
    for id, event in event_df.groupby(event_df.start_time.dt.floor('d')):
        daily_audibility_df['time_audible'].loc[event.start_time.dt.floor('d')] = 100 * event.duration.sum() / (60*60*24)
    daily_audibility_df.fillna(0, inplace=True)  # Fill in any days that had no audibile events with 0
    # Generate a stats summary dataframe for daily audibility
    daily_audibility_stats = daily_audibility_df.agg(['min', 'max','mean', 'std', stats.median_abs_deviation, 'sem'])
    # Add in quantiles as specified by user
    quantile_string = [str(int(100*q))+'% quantile' for q in quantiles]  # Generates index label for each quantile
    for i, q in enumerate(quantiles):
        daily_audibility_stats.loc[quantile_string[i]] = daily_audibility_df.quantile(q)  # Calculates each quantile
        
    # Create a datetimeindex that spans all hours within the start and end date; add column time_audible
    hourly_audibility_df = pd.DataFrame(index=pd.date_range(start_date, end_date, freq='h'), columns=['time_audible'])
    hourly_audibility_df = hourly_audibility_df[[date.month in months for date in hourly_audibility_df.index]]
    for id, event in event_df.groupby(event_df.start_time.dt.floor('h')):
        hourly_audibility_df['time_audible'].loc[event.start_time.dt.floor('h')] = 100 * event.duration.sum() / (60*60)
    hourly_audibility_df.fillna(0, inplace=True) # Fill in any hours that had no audibile events with 0
    # Generate a stats summary dataframe for hourly audibility
    hourly_audibility_stats = hourly_audibility_df.agg(['min', 'max', 'mean', 'std', stats.median_abs_deviation, 'sem'])
    # Add in quantiles as specified by user
    quantile_string = [str(int(100*q))+'% quantile' for q in quantiles]  # Generates index label for each quantile
    for i, q in enumerate(quantiles):
        hourly_audibility_stats.loc[quantile_string[i]] = hourly_audibility_df.quantile(q)  # Calculates each quantile
    
    return daily_audibility_stats, hourly_audibility_stats

def compute_NFI_stats(NFI_dataframe, start_date, end_date, months=list(range(1,13)), quantiles=.5):
    '''
    Computes noise-free interval (NFI) stats using the event dataframe from tracks2events(). 
    The output includes the following NFI stats in timestamp format:
                                                                         - min
                                                                         - max
                                                                         - mean
                                                                         - standard deviation
                                                                         - median absolute deviation
                                                                         - desired quantiles                                                  
    
    Parameters
    ----------
    NFI_dataframe : GeoDataFrame
        A GeoDataFrame containing each noise-free interval. Looks like:
            start_time | end_time | duration
    start_date : string
        The start date to begin calculating NFI stats, formatted as 'yyyy-mm-dd'. Partitions by event start date
    end_date : string
        The end date to stop calculating NFI stats, formatted as 'yyyy-mm-dd'. Partitions by event start date
    months : int or list of ints (between 1 and 12)
        Default is the full year, an optional input to specify the months of interest as a list of integers, 1-12. 
        This is helpful for highly seasonal flight patterns, such as Denali's summer vs winter splits.
    quantiles : float or list of floats (between 0 and 1)
        Default is .5 (the median), specifies which quantiles to output. E.g., [.1, .5., .9] will output 10th, 50th, and 90th quantiles
    
        
    Returns
    -------
    NFI_stats : GeoDataFrame
        A GeoDataFrame containing the NFI duration min, max, mean, standard deviation, median absolute deviation, and quantiles
    '''    
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
    print("Computing noise-free interval stats (in timedelta format)")
    
    # Extract events within the start and end date
    NFI_df = NFI_dataframe[(NFI_dataframe.start_time >= start_date) & (NFI_dataframe.start_time <= end_date)]
    # Filter just the months specified in the 'months' parameter
    NFI_df = NFI_df[[date.month in months for date in NFI_df.start_time]]
    # Generate a stats summary dataframe for noise-free interval duration
    NFI_stats = NFI_df.duration.agg(['min', 'max', 'mean', 'std', stats.median_abs_deviation, 'sem'])
    
    # Add in quantiles as specified by user
    quantile_string = [str(int(100*q))+'% quantile' for q in quantiles]  # Generates index label for each quantile
    for i, q in enumerate(quantiles):
       NFI_stats.loc[quantile_string[i]] = NFI_df.duration.quantile(q)   # Calculates each quantile
    NFI_stats = NFI_stats * np.timedelta64(1,'s')  # Convert to timedelta64 format (more human readable)
    return NFI_stats

def compute_event_stats(event_dataframe, start_date, end_date, months=list(range(1,13)), quantiles=.5):
    '''
    Computes event count stats using the event dataframe from tracks2events(). 
    The output includes the following, both hourly and daily in # of events:
                                                                             - min
                                                                             - max
                                                                             - mean
                                                                             - standard deviation
                                                                             - median absolute deviation
                                                                             - desired quantiles                                                  
    
    Parameters
    ----------
    event_dataframe : GeoDataFrame
        A GeoDataFrame containing each noise event. Looks like:
            start_time | end_time | duration
    start_date : string
        The start date to begin calculating event stats, formatted as 'yyyy-mm-dd'. Partitions by event start date
    end_date : string
        The end date to stop calculating event stats, formatted as 'yyyy-mm-dd'. Partitions by event start date
    months : int or list of ints (between 1 and 12)
        Default is the full year, an optional input to specify the months of interest as a list of integers, 1-12. 
        This is helpful for highly seasonal flight patterns, such as Denali's summer vs winter splits.
    quantiles : float or list of floats (between 0 and 1)
        Default is .5 (the median), specifies which quantiles to output. E.g., [.1, .5., .9] will output 10th, 50th, and 90th quantiles
    
        
    Returns
    -------
    daily_event_stats : GeoDataFrame
        A GeoDataFrame containing the daily event count min, max, mean, standard deviation, median absolute deviation, and quantiles
    hourly_event_stats : GeoDataFrame
        A GeoDataFrame containing the hourly event count min, max, mean, standard deviation, median absolute deviation, and quantiles
    '''
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
    print("Computing daily and hourly event metrics (# of events)")
            
    # Extract events within the start and end date
    event_df = event_dataframe[(event_dataframe.start_time >= start_date) & (event_dataframe.start_time <= end_date)]
    # Filter just the months specified in the 'months' parameter
    event_df = event_df[[date.month in months for date in event_df.start_time]]

    # Create a datetimeindex that spans all days within the start and end date; add column event_count
    daily_event_df = pd.DataFrame(index=pd.date_range(start_date, end_date, freq='d'), columns=['event_count'])
    daily_event_df = daily_event_df[[date.month in months for date in daily_event_df.index]]  # Filter by 'months'
    for id, event in event_df.groupby(event_df.start_time.dt.floor('d')):
        daily_event_df['event_count'].loc[event.start_time.dt.floor('d')] = len(event)  # event count is just length of dataframe group
    daily_event_df.fillna(0, inplace=True)    # Fill in any days that had no audibile events with 0
    # Generate a stats summary dataframe for daily event count
    daily_event_stats = daily_event_df.agg(['min', 'max','mean', 'std', stats.median_abs_deviation, 'sem'])
    
    # Add in quantiles as specified by user
    quantile_string = [str(int(100*q))+'% quantile' for q in quantiles]  # Generates index label for each quantile
    for i, q in enumerate(quantiles):
        daily_event_stats.loc[quantile_string[i]] = daily_event_df.quantile(q)  # Calculates each quantile

    # Create a datetimeindex that spans all hours within the start and end date; add column event_count
    hourly_event_df = pd.DataFrame(index=pd.date_range(start_date, end_date, freq='h'), columns=['event_count'])
    hourly_event_df = hourly_event_df[[date.month in months for date in hourly_event_df.index]]  # Filter by 'months'
    for id, event in event_df.groupby(event_df.start_time.dt.floor('h')):
        hourly_event_df['event_count'].loc[event.start_time.dt.floor('h')] = len(event)  # event count is just length of dataframe group
    hourly_event_df.fillna(0, inplace=True)   # Fill in any hours that had no audibile events with 0
    # Generate a stats summary dataframe for hourly event count
    hourly_event_stats = hourly_event_df.agg(['min', 'max', 'mean', 'std', stats.median_abs_deviation, 'sem'])
    # Add in quantiles as specified by user    
    quantile_string = [str(int(100*q))+'% quantile' for q in quantiles]   # Generates index label for each quantile
    for i, q in enumerate(quantiles):
        hourly_event_stats.loc[quantile_string[i]] = hourly_event_df.quantile(q)  # Calculates each quantile
    
    return daily_event_stats, hourly_event_stats

def compute_duration_stats(event_dataframe, start_date, end_date, months=list(range(1,13)), quantiles=.5):
    '''
    Computes event duration stats using the event dataframe from tracks2events(). 
    The output includes the following duration stats in timestamp format:
                                                                         - min
                                                                         - max
                                                                         - mean
                                                                         - standard deviation
                                                                         - median absolute deviation
                                                                         - desired quantiles                                                  
    
    Parameters
    ----------
    event_dataframe : GeoDataFrame
        A GeoDataFrame containing each noise event. Looks like:
            start_time | end_time | duration
    start_date : string
        The start date to begin calculating duration stats, formatted as 'yyyy-mm-dd'. Partitions by event start date
    end_date : string
        The end date to stop calculating duration stats, formatted as 'yyyy-mm-dd'. Partitions by event start date
    months : int or list of ints (between 1 and 12)
        Default is the full year, an optional input to specify the months of interest as a list of integers, 1-12. 
        This is helpful for highly seasonal flight patterns, such as Denali's summer vs winter splits.
    quantiles : float or list of floats (between 0 and 1)
        Default is .5 (the median), specifies which quantiles to output. E.g., [.1, .5., .9] will output 10th, 50th, and 90th quantiles
    
        
    Returns
    -------
    duration_stats : GeoDataFrame
        A GeoDataFrame containing the event duration min, max, mean, standard deviation, median absolute deviation, and quantiles
    '''
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
    print("Computing event duration stats (in timedelta format)")
    
    # Extract events within the start and end date
    event_df = event_dataframe[(event_dataframe.start_time >= start_date) & (event_dataframe.start_time <= end_date)]
    # Filter just the months specified in the 'months' parameter
    event_df = event_df[[date.month in months for date in event_df.start_time]]
    # Generate a stats summary dataframe for event duration
    duration_stats = event_df.duration.agg(['min', 'max', 'mean', 'std', stats.median_abs_deviation, 'sem'])
    
    # Add in quantiles as specified by user
    quantile_string = [str(int(100*q))+'% quantile' for q in quantiles]  # Generates index label for each quantile
    for i, q in enumerate(quantiles):
       duration_stats.loc[quantile_string[i]] = event_df.duration.quantile(q)   # Calculates each quantile
    duration_stats = duration_stats * np.timedelta64(1,'s')  # Convert to timedelta64 format (more human readable)
    return duration_stats

def get_all_stats(event_df, NFI_df, start_date, end_date, months=list(range(1,13)), quantiles=.5):

    start_date = np.datetime64(start_date)  # Convert to datetime64
    end_date = np.datetime64(end_date)      # Convert to datetime64

    duration_stats = compute_duration_stats(event_df, start_date, end_date, months=months, quantiles=quantiles)
    daily_event_stats, hourly_event_stats = compute_event_stats(event_df, start_date, end_date, months=months, quantiles=quantiles)
    NFI_stats = compute_NFI_stats(NFI_df, start_date, end_date, months=months, quantiles=quantiles)
    daily_audibility_stats, hourly_audibility_stats = compute_audibility_stats(event_df, start_date, end_date, months=months, quantiles=quantiles)
    return duration_stats, daily_event_stats, hourly_event_stats, daily_audibility_stats, hourly_audibility_stats, NFI_stats