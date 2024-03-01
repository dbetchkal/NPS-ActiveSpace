#Helper Functions
import imports

def contiguous_regions(condition):

    """
    Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    """

    # Find the indicies of changes in "condition"
    d = imports.np.diff(condition)
    idx, = d.nonzero() 

    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = imports.np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = imports.np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)

    return idx

def calculate_duration_summary(noise_intervals):

    '''
    This function computes the duration of noise event in `noise_intervals`.
    It's called as part of `add_ambience()` and works behind the scenes.
    inputs
    ------
    self
    outputs
    -------
    a 2D numpy array of sound pressure level metrics as 'observed' in `self.full_record`:
        [0] a list of each event's duration
        [1] the mean duration
        [2] the standard deviation of the durations
        [3] the median duration
        [4] the median absolute deviation of the durations
    '''

    # the durations, themselves
    duration_list = noise_intervals.T[1] - noise_intervals.T[0]

    # mean duration
    mean = imports.np.mean(duration_list)

    # standard deviation duration
    stdev = imports.np.std(duration_list)

    # median duration
    median = imports.np.percentile(duration_list, 50)

    # median absolute deviation of duration
    mad = imports.np.percentile(imports.np.absolute(duration_list - median), 50)

    # combine the results and update the class attribute
    out = imports.np.array([duration_list, mean, stdev, median, mad])

    # it's convenient to return the results
    return out


def adjust_noise_free_intervals(noise_free_intervals, noise_intervals):

    '''
    In this simulation our convention will be to have closed noise intervals.
    To achieve this, we need to bound our noise free intervals. 
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

        # ...the first noise free interval stays the same, and equals zero
        # the rest are + 1
        nfi_starts = nfi_starts + 1
        nfi_starts[0] = 0


    # ------- Account for different ending conditions -----------------------------------------

        # the record ends with noise...
    if(noise_intervals[-1, 0] == 0):

        # ...the last noise free interval (and thus ALL intervals) need to end one second earlier
        nfi_ends = nfi_ends - 1


    # the record ends with quietude...
    else:

        # ...the last noise free interval stays the same, and equals zero
        # the rest are - 1
        save = nfi_ends[-1]
        nfi_ends = nfi_ends - 1
        nfi_ends[-1] = save


    # reset attribute to these updated, correct values
    new_noise_free_interval = imports.np.array([nfi_starts, nfi_ends]).T
    
    return new_noise_free_interval

def load_activespace(u, s, y, gain, third_octave=True, crs=None, PROJ_DIR=r"V:\NMSim\01 SITES"):

    if gain < 0:
        sign = "-"
    else:
        sign = "+"
        
    if gain < 10:
        gain_string = "0" + str(imports.np.abs(int(10*gain)))
    else:
        gain_string = str(imports.np.abs(int(10*gain)))

    path = PROJ_DIR + imports.os.sep + u+s + imports.os.sep + u+s+str(y) + '_O_' + sign + gain_string + ".geojson"
    active_space = imports.gpd.read_file(path)
    
    if crs is not None:
        active_space = active_space.to_crs(crs)
    
    return active_space

def load_studyarea(u, s, y, crs=None, PROJ_DIR=r"V:\NMSim\01 SITES"):

    # load in the study area as well
    study_area_path = imports.glob.glob(PROJ_DIR + os.sep + u+s + imports.os.sep + u+s + '*study*area*.shp')[0]
    study_area = imports.gpd.read_file(study_area_path)
    
    if crs is not None:
        study_area = study_area.to_crs(crs)
        
    return study_area

def cosdir_azim(azim):
    
    az = imports.np.deg2rad(azim)
    cosa = imports.np.sin(az)
    cosb = imports.np.cos(az)
    
    return cosa,cosb

def circular_offset(WindDir, offset=None):
    '''Applies wind direction offset to North'''
    return (WindDir + offset)%360

def compute_HeadingVectors(start_point, heading, speed):
    
    cosa, cosb = cosdir_azim(heading)
    norm = imports.mcolor.Normalize(vmin=0, vmax=360)
    HSV_color = imports.cm.hsv(norm(heading))
    
    m = imports.interpolate.interp1d([0.1, 11.0],[21, 2300])
    
    if(type(start_point) == 'shapely.geometry.point.Point'):
        start_point = [p[0] for p in start_point.coords.xy]
        
    end_point = imports.Point(start_point[0]+(m(speed)*cosa), start_point[1]+(m(speed)*cosb))
    
    vector = imports.LineString((start_point, end_point))
    
    return vector, HSV_color

def circular_median(WindDir):
    
    WindDir_R = imports.np.deg2rad(WindDir)
    
    # calculate the median of the cosine and sine components
    median_COS = imports.np.nanmedian(imports.np.cos(WindDir_R))
    median_SIN = imports.np.nanmedian(imports.np.sin(WindDir_R))
    
    # Take the arctangent
    out = imports.np.arctan2(median_SIN, median_COS)
    
    return 180 + imports.np.rad2deg(out)

import math

def round_values(array):
    rounded_array = []
    first_rounded = math.floor(array[0])
    rounded_array.append(first_rounded)
    for idx, value in enumerate(array[1:]):
        rounded_value = math.floor(value) + (round(value) - first_rounded)
        rounded_array.append(rounded_value)
    return rounded_array