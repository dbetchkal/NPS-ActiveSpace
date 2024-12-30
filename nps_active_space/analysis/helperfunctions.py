#Helper Functions
import imports


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
