
"""
Script for annotating overflight noise impacts at sound monitoring stations in Denali National Park.

usage: flight_track_truthing.py [-h] UNIT SITE YEAR. If no positional arguments are given, they will be asked for via
user input through the command line.

positional arguments:
  UNIT        Four letter park service unit code
  SITE        Site character code
  YEAR        Year to look for overflights and an acoustic record

optional arguments:
  -h, --help  show this help message


Created 2021-07-15 14:22 by Kirby Heck

"""
import argparse
import os
import sys

import matplotlib.dates as mdates
import rasterio.plot

cd = os.getcwd()
# sys.path.append(os.path.join(cd, "iyore-master"))
# sys.path.append(os.path.join(cd, "soundDB-master"))
from soundDB import *
import iyore

# DENA RDS computer
RDS = r"\\INPDENATERM01\overflights"

sys.path.append(os.path.join(RDS, "scripts"))
from query_tracks import query_tracks

from active_space_utils import *

# ==================== global variables, parse command line arguments ==================

if len(sys.argv) == 1:
    u = input("Unit name (e.g. DENA): ")
    s = input("Site name (e.g. KAHP): ")
    y = int(input("Year: "))

else:
    parser = argparse.ArgumentParser(description="Annotating overflight noise impacts at sound monitoring stations in "
                                                 "Denali National Park")
    parser.add_argument('UNIT', help="Four letter park service unit code")
    parser.add_argument('SITE', help="Site character code")
    parser.add_argument('YEAR', help="Year to look for overflights and an acoustic record")
    args = parser.parse_args()

    u, s, y = args.UNIT, args.SITE, int(args.YEAR)

print("Beginning truthing for " + u + s + str(y) + "...\n")

# ============== Load global variables and flights/acoustic record ==================

PROJ_DIR = r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_ProjectDirectories" + os.sep + u + s
annotation_filename = u+s+str(y)+'_saved_annotations.csv'

# we'll also load the site metadata sheet
site_info = get_site_info(unit_name=u, site_name=s, year=y)

# the year may need to be updated if there was no acoustic record from that year
y = int(site_info.year.values[0])

# NVSPL archive
# TODO: this will inevitably change once the data are archived
if y == 2021:
    archive = iyore.Dataset(r'V:\Data 2021')
else:
    archive = iyore.Dataset(r"E:")

# find the study area and active space shape files
study_area = gpd.read_file(glob.glob(PROJ_DIR + os.sep + u + s + "*study*.shp")[0])

# worth while checking coordinates are in UTM
CRS = site_info.utm_zone.values[0]
if study_area.crs != CRS:
    study_area = study_area.to_crs(CRS)

# ==================== query tracks =======================

# return `datetime` objects representing days for which at least some NVSPL data exist
sampled_dates = sorted(list(set([dt.datetime(int(e.year), int(e.month), int(e.day)) for
                                 e in archive.nvspl(unit=u, site=s, year=y)])))

# load tracks from the database over NVSPL dates, using the study area as a mask
tracks = query_tracks(connection_txt=r"\\INPDENATERM01\overflights\config\connection_info.txt",
                      start_date=sampled_dates[0], end_date=sampled_dates[-1],
                      mask=study_area, clip_output=False,
                      aircraft_info=False)

# =========== load acoustic data from NVSPL ===============

# sort chronologically
# this is necessary for grouping in the next cell and then locating temporal boundaries
tracks = tracks.sort_values("ak_datetime")

# glean the unique dates from the GPS points in a numpy array
unique_hours = tracks["ak_datetime"].apply(lambda t: (t.year, t.month, t.day, t.hour)).unique()

# use `soundDB` to load only the relevant NVSPL files
# iterate through unique_hours and reformat into a list of dictionaries
items = [{"year": str(y),
          "month": "{:02d}".format(m),
          "day": "{:02d}".format(d),
          "hour": "{:02d}".format(h)} for y, m, d, h in unique_hours]

# load the acoustic record
nv = nvspl(archive, unit=u, site=s, items=items, year=y).combine()

# if there are more than one hour of data loaded, drop the outer (filename) index
if (len(nv.index.names) > 1):
    nv.index = nv.index.droplevel()

temporal_boundaries = [(_id, track["ak_datetime"].iloc[0], track["ak_datetime"].iloc[-1]) for
                       _id, track in tracks.groupby("id")]

site = gpd.GeoDataFrame(geometry=[Point(site_info["long"], site_info["lat"])], crs='epsg:4326')
site_utm = site.to_crs(CRS)


####################################################################
# iterate through each track and prompt the user's input on if the #
# flight is audible in the spectrogram                             #
####################################################################


# convert all flight tracks to the correct coordinate system
tracks = tracks.to_crs(CRS)
valid_points = gpd.GeoDataFrame()
end_truthing = False

annotations = pd.DataFrame(columns={'_id', 'valid', 'audible', 'starttime', 'endtime'})

# ================== try and load previous annotations ==================

if os.path.exists(PROJ_DIR + os.sep + annotation_filename):  # check if this file even exists
    if yes_no_button('Found prior annotations, load them?'):  # if it is exists, ask to load it
        annotations = pd.read_csv(PROJ_DIR + os.sep + annotation_filename, index_col=0)

# iterate through all IDs and annotate or load prior annotations

print(len(temporal_boundaries))
for _id, start, end in temporal_boundaries:

    print(_id, start, end)
    flight = tracks.loc[tracks["id"] == _id, :]

    # because the study area is subject to change, all saved times should be relative to the departure time,
    # not the `start` time
    takeoff = flight['departure_datetime'].iloc[0]
    offset = pd.Timedelta(start-takeoff)

    try:  # this is a long `try` block -- sloppy :(
        duration = (end - start).total_seconds()
        # pad = dt.timedelta(seconds=180)  # TODO: how to get the sliders to line up with NVSPL record if padding > 0?
        spect = nv.loc[start:end, "12.5":"20000"]  # time slice of SPL record

        # convert the NVSPL's nice datetime axis to numbers
        x_lims = mdates.date2num(spect.index)

        # compute the aspect ratio of the spectrogram...
        # this will throw an IndexError if x_lims is empty (i.e. no spectrogram/audio record)
        aspect = (x_lims[-1] - x_lims[0]) / (8 * 33)

        # convert flights to UTM coordinates
        if flight.crs != CRS:
            flight_utm = flight.to_crs(CRS)  # this is slow
        else:
            flight_utm = flight

        flight_spline, closest_time = compute_spline(flight_utm, site_info)
        intervals = compute_intervals(flight_utm, site_info)

    except IndexError as e:
        print(e)  # prints error message
        plt.close()

        if len(spect.index) == 0:
            print("\tNo spectrogram to accompany flight path, dropped flight ID " + str(_id))
            # remove tracks without relevant spectrogram
            tracks = tracks[tracks["id"] != _id]

        else:
            print("\tUnknown error has occured, proceeding with other data...")

        continue  # skip the rest of this track

    except TypeError as e:
        # this may occur if the flight through the study area only has 1 point in it
        print(e)
        plt.close()

        continue  # skip the rest of this track

    # ============= check for prior saved annotations ================

    prior_data = annotations.loc[annotations._id==_id]  # previously made one-row annotation

    if len(prior_data) > 0:  # keep track?
        if prior_data.valid.values:  # valid track?
            if prior_data.audible.values:  # audible?
                flight_spline['audible'] = np.all([flight_spline.time_audible >= prior_data.starttime.values[0],
                                                   flight_spline.time_audible <= prior_data.endtime.values[0]], axis=0)

            else:  # inaudible flight
                flight_spline['audible'] = False

            if len(valid_points) == 0:
                valid_points = flight_spline  # for some reason this is necessary to transfer the CRS over?
            else:
                valid_points = valid_points.append(flight_spline, ignore_index=True)

            print('\tAnnotation found for flight ID ' + str(_id) + ', loaded data...')

        else:
            pass

        continue  # don't plot; already annotated. Proceed with other data

    # ===================== plotting with interactive plots =======================

    fig, ax = plt.subplots(1, 2, figsize=(15, 5),
                           gridspec_kw={'width_ratios': [1, 2]})
    plt.subplots_adjust(bottom=0.25)

    ax[1].set_title("Audio record, Flight ID: " + str(_id), loc="left")

    # origin "lower" flips y-axis; plot spectrogram data
    ax[1].imshow(spect.T, origin="lower", aspect=aspect, cmap="plasma",
                 extent=[x_lims[0], x_lims[-1], 0, spect.shape[1]],
                 interpolation=None, vmin=-10, vmax=80)

    # add in the time of closest approach in red
    ax[1].axvline(mdates.date2num(closest_time), alpha=0.7, color="red", zorder=2, linewidth=3)

    ax[1].set_yticks(np.arange(spect.shape[1])[::6])
    ax[1].set_yticklabels(spect.columns.astype('float')[::6])
    ax[1].set_ylabel("Freq. (Hz)", labelpad=15)

    # tell matplotlib that the numeric axis should be formatted as dates
    ax[1].xaxis_date()
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%b-%d\n%H:%M"))  # tidy them!

    # ================== finished with spectrogram, now for map annotation ===================

    if flight.crs != CRS:
        flight = flight.to_crs(CRS)
    flight.plot(ax=ax[0], color="blue", zorder=1, markersize=3)

    # we could plot with geopandas, but using pyplot will allow for a dynamic plot that reacts to the range slider input
    # that is, xy_UTM is a copy of the flight_spline geometry for plotting ONLY
    xy_UTM = flight_spline.geometry.bounds.iloc[:, 0:2]
    xy_UTM['time_audible'] = flight_spline['time_audible']  # when we can actually HEAR the noise, mic time

    # `highlight` is a dynamic line object that will update with the range slider
    highlight, = ax[0].plot(xy_UTM.minx, xy_UTM.miny, lw=5, color='deepskyblue', ls='-', zorder=1, alpha=0.4)

    # draw a box around the study area
    study_area.geometry.boundary.plot(ax=ax[0], ls="--", lw=0.5, color="blue")
    ax[0].plot(site_utm.geometry.x, site_utm.geometry.y, ls="", marker="x", ms=10, color="red", zorder=10)
    # rasterio.plot.show(raster, ax=ax[0], alpha=0.6)  # show basemap? - makes annotating much slower

    # glean the spatial extent of the points
    xmin, ymin, xmax, ymax = study_area.total_bounds

    # this will result in a square map
    ax[0].set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax[0].set_aspect((xmax - xmin) / (ymax - ymin))
    ax[0].set_title("Point positions within Study Area for " + flight["registration"].iloc[0], loc="left")
    ax[0].tick_params(axis='both', labelsize=6)
    ax[0].ticklabel_format(style='plain')  # disable scientific notation
    # plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)

    # =============== range slider & interactive buttons ===================

    # Create the RangeSlider object
    slider_ax = plt.axes([0.45, 0.1, 0.45, 0.03])
    slider = RangeSlider(slider_ax, "Audible Window [flight duration, s]", 0, duration)
    slider.set_val([-60, 60] + (closest_time - start).total_seconds())

    # Create the moving vertical lines on the histogram with axvline()
    lower_limit_line = ax[1].axvline(mdates.date2num(start + dt.timedelta(seconds=slider.val[0])),
                                     ls="--", alpha=0.7, color="white", zorder=2, linewidth=1)
    upper_limit_line = ax[1].axvline(mdates.date2num(start + dt.timedelta(seconds=slider.val[1])),
                                     ls="--", alpha=0.7, color="white", zorder=2, linewidth=1)


    def slider_update(val):
        # The val passed to a callback by the RangeSlider will
        # be a tuple of (min, max)

        # read the slider values
        lower_t = start + dt.timedelta(seconds=slider.val[0])
        upper_t = start + dt.timedelta(seconds=slider.val[1])
        lower_date = mdates.date2num(lower_t)
        upper_date = mdates.date2num(upper_t)

        # update the vertical lines on the spectrogram
        lower_limit_line.set_xdata([lower_date, lower_date])
        upper_limit_line.set_xdata([upper_date, upper_date])

        subset = xy_UTM.loc[np.all([flight_spline.time_audible >= lower_t,
                                    flight_spline.time_audible <= upper_t],
                                   axis=0)]  # subset of entire flight spline
        highlight.set_data(subset.minx, subset.miny)

        # Redraw the figure to ensure it updates
        fig.canvas.draw_idle()


    slider_update(None)  # update highlight
    slider.on_changed(slider_update)

    # add radio buttons (selection) too:
    ax_buttons = plt.axes([0.125, 0.02, 0.1, 0.18])
    radio_buttons = RadioButtons(ax_buttons, ['Audible', 'Inaudible', 'Unknown', 'Exit'], active=None)


    def radio_clicked(event):
        # for radio buttons, `event` is a string of the label selected
        # This is the audibility selector panel on the annotations screen

        global valid_points  # need the global `version` of this (?)
        global end_truthing  # there must be a better way than this...
        global annotations

        if (event == "Audible"):
            # these are valid points; mark 'true' or 'false' to their audibility
            lower_t = start + dt.timedelta(seconds=slider.val[0])
            upper_t = start + dt.timedelta(seconds=slider.val[1])

            flight_spline['audible'] = np.all([flight_spline.time_audible >= lower_t,
                                               flight_spline.time_audible <= upper_t], axis=0)
            if len(valid_points) == 0:
                valid_points = flight_spline
            else:
                valid_points = valid_points.append(flight_spline, ignore_index=True)

            annotations = annotations.append({'_id': _id,
                                              'valid': True,
                                              'audible': True,
                                              'starttime': lower_t,
                                              'endtime': upper_t}, ignore_index=True)

        elif (event == 'Inaudible'):
            # these are also valid points; they are all simply inaudible
            flight_spline['audible'] = False

            if len(valid_points) == 0:
                valid_points = flight_spline
            else:
                valid_points = valid_points.append(flight_spline, ignore_index=True)

            annotations = annotations.append({'_id': _id,
                                              'valid': True,
                                              'audible': False}, ignore_index=True)

        elif (event == "Unknown"):
            annotations = annotations.append({'_id': _id,
                                              'valid': False}, ignore_index=True)

        else:
            end_truthing = True

        plt.close('all')


    radio_buttons.on_clicked(radio_clicked)
    plt.show()

    if end_truthing:
        break
    # ============ end for loop ==============

# ask to save annotations
if yes_no_button('Save annotations?'):
    annotations.to_csv(PROJ_DIR + os.sep + annotation_filename)
    print('Saved annotations to ' + PROJ_DIR + os.sep + annotation_filename)

print("\nFINISHED.")

if yes_no_button("Plot results?"):
    # ==================== plot the annotated flights =======================

    # plotting CRS, lat/long looks nicer here...
    plot_CRS = 'epsg:4326'

    fig, ax = plt.subplots(1, 1, figsize=(6, 9))

    # show basemap in lat/long
    # with rasterio.open(r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\DENA Park Brochure Map wi.tif") as raster:
    #     rasterio.plot.show(raster, ax=ax, alpha=0.6)

    # plot study area bounding box
    study_area = study_area.to_crs(plot_CRS)
    study_area.geometry.boundary.plot(ax=ax, ls="--", color="navy")

    # # plot microphone position
    ax.plot(site_info["long"], site_info["lat"], ls="", marker="*", ms=6, color="black",
            zorder=3, label="microphone\n" + u + s + str(y))

    # audible vs inaudible points:
    valid_points = valid_points.to_crs(plot_CRS)
    valid_points.loc[valid_points.audible].plot(ax=ax, color='deepskyblue', alpha=0.05, markersize=3, zorder=3,
                                                label="Audible points")
    valid_points.loc[~valid_points.audible].plot(ax=ax, color='red', alpha=0.05, markersize=3, zorder=2,
                                                 label="Inaudible points")

    # create a legend at the bottom center
    leg = plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), markerscale=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)  # set legend handles to opaque for visibility

    xmin, ymin, xmax, ymax = study_area.total_bounds
    pad = np.array([(xmax - xmin) * 0.1, (ymax - ymin) * 0.1])

    # this will result in a square map
    ax.set(xlim=(xmin - pad[0], xmax + pad[0]), ylim=(ymin - pad[1], ymax + pad[1]))
    ax.set_title("Annotated audible overflights segments")
    ax.tick_params(axis='both', labelsize=6)

    plt.show()
