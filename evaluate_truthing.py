# -*- coding: utf-8 -*-
"""
Loads flight truthing annotations and generates/saves active spaces for various omnidirectional sources to best fit
the audible points.

Created 2021-08-13 13:29

@author: Kirby Heck
"""

import os
import sys
import argparse
import time

import rasterio.plot

cd = os.getcwd()
sys.path.append(os.path.join(cd, "iyore-master"))
sys.path.append(os.path.join(cd, "soundDB-master"))
from soundDB import *
import iyore

# DENA RDS computer
RDS = r"\\inpdenards\overflights"

sys.path.append(os.path.join(RDS, "scripts"))
from query_tracks import query_tracks

# need the active space scripts as well
cd = os.getcwd()
sys.path.append(os.path.join(cd, "NMSIN-Python"))
from ActiveSpace import CreateActiveSpace
from active_space_utils import *

# ==================== global variables, parse command line arguments ==================

if len(sys.argv) == 1:
    u = input("Unit name (e.g. 'DENA'): ")
    s = input("Site name (e.g. 'KAHP'): ")
    y = int(input("Year: "))

    # maximum number of tracks to look through
    max_tracks = None

else:
    parser = argparse.ArgumentParser(description="Truthing active space calculations at sound monitoring stations in "
                                                 "Denali National Park")
    parser.add_argument('UNIT', help="Four letter park service unit code")
    parser.add_argument('SITE', help="Site character code")
    parser.add_argument('YEAR', help="Year to look for field data")
    parser.add_argument('--maxTracks', default=None,
                        help="Maximum number of tracks to load")
    parser.add_argument('--showplots', action='store_true', default=False,
                        help="Show plots at the end of truthing (saves automatically regardless)")

    args = parser.parse_args()

    u, s, y = args.UNIT, args.SITE, int(args.YEAR)
    max_tracks = args.maxTracks
    if max_tracks is not None:
        max_tracks = int(max_tracks)

print("Beginning evaluation for " + u + s)

PROJ_DIR = r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_ProjectDirectories" + os.sep + u + s
annotation_filename = u + s + str(y) + '_saved_annotations.csv'

# we'll also load the site metadata sheet
site_info = get_site_info(unit_name=u, site_name=s, year=y)

study_area = gpd.read_file(glob.glob(PROJ_DIR + os.sep + u + s + "*study*.shp")[0])
CRS = site_info.utm_zone.values[0]

site_utm = gpd.GeoDataFrame(geometry=[Point(site_info["x"], site_info["y"])], crs=site_info["utm_zone"].values[0])

# ================== find available years of annotations ==================

annotation_paths = glob.glob(PROJ_DIR + os.sep + u + s + "*_saved_annotations.csv")

if len(annotation_paths) == 0:
    print("No annotations found for " + u + s + ".")
    exit(-1)

years = [int(os.path.basename(filename)[8:12]) for filename in annotation_paths]

annotations = pd.DataFrame()

# merge all annotations from all years into one DataFrame
for i, yr in enumerate(years):
    annotations = pd.concat([annotations, pd.read_csv(annotation_paths[i], index_col=0)])

# ==================== query tracks =======================

tracks = None

for yr in years:
    if yr == 2021:
        archive_path = r'V:\Data 2021'
        archive = iyore.Dataset(archive_path)
        print("Switched archive to", archive_path)
    else:
        archive = iyore.Dataset(r"E:")

    print("Loading tracks for", yr)
    # return `datetime` objects representing days for which at least some NVSPL data exist
    sampled_dates = sorted(list(set([dt.datetime(int(e.year), int(e.month), int(e.day)) for
                                     e in archive.nvspl(unit=u, site=s, year=yr)])))

    if tracks is None:
        # load tracks from the database over NVSPL dates, using the study area as a mask
        tracks = query_tracks(connection_txt=r"\\inpdenards\overflights\config\connection_info.txt",
                              start_date=sampled_dates[0], end_date=sampled_dates[-1],
                              mask=study_area, clip_output=False,
                              aircraft_info=False)
    else:
        tracks = tracks.append(query_tracks(connection_txt=r"\\inpdenards\overflights\config\connection_info.txt",
                                            start_date=sampled_dates[0], end_date=sampled_dates[-1],
                                            mask=study_area, clip_output=False,
                                            aircraft_info=False))

ids = tracks.id.unique()

tracks = tracks.to_crs(CRS)  # make sure tracks are all in the same crs
valid_points = gpd.GeoDataFrame()

if max_tracks is not None:
    total_tracks = min(len(ids), max_tracks)
else:
    total_tracks = len(ids)

# at this point, iterate through the flight id's and make note of audible points by
# reading the annotations (no annotations = no audio record)
for k, _id in enumerate(ids):

    print('[{}/{}] Flight ID {}: '.format(k+1, total_tracks, _id))
    flight = tracks.loc[tracks["id"] == _id, :]

    # we need to sort the track chronologically for the splines to fit properly
    flight = flight.sort_values('ak_datetime')

    # ============= check for prior saved annotations ================

    prior_data = annotations.loc[annotations._id == _id]  # previously made one-row annotation

    if len(prior_data) > 0:  # keep track? Tracks with no data have no audio record and are thrown out
        # convert flights to UTM coordinates
        if flight.crs != CRS:
            flight_utm = flight.to_crs(CRS)  # this is slow
        else:
            flight_utm = flight

        try:
            flight_spline, closest_time = compute_spline(flight_utm, site_info)

        except TypeError as e:  # this may occur if the flight through the study area only has 1 point in it
            print(e)
            plt.close()
            continue  # skip the rest of this track

        except ValueError as e:  # STILL not sure why this happens... ?
            print(e)
            print(flight_spline)

        # ======================== now load the annotations ===========================

        if prior_data.valid.values:  # valid track? T/F
            if prior_data.audible.values:  # audible? T/F
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

    else:
        print('\tNo accompanying annotation or audio record for flight.')

    print('')  # newline

    if max_tracks is not None and k+1 >= max_tracks:
        break

#  ========= compute metrics based on the valid points identified ========

# check if there already exists an active spaces computation file
z_data = valid_points.geometry.z.values  # in [m]
alt_ft = int(np.mean(z_data) * 3.28084)  # [m] to [ft]

results_path = glob.glob(PROJ_DIR + os.sep + u + s + '_active_spaces_3octband_*ft.txt')

if len(results_path) > 0:
    print('Results found for ' + u + s +'. Skipping to plotting... ')

    res = pd.read_csv(results_path[0])
    res['geometry'] = gpd.GeoSeries.from_wkt(res['geometry'])
    results = gpd.GeoDataFrame(res, geometry='geometry')
    results.crs = CRS

else:  # create several different active spaces and produce a precision-recall plot
    source_dir = r'T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_TuningSources'
    n_sources = 51
    source_df = get_omni_sources(source_dir, n_sources=n_sources, upper=30, lower=-20)
    print('Average flight altitude for valid points: ', alt_ft, 'ft')

    results = gpd.GeoDataFrame(columns=['f1', 'precision', 'recall'])
    make_AS = CreateActiveSpace(unit_name=u, site_name=s, yr=y, project_dir=PROJ_DIR,
                                freq_bands=True)  #, ambience_source='Mennitt')

    print(f'Beginning {n_sources} active space calculations. \n')

    starttime = time.perf_counter()

    # loop through source file and create an active space for each one
    for i, row in enumerate(source_df.iterrows()):
        make_AS.set_source(row[1].full_path)  # set the source

        act_space = make_AS.run_model(altitude_ft=alt_ft, n_tracks=1, out_crs=CRS)
        if act_space is None:
            print(f"\n\t[{i + 1}/{n_sources}] Skipped active space for source", row[1].filename)
            continue  # noise source was too quiet (empty active space) or the script failed to produce a polygon

        f1, precision, recall, n_tot = compute_f1(valid_points, act_space)
        act_space['f1'] = f1
        act_space['precision'] = precision
        act_space['recall'] = recall
        act_space.index = [row[1].gain]  # index

        results = results.append(act_space)

        print(f"\n\t[{i + 1}/{n_sources}] Finished active space for source", row[1].filename)
        print(
            "\tF1 score: {:1.2}, Precision: {:1.2}, Recall: {:1.2} for {:0} points.\n".format(f1, precision,
                                                                                              recall, n_tot))

    # this needs to come after there is at least one row
    results.index.name = 'gain'
    results.to_wkt().to_csv(PROJ_DIR + os.sep + u + s + '_active_spaces_3octband_{:.0f}ft.txt'.format(alt_ft))

    endtime = time.perf_counter()

    print("\nFINISHED COMPUTING ACTIVE SPACE, average time {:.2f} s per polygon".format((endtime-starttime)/n_sources))


# ==================== plot the annotated flights =======================

# plotting CRS, lat/long looks nicer here...
plot_CRS = 'epsg:4326'

fig, ax = plt.subplots(1, 1, figsize=(6, 8))

# show basemap in lat/long

with rasterio.open(r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\DENA Park Brochure Map wi.tif") as raster:
    rasterio.plot.show(raster, ax=ax, alpha=0.6)

# plot active space estimate
max_f1 = results.loc[results.f1==results.f1.max()]
active_space = max_f1.geometry
active_space = active_space.to_crs(plot_CRS)
active_space.geometry.boundary.plot(ax=ax, color="blue", facecolor="blue", alpha=0.4, zorder=10,
                                    label="active space estimate")

# plot study area bounding box
study_area = study_area.to_crs(plot_CRS)
study_area.geometry.boundary.plot(ax=ax, ls="--", color="navy")

# # plot microphone position
ax.plot(site_info["long"], site_info["lat"], ls="", marker="*", ms=6, color="black",
        zorder=3, label="microphone\n" + u + s)

# audible vs inaudible points:
plot_points = valid_points.to_crs(plot_CRS)  # make a copy of these for plotting in WGS84
plot_points.loc[plot_points.audible].plot(ax=ax, color='deepskyblue', alpha=0.05, markersize=3, zorder=3,
                                          label="Audible points")
plot_points.loc[~plot_points.audible].plot(ax=ax, color='red', alpha=0.05, markersize=3, zorder=2,
                                           label="Inaudible points")

# create a legend at the bottom center
leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.4), markerscale=2)
for lh in leg.legendHandles:
    lh.set_alpha(1)  # set legend handles to opaque for visibility

xmin, ymin, xmax, ymax = study_area.total_bounds
pad = np.array([(xmax - xmin) * 0.1, (ymax - ymin) * 0.1])

# this will result in a square map
ax.set(xlim=(xmin - pad[0], xmax + pad[0]), ylim=(ymin - pad[1], ymax + pad[1]))
ax.set_title(u+s+" Annotated audible overflights")
ax.tick_params(axis='both', labelsize=6)

fig.savefig(PROJ_DIR + os.sep + u + s + '_active_truthing_3octband_{:}ft.png'.format(alt_ft),
            dpi=300, bbox_inches="tight", facecolor="white")

if args.showplots:
    plt.show()
else:
    plt.close()
