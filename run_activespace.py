# -*- coding: utf-8 -*-
"""
Run a park-wide raster computation of active spaces for Denali National Park and Preserve.

Created 2021-08-16 13:10

@author: Kirby Heck
"""

import os
import sys
from active_space_utils import *
from ModelManager import ModelManager
import argparse
import time


def main():

    # ========== parse command line arguments ============

    if len(sys.argv) == 1:  # no command line args given, prompt user for input
        out_dir = input("Output directory path: ")
    else:
        parser = argparse.ArgumentParser(
            description="Compute active spaces for Denali National Park")
        parser.add_argument('out_dir', help='Output directory for all files, metadata, etc. ')
        parser.add_argument('--cellsize', default=10000, help='Cell size (in [m]) for the rasterized data')
        parser.add_argument('--n_jobs', default=0.8, help='Number of multiprocessing processes. If n_jobs is between '
                                                      '0 and 1, will use that fraction of total computational cores.')
        parser.add_argument('--out_filename', default='active_spaces_{}m.txt',
                            help='Output filename for the written polygons, default: active_spaces_{cellsize}m.txt')
        parser.add_argument('--plot_results', action='store_true', default=False,
                            help='Toggles plotting results')
        parser.add_argument('--load_previous', action='store_true', default=True,
                            help='Attempts to load runs from a previous active space calculation.')
        parser.add_argument('--group_runs', default=100,
                            help='Groups parallelized runs into chunks to periodically save progress; default 100')
        parser.add_argument('--omni_gain', default=23,
                            help='Gain applied to omnidirectional noise source for a Cessna 206; default +23 [dB]')
        # other args go here

        args = parser.parse_args()
        out_dir = args.out_dir
        cellsize = int(args.cellsize)
        out_filename = args.out_filename.format(cellsize)
        group_runs = int(args.group_runs)
        n_jobs = float(args.n_jobs)
        gain = float(args.omni_gain)

    # ========== create park raster ===========

    # load DENA boundary:
    DENA_shp = r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NIMSIM_DEMSs\DENA_outline.shp"
    DENA_poly = gpd.read_file(DENA_shp).geometry

    raster_path = os.path.join(out_dir, 'active_space_blank_{:d}m.tif'.format(cellsize))

    if args.load_previous:
        # try to load previous results
        try:
            print("Attempting to load previous results from " + os.path.join(out_dir, out_filename))
            poly_list = get_run_results(os.path.join(out_dir, out_filename))
            coords_ls = get_coord_list(poly_list)
            print("{} results found and loaded... \n".format(len(poly_list)-len(coords_ls)))

        except FileNotFoundError as e:
            # no results found, either the path has changed or no results exist
            print("Could not find file: " + os.path.join(out_dir, out_filename))
            args.load_previous = False

    if not args.load_previous:  # kind of a strange way to enter this block...
        poly_list = None
        coord_key = create_blank_raster(DENA_poly, raster_path, cellsize=cellsize)
        coords_ls = get_coord_list(coord_key)

    # ========== Initiate ModelManager object(s) ===========

    sound_source_dir = r'T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_TuningSources'
    sound_sources = get_omni_sources(sound_source_dir, upper=gain, lower=gain)
    sound_src = sound_sources.full_path.values[0]
    init_args = {'project_dir': 'tmp',
                 'ambience_source': 'Mennitt',
                 'source_path': sound_src}

    ak_albers = 'epsg:6393'

    n_runs = int(np.ceil(len(coords_ls)/group_runs))

    print("=============== Creating {} active spaces ===============".format(len(coords_ls)))

    starttime = time.perf_counter()

    # for each 'chunk', run the ModelManager and save the (ongoing) polygon list
    for i in range(n_runs):

        # Grab the subset of coordinates:
        if i == n_runs-1:  # last run
            coords_sub = coords_ls[i*group_runs:]
        else:
            coords_sub = coords_ls[i*group_runs:(i+1)*group_runs]

        # create modelmanager for the given coordinates list
        mm = ModelManager(n_jobs=n_jobs, coords_ls=coords_sub, crs=ak_albers, init_args=init_args)
        # run the model!
        polys = mm.run_all()

        # because we are saving as a .wkt, let's be explicit about saving the coordinate system too
        polys = polys.assign(crs=ak_albers)

        print("\n\n End of run group {}, merging polygons to the existing record.".format(i+1))

        if poly_list is None:  # on first run, create the list of polygons
            poly_list = polys.merge(coord_key, how='outer', left_index=True, right_index=True)
        else:  # start with the existing list and merge these in
            poly_list.loc[polys.index, polys.columns] = polys

        # save progress
        print("Saving polygons to well-known text... ")
        poly_list.to_wkt().to_csv(os.path.join(out_dir, out_filename))
        print("Save complete, backup created at '{}'\n\n".format(os.path.join(out_dir, out_filename)))

    endtime = time.perf_counter()
    print('\n\n====================== FINISHED ======================')
    print('\tTotal computations: {:d}'.format(len(coords_ls)))
    print('\tAverage time per computation: {:.2f} seconds'.format((endtime-starttime)/(len(coords_ls))))

    # ================ Plot results ==================

    if args.plot_results:
        coords_wgs = gpd.GeoDataFrame(geometry=[Point(row[1]) for row in coords_ls], crs=ak_albers).to_crs('epsg:4326')
        poly_list_wgs = poly_list.to_crs('epsg:4326')
        ras_src = r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\DENA Park Brochure Map wi.tif"

        fig, ax = plt.subplots(figsize=(6, 7))
        with rasterio.open(ras_src) as ras:
            rasterio.plot.show(ras, ax=ax)
        poly_list_wgs.geometry.boundary.plot(ax=ax, alpha=0.4, color='blue', facecolor='blue', label='Active Space')
        coords_wgs.plot(ax=ax, marker='*', color='black', label='Site Location')

        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25))

        fig.savefig(os.path.join(out_dir, 'active_space_computations_{:d}m.png'.format(cellsize)),
                    dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()


if __name__ == '__main__':
    main()  # run main script
