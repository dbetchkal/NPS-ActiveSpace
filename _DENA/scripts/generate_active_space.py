import glob
import time
import os
from argparse import ArgumentParser

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

import iyore

import _DENA.resource.config as cfg
from _DENA import _DENA_DIR
from _DENA.resource.helpers import get_deployment, get_logger, get_omni_sources
from nps_active_space.utils import Annotations, coords_to_utm, Nvspl, NMSIM_bbox_utm
from nps_active_space.active_space import ActiveSpaceGenerator


if __name__ == '__main__':

    argparse = ArgumentParser()

    argparse.add_argument('-e', '--environment', required=True,
                          help="The configuration environment to run the script in.")
    argparse.add_argument('-u', '--unit', required=True,
                          help="Four letter unit code. E.g. DENA")
    argparse.add_argument('-s', '--site', required=True,
                          help="Four letter site code. E.g. TRLA")
    argparse.add_argument('-y', '--year', type=int, required=True,
                          help="Four digit year. E.g. 2018")
    argparse.add_argument('-a', '--ambience', default='nvspl', choices=['nvspl', 'mennitt'],
                          help='What type of ambience to use in NMSIM calculations.')
    argparse.add_argument('--headings', nargs='+', type=int, default=[0, 120, 240],
                          help='Headings of active spaces to dissolve.')
    argparse.add_argument('--omni-min', type=int, default=-20,
                          help='The minimum omni source to run the mesh for.')
    argparse.add_argument('--omni-max', type=int, default=30,
                          help='The maximum omni source to run the mesh for.')
    # argparse.add_argument('-n', '--max-tracks', type=int, default=None, # TODO
    #                       help='Maximum number of tracks to load')
    args = argparse.parse_args()

    cfg.initialize(f"{_DENA_DIR}/config", environment=args.environment)
    project_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    logger = get_logger(f"ACTIVE-SPACE: {args.unit}{args.site}{args.year}")

    # # Verify that annotation files exist
    # logger.info("Locating unit/site annotations...")
    # annotation_files = glob.glob(f"{project_dir}/{args.unit}{args.site}*_saved_annotations.geojson")
    # if len(annotation_files) == 0:
    #     logger.info(f"No track annotations found for {args.unit}{args.site}{args.year}. Exiting...")
    #     exit(-1)
    #
    # # If annotations do exist, load them into memory.
    # annotations = Annotations()
    # for file in tqdm(annotation_files, desc='Loading annotation files', unit='files', colour='blue'):
    #     annotations = annotations.append(Annotations(file), ignore_index=True)

    # TODO: altutide calucalte

    # Load the microphone deployment site metadata and the study area shapefile.
    microphone = get_deployment(args.unit, args.site, args.year, cfg.read('data', 'site_metadata'), elevation=False)
    study_area = gpd.read_file(glob.glob(f"{project_dir}/*study*.shp")[0])

    # Load NVSPL data or the mennitt raster depending on the user input.
    if args.ambience == 'nvspl':
        archive = iyore.Dataset(cfg.read('data', 'nvspl_archive'))
        nvspl_files = [e.path for e in archive.nvspl(unit=args.unit, site=args.site, year=str(args.year), n=2)] # TODO remove this....
        ambience = Nvspl(nvspl_files)
    else:
        ambience = cfg.read('data', 'mennitt')

    n_sources = 51
    omni_sources = get_omni_sources(lower=args.omni_min, upper=args.omni_max)
    results = gpd.GeoDataFrame(columns=['f1', 'precision', 'recall'])

    active_space_generator = ActiveSpaceGenerator(
        NMSIM=cfg.read('project', 'nmsim'),
        root_dir=project_dir,
        study_area=study_area,
        ambience_src=ambience,
        dem_src=cfg.read('data', 'dena_dem'),

    )

    logger.info(f"Generating active space for: {args.unit}{args.site}{args.year}...")

    # loop through source file and create an active space for each one
    for i, omni_source in enumerate(omni_sources):

        active_spaces = gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs='epsg:4269')

        for heading in args.headings:

            active_space = active_space_generator.generate(
                omni_source=omni_source,
                mic=microphone,
                heading=heading
            )
            active_spaces = active_spaces.append(active_space, ignore_index=True)

        dissolved_active_space = active_spaces.dissolve()

        dissolved_active_space.to_file(f'C:/Users/azucker/Desktop/{args.unit}{args.site}{args.year}_{os.path.basename(omni_source)}.geojson', driver='GeoJSON', mode='w', index=False)  # TODO: MOVE OUT ONE PLACE and change file name

    #
    #     f1, precision, recall, n_tot = compute_f1(valid_points, active_space)
    #     active_space['f1'] = f1
    #     active_space['precision'] = precision
    #     active_space['recall'] = recall
    #     active_space.index = [row[1].gain]  # index
    #
    #     results = results.append(active_space)
    #
    #     print(f"\n\t[{i + 1}/{n_sources}] Finished active space for source", row[1].filename)
    #     print(
    #         "\tF1 score: {:1.2}, Precision: {:1.2}, Recall: {:1.2} for {:0} points.\n".format(f1, precision,
    #                                                                                           recall, n_tot))
    #
    # # this needs to come after there is at least one row
    # results.index.name = 'gain'
    # results.to_wkt().to_csv(PROJ_DIR + os.sep + u + s + '_active_spaces_3octband_{:.0f}ft.txt'.format(alt_ft))
    #
    # endtime = time.perf_counter()
    #
    # print("\nFINISHED COMPUTING ACTIVE SPACE, average time {:.2f} s per polygon".format(
    #     (endtime - starttime) / n_sources))




