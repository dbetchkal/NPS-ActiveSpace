import glob
import time
import os
from argparse import ArgumentParser

import geopandas as gpd
import pandas as pd

import iyore

import _DENA.resource.config as cfg
from _DENA import _DENA_DIR
from _DENA.resource.helpers import get_deployment, get_logger, get_omni_sources
from nps_active_space.utils import Annotations
from legacy_code.NMSIM-Python import ActiveSpace import CreateActiveSpace


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
    argparse.add_argument('-o', '--outfile', required=False,
                          help='A specific file to output active space polygon to.')
    argparse.add_argument('-n', '--max-tracks', type=int, default=None,
                          help='Maximum number of tracks to load')
    args = argparse.parse_args()

    cfg.initialize(f"{_DENA_DIR}/config", environment=args.environment)
    logger = get_logger(f"ACTIVE-SPACE: {args.unit}{args.site}")
    project_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"

    # Verify that annotation files exist
    logger.info("Locating unit/site annotations...")
    annotation_files = glob.glob(f"{project_dir}/{args.unit}{args.site}*_saved_annotations")
    if len(annotation_files) == 0:
        logger.info(f"No track annotations found for {args.unit}{args.site}{args.year}. Exiting...")
        exit(-1)

    # If annotations do exist, load them into memory.
    years = []
    logger.info(f"Annotations found for: {years}")
    annotations = Annotations()
    for file in annotation_files:
        annotations = pd.concat([annotations, Annotations(file)], ignore_index=True)

    # TODO outfile = args.outfile if args.outfile else f"{project_dir}/{args.unit}{args.site}{args.year}_saved_annotations"

    # Load the microphone deployment site metadata and the study area shapefile.
    microphone = get_deployment(args.unit, args.site, args.year, cfg.read('data', 'site_metadata'))
    study_area = gpd.read_file(glob.glob(f"{project_dir}/*study*.shp")[0])

    n_sources = 51
    source_df = get_omni_sources(lower=-20, upper=30)
    results = gpd.GeoDataFrame(columns=['f1', 'precision', 'recall'])

    active_space_generator = CreateActiveSpace(
        unit_name=args.unit,
        site_name=args.site,
        yr=args.year,
        project_dir=project_dir,
        freq_bands=True
    )
    # , ambience_source='Mennitt')

    logger.info(f"Generating active space for: {args.unit}{args.site}{args.year}...")

    # loop through source file and create an active space for each one
    for i, row in enumerate(source_df.iterrows()):
        active_space_generator.set_source(row[1].full_path)  # set the source

        active_space = active_space_generator.run_model(altitude_ft=alt_ft, n_tracks=1, out_crs=microphone.crs)
        if not active_space:
            print(f"\n\t[{i + 1}/{n_sources}] Skipped active space for source", row[1].filename)
            continue  # noise source was too quiet (empty active space) or the script failed to produce a polygon

        f1, precision, recall, n_tot = compute_f1(valid_points, active_space)
        active_space['f1'] = f1
        active_space['precision'] = precision
        active_space['recall'] = recall
        active_space.index = [row[1].gain]  # index

        results = results.append(active_space)

        print(f"\n\t[{i + 1}/{n_sources}] Finished active space for source", row[1].filename)
        print(
            "\tF1 score: {:1.2}, Precision: {:1.2}, Recall: {:1.2} for {:0} points.\n".format(f1, precision,
                                                                                              recall, n_tot))

    # this needs to come after there is at least one row
    results.index.name = 'gain'
    results.to_wkt().to_csv(PROJ_DIR + os.sep + u + s + '_active_spaces_3octband_{:.0f}ft.txt'.format(alt_ft))

    endtime = time.perf_counter()

    print("\nFINISHED COMPUTING ACTIVE SPACE, average time {:.2f} s per polygon".format(
        (endtime - starttime) / n_sources))




