import glob
import multiprocessing as mp
import os
import pprint
from argparse import ArgumentParser
from copy import deepcopy
from statistics import mean
from typing import List, TYPE_CHECKING

import geopandas as gpd
from shapely.geometry import Point
from tqdm import tqdm

import iyore

import _DENA.resource.config as cfg
from _DENA import DENA_DIR
from _DENA.resource.helpers import get_deployment, get_logger, get_omni_sources
from nps_active_space.utils import Annotations, compute_f1, Nvspl
from nps_active_space.active_space import ActiveSpaceGenerator

if TYPE_CHECKING:
    from nps_active_space.utils import Microphone


def _update_pbar(_):
    pbar.update()


def _run_active_space(generator: ActiveSpaceGenerator, headings: List[int], omni_source: str, microphone: 'Microphone',
                      altitude: int, outfile: str = None) -> gpd.GeoDataFrame:
    """
    Function to be multiprocessed to generate active spaces for mulitple omni sources.

    Parameters
    ----------
    generator : ActiveSpaceGenerator
    headings : List[int]
        List of directional headings to generate active spaces for. Active spaces for each heading will be dissolved
        to create a single all encompassing active space.
    omni_source : str
        Tuning source to generate the active space with.
    microphone : Microphone
        Microphone location to generate the active space around.
    altitude : int
        Altitude (in meters) to generate the active space at.

    # TODO fix outfile

    Returns
    -------
    dissolved_active_space : gpd.GeoDataFrame
        The final generated active space for the given parameters.
    """
    mic = deepcopy(microphone)
    mic.name = f"{microphone.name}{os.path.basename(omni_source)}"
    active_spaces = gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs='epsg:4269')
    for heading in headings:
        active_space = generator.generate(
            omni_source=omni_source,
            mic=mic,
            heading=heading,
            altitude_m=altitude
        )
        active_spaces = active_spaces.append(active_space, ignore_index=True)

    dissolved_active_space = active_spaces.dissolve()
    dissolved_active_space.to_file(outfile, driver='GeoJSON', mode='w', index=False)
    return dissolved_active_space


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
    argparse.add_argument('--omni-min', type=float, default=-20,
                          help='The minimum omni source to run the mesh for.')
    argparse.add_argument('--omni-max', type=float, default=30,
                          help='The maximum omni source to run the mesh for.')
    argparse.add_argument('-l', '--altitude', type=int, required=False,
                          help='Altitude to run NSMIM with in meters.')
    # argparse.add_argument('-n', '--max-tracks', type=int, default=None, # TODO: do we need this...?
    #                       help='Maximum number of tracks to load')
    args = argparse.parse_args()

    cfg.initialize(f"{DENA_DIR}/config", environment=args.environment)
    project_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    logger = get_logger(f"ACTIVE-SPACE: {args.unit}{args.site}{args.year}")

    # Verify that annotation files exist
    logger.info("Locating unit/site annotations...")
    annotation_files = glob.glob(f"{project_dir}/{args.unit}{args.site}*_saved_annotations.geojson")
    if len(annotation_files) == 0:
        logger.info(f"No track annotations found for {args.unit}{args.site}{args.year}. Exiting...")
        exit(-1)

    # If annotations do exist, load them into memory.
    annotations = Annotations()
    for file in tqdm(annotation_files, desc='Loading annotation files', unit='files', colour='blue'):
        annotations = annotations.append(Annotations(file), ignore_index=True)

    annotations = annotations[annotations.valid == True]
    if not args.altitude:
        # Extract the altitudes from each linestring to get the average height (in meters) of audible flight segments.
        logger.info("Calculating average altitude (in meters)...")
        annotations['z_vals'] = (annotations['geometry'].apply(lambda geom: mean([coords[-1] for coords in geom.coords])))
        altitude_ = int(mean(annotations[(annotations.valid == True) & (annotations.audible == True)].z_vals.tolist()))
        logger.info(f"Average altitude is: {altitude_}m")
    else:
        altitude_ = args.altitude

    # Extract the valid points
    logger.info('Building valid points gdf...')
    valid_points_lst = []
    for idx, row in annotations[(annotations.valid == True)].iterrows():
        valid_points_lst.extend([{'audible': row.audible, 'geometry': Point(coords)} for coords in row.geometry.coords]) # TODO is there a better way to do this...
    valid_points = gpd.GeoDataFrame(data=valid_points_lst, geometry='geometry', crs=annotations.crs)

    # Load the microphone deployment site metadata and the study area shapefile.
    microphone_ = get_deployment(args.unit, args.site, args.year, cfg.read('data', 'site_metadata'), elevation=False)
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

    generator_ = ActiveSpaceGenerator(
        NMSIM=cfg.read('project', 'nmsim'),
        root_dir=project_dir,
        study_area=study_area,
        ambience_src=ambience,
        dem_src=cfg.read('data', 'dem'),
    )

    logger.info('Setting dem...')
    generator_.set_dem(microphone_)

    # Create active space for each omni source.
    logger.info(f"Generating active spaces for: {args.unit}{args.site}{args.year}...")
    pbar = tqdm(desc='Omni Sources', unit='omni source', colour='red', total=len(omni_sources), position=0, leave=True)
    with mp.Pool(mp.cpu_count() - 1) as pool:
        processes = []
        for omni_source_ in omni_sources:
            outfile = f'{project_dir}/{args.unit}{args.site}{args.year}_{os.path.basename(omni_source_)}.geojson'
            processes.append(pool.apply_async(_run_active_space,
                                              args=(generator_, args.headings, omni_source_, microphone_, altitude_, outfile),
                                              callback=_update_pbar))
        results = [p.get() for p in processes]

    active_space_scores = {}
    for res in results:
        f1, precision, recall, n_tot = compute_f1(valid_points, res)
        active_space_scores[f1] = {'omni': omni_source_, 'precision': precision, 'recall': recall}
    pprint.pprint(active_space_scores)
    print(max(active_space_scores.keys()), active_space_scores[max(active_space_scores.keys())])
    logger.info('Complete')
