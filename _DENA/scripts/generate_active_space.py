import glob
import multiprocessing as mp
import os
import numpy as np
from argparse import ArgumentParser
from copy import deepcopy
from functools import partial
from pathlib import Path
from typing import List, Tuple, TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point
from tqdm import tqdm

# for some users relative imports are prohibitive
# we simplify imports by adding three directories to the path environment variable
import sys
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
config_dir = os.path.join(repo_dir, "_DENA")
script_dir = os.path.join(repo_dir, "nps_active_space")
sys.path.append(repo_dir)
sys.path.append(config_dir)
sys.path.append(script_dir)

import iyore

import _DENA.resource.config as cfg
from _DENA import DENA_DIR
from _DENA.resource.helpers import get_deployment, get_logger, get_omni_sources, load_annotations
from nps_active_space.utils import Annotations, Nvspl
from nps_active_space.utils.computation import select_optimal, ambience_from_nvspl, ambience_from_raster, normalize_point_density
from nps_active_space.active_space import ActiveSpaceGenerator

if TYPE_CHECKING:
    from nps_active_space.utils import Microphone


# Callback functions for multiprocessing.
_handle_error = lambda error: print(f'Error: {error}', flush=True)
_update_pbar = lambda _: pbar.update()


def _run_active_space(outfile: str, omni_source: str, generator: ActiveSpaceGenerator, headings: List[int],
                      microphone: 'Microphone', altitude: int) -> Tuple[str, gpd.GeoDataFrame]:
    """
    Function to be multiprocessed to generate active spaces for multiple omni sources.

    Parameters
    ----------
    outfile : str
        Name of the file where the final active space should be output.
    omni_source : str
        Tuning source to generate the active space with.
    generator : ActiveSpaceGenerator
        The active space generator to use to create active spaces.
    headings : List[int]
        List of directional headings to generate active spaces for. Active spaces for each heading will be dissolved
        to create a single all encompassing active space.
    microphone : Microphone
        Location to generate the active space around.
    altitude : int
        Altitude (in meters) to generate the active space at.

    Returns
    -------
    omni_source : str
        The path to the omni source file that was used to create the active space.
    dissolved_active_space : gpd.GeoDataFrame
        The final generated active space for the given parameters.
    """
    # NOTE: Since the microphone is being used in multiple processes and in those processes is altered, it's safer to
    #  make copies of the microphone with unique names to avoid any issues with shared resources.
    mic_copy = deepcopy(microphone)
    mic_copy.name = f"{microphone.name}{Path(omni_source).stem}"

    active_spaces = None
    for heading in headings:
        active_space = generator.generate(
            omni_source=omni_source,
            mic=mic_copy,
            heading=heading,
            altitude_m=altitude
        )

        if active_spaces is None:
            active_spaces = active_space
        else:
            active_spaces = pd.concat([active_spaces, active_space], ignore_index=True)

    # Combine the active spaces from each heading into a single active space and write it to a geojson file.
    dissolved_active_space = active_spaces.dissolve()
    dissolved_active_space.to_file(outfile, driver='GeoJSON', mode='w', index=False)

    return Path(omni_source).stem, dissolved_active_space


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
    argparse.add_argument('-a', '--ambience', default='nvspl',
                          help="What type of ambience to use in NMSIM calculations. Choose from ['nvspl', 'mennitt', or a path to an ambience .pkl file]")
    argparse.add_argument('--headings', nargs='+', type=int, default=[0, 120, 240],
                          help="Headings of active spaces to dissolve. Accepts one or more values.")
    argparse.add_argument('--omni-min', type=float, default=-20,
                          help="The minimum omni source to run the mesh for.")
    argparse.add_argument('--omni-max', type=float, default=30,
                          help="The maximum omni source to run the mesh for.")
    argparse.add_argument('-l', '--altitude', type=int, required=False,
                          help="Altitude to run NSMIM with in meters.")
    argparse.add_argument('-b', '--beta', nargs='+', type=float, default=[1.0, 2.0], 
                          help="Beta value(s) to use when calculating fbeta. Accepts one or more values.")
    argparse.add_argument('--cleanup', action='store_true',
                          help="Remove intermediary control and batch files.")
    args = argparse.parse_args()

    ambience_valid = (args.ambience == "nvspl") or (args.ambience == "mennitt") or (args.ambience.endswith(".pkl") and os.path.exists(args.ambience))
    assert ambience_valid, "Ambience argument must be 'nvspl', 'mennitt', or a .pkl file"

    # --------------- INIT --------------- #

    cfg.initialize(f"{DENA_DIR}/config", environment=args.environment)
    project_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    logger = get_logger(f"ACTIVE-SPACE: {args.unit}{args.site}{args.year}")

    omni_sources = get_omni_sources(lower=args.omni_min, upper=args.omni_max)

    # --------------- DATA SELECTION --------------- #

    # Load the microphone deployment site metadata and the study area shapefile.
    mic_ = get_deployment(cfg.read('project', 'dir'), args.unit, args.site, args.year, elevation=False)
    study_area = gpd.read_file(glob.glob(f"{project_dir}/*study*.shp")[0])

    # Compute ambience
    # Load NVSPL data or the mennitt raster depending on the user input.
    if args.ambience == 'nvspl':
        archive = iyore.Dataset(cfg.read('data', 'nvspl_archive'))
        nvspl_files = [e.path for e in archive.nvspl(unit=args.unit, site=args.site, year=str(args.year))]
        nvspl = Nvspl(nvspl_files)
        ambience_quantile = 90  # L90 = 90% exceedance = 10% quantile sound level
        ambience = ambience_from_nvspl(nvspl, ambience_quantile, broadband=False)
    elif args.ambience == 'mennitt':
        ambience = ambience_from_raster(cfg.read('data', 'mennitt'), mic_)
    else:
        # should be a .pkl filename
        ambience = pd.read_pickle(args.ambience)
        print(f"Read ambience from {args.ambience}")

    # --------------- ANNOTATION LOGIC --------------- #

    # Verify that annotation files exist for the unit/site location. If they do exist, load them into memory.
    logger.info("Locating unit/site annotations...")
    annotations = load_annotations(cfg.read("project", "dir"), args.unit, args.site, args.year)
    if annotations.empty:
        logger.info(f"No track annotations found for {args.unit}{args.site}{args.year}. Exiting...")
        exit(-1)
    print(f"{annotations.shape[0]} valid annotated segments found")

    # Extract all valid points from their LineStrings. These will be needed for calculating fbeta scores later.
    valid_points_lst = []
    for idx, row in tqdm(annotations.iterrows(), total=annotations.shape[0], desc='Extracting valid points', unit='valid track', colour='white'):
        valid_points_lst.extend([{'annotation_idx': idx, 'audible': row.audible, 'geometry': Point(coords)} for coords in row.geometry.coords])
    valid_points = gpd.GeoDataFrame(data=valid_points_lst, geometry='geometry', crs=annotations.crs)

    # Reduce point density to median density, so very dense areas (e.g. airports) don't skew the fit
    valid_points = normalize_point_density(valid_points, study_area, random_seed=679)

    # If the user does not pass an altitude, calculate the average altitude of all valid tracks. Extract the altitudes
    #  from each linestring to get the average height (in meters) of audible flight segments.
    # We just use the audible segments so that we represent the typical altitude in the local area.
    #  (some inaudible segments are very far away / at different altitudes)
    if not args.altitude:
        logger.info("Calculating average altitude (in meters)...")
        annotation_altitudes = []
        relevant_mask = (valid_points["audible"]) & (valid_points.geometry.z > 0) & (valid_points.geometry.z <= 10000)
        # NOTE we only apply the altitude filter for the purposes of calculating mean altitude, instead of removing
        # them altogether, because removing the negative values could be severe for some ADS-B loggers
        for idx, group in valid_points.loc[relevant_mask].groupby("annotation_idx"):
            annotation_altitudes.append(group.geometry.z.mean())
        altitude_ = int(np.mean(annotation_altitudes))
        # annotations['z_val'] = (annotations['geometry'].apply(lambda geom: mean([coords[-1] for coords in geom.coords])))
        # altitudes_ = annotations[annotations.audible == True].z_val
        # altitudes_ = altitudes_[(altitudes_ > 0)&(altitudes_ <= 10000)] # NOTE removing the negative values could be severe for some ADS-B loggers
        # altitude_ = int(mean(altitudes_.tolist()))
        logger.info(f"Average altitude is: {altitude_}m")
    else:
        altitude_ = args.altitude

    # --------------- ACTIVE SPACE GENERATION --------------- #

    # Create an ActiveSpaceGenerator instance and set the DEM data for the microphone location since we will be using
    #  the same location for every active space. This is a MAJOR time saver!
    generator_ = ActiveSpaceGenerator(
        NMSIM=cfg.read('project', 'nmsim'),
        root_dir=project_dir,
        study_area=study_area,
        ambience=ambience,
        dem_src=cfg.read('data', 'dem'),
    )
    logger.info('Setting dem...')
    generator_.set_dem(mic_)

    # Create active space for each omni source. Active spaces are created in parallel asynchronously for maximum
    #  speed benefits.
    logger.info(f"Generating active spaces for: {args.unit}{args.site}{args.year}...")
    _run = partial(_run_active_space, generator=generator_, headings=args.headings, microphone=mic_, altitude=altitude_)
    with mp.Pool(mp.cpu_count() - 1) as pool:
        with tqdm(desc='Omni Sources', unit='omni source', colour='green', total=len(omni_sources), leave=True) as pbar:
            processes = []
            for omni_source_ in omni_sources:
                outfile_ = f'{project_dir}/{args.unit}{args.site}{args.year}_{Path(omni_source_).stem}.geojson'
                processes.append(pool.apply_async(_run, kwds={'outfile': outfile_, 'omni_source': omni_source_},
                                                  callback=_update_pbar, error_callback=_handle_error))
            results = [p.get() for p in processes]

    valid_results = [result for result in results if result is not None]

    # Clean up intermediary files if the user requests.
    if args.cleanup:
        for file in glob.glob(f"{project_dir}/control*"):
            os.remove(file)
        for file in glob.glob(f"{project_dir}/batch*"):
            os.remove(file)

    # --------------- ANALYSIS --------------- #

    for beta_ in args.beta:

        plot_savepath = f'{project_dir}/PrecisionRecallPlot_{args.unit}{args.site}{args.year}_{str(beta_).replace(".","p")}.png'
        best_omni, max_fbeta, best_precision, best_recall, detection_results = select_optimal(unit=args.unit,
                                                                                              site=args.site,
                                                                                              year=args.year,
                                                                                              valid_points=valid_points,
                                                                                              active_space_polygons=results,
                                                                                              beta_=beta_,
                                                                                              plot=True,
                                                                                              plot_savepath=plot_savepath)

        logger.info(f"The best performing omni source for F-{beta_} is: {best_omni} (fbeta: {max_fbeta})")
