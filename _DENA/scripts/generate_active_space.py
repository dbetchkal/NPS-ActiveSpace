import glob
import multiprocessing as mp
import os
import re
import numpy as np
from argparse import ArgumentParser
from copy import deepcopy
from functools import partial
from pathlib import Path
from typing import List, Optional, Tuple, TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point
from tqdm import tqdm
import pickle

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


def _run_active_space(outfile: str, omni_source: str, generator: ActiveSpaceGenerator, headings: List[int],
                      microphone: 'Microphone', altitude: int, pretested_pts_dict : Optional[dict] = None
                      ) -> Tuple[str, gpd.GeoDataFrame, dict]:
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
    pretested_pts_dict : dict, default None
        Dictionary storing points we know are audible/inaudible already. This can happen when a quieter source
        is determined to be audible somewhere; a louder source will still be audible there.
        Keys are headings, values are GeoDataFrames of 3D points with a boolean "audible" field.

    Returns
    -------
    omni_source : str
        The path to the omni source file that was used to create the active space.
    dissolved_active_space : gpd.GeoDataFrame
        The final generated active space for the given parameters.
    tested_pts_dict : dict
        Dictionary containing the points that were tested for audibility.
        Keys are headings, values are GeoDataFrames of 3D points with a boolean "audible" field.
    """
    # NOTE: Since the microphone is being used in multiple processes and in those processes is altered, it's safer to
    #  make copies of the microphone with unique names to avoid any issues with shared resources.
    mic_copy = deepcopy(microphone)
    mic_copy.name = f"{microphone.name}{Path(omni_source).stem}"

    active_space_list = []
    tested_pts_dict = {}
    for heading in headings:
        predetermined_audibility_pts = None if pretested_pts_dict is None else pretested_pts_dict[heading]

        active_space, tested_pts = generator.generate(
            omni_source=omni_source,
            mic=mic_copy,
            heading=heading,
            altitude_m=altitude,
            predetermined_audibility_pts=predetermined_audibility_pts
        )
        active_space_list.append(active_space)
        tested_pts_dict[heading] = tested_pts

    active_spaces = pd.concat(active_space_list, ignore_index=True)

    # Combine the active spaces from each heading into a single active space and write it to a geojson file.
    dissolved_active_space = active_spaces.dissolve()
    dissolved_active_space.to_file(outfile, driver='GeoJSON', mode='w', index=False)

    tested_pts_file = os.path.join(os.path.dirname(outfile), f"tested_pts_{Path(omni_source).stem}.pkl")
    with open(tested_pts_file, "wb") as f:
        pickle.dump(tested_pts_dict, f)

    return Path(omni_source).stem, dissolved_active_space, tested_pts_dict


def omni_to_gain(omni_source):
    """
    TODO - document
    """
    match = re.search(r"O_(....)", omni_source)
    return int(match.group(1)) / 10


def group_omni_sources(omnis):
    gains = list(map(omni_to_gain, omnis))
    df = pd.DataFrame({"omni": omnis}, index=gains)
    df = df.sort_index()

    def recurse(idxs, group=0):
        if not idxs:
            return
        mid = len(idxs) // 2
        df.loc[idxs[mid], "group"] = group
        recurse(idxs[:mid], group+1)   # left half
        recurse(idxs[mid+1:], group+1) # right half

    recurse(df.index.tolist())

    # collapse the early groups together, multiprocessing will be able to handle it
    df.loc[df["group"] < 2, "group"] = 2

    groups = []
    for _, df_group in df.groupby("group"):
        groups.append(df_group["omni"].tolist())

    return groups


def get_pretested_pts(tested_pts_record, gain, headings):
    """
    TODO - document
    """
    # search for previous gain(s) most closely lower and upper bounding this gain

    prev_gains = np.array(list(tested_pts_record.keys()))
    
    smaller_gain = None
    smaller_gains = prev_gains[prev_gains < gain]
    if smaller_gains.shape[0] > 0:
        smaller_gain = smaller_gains.max()
    
    larger_gain = None
    larger_gains = prev_gains[prev_gains > gain]
    if larger_gains.shape[0] > 0:
        larger_gain = larger_gains.min()

    tqdm.write(f"smaller_gain: {smaller_gain}, larger_gain: {larger_gain}")

    if smaller_gain is None and larger_gain is None:
        return None

    # see if we can predetermine any audibility points
    pts_dict = {h: [] for h in headings}

    # any audible points from a smaller gain will still be audible
    if smaller_gain is not None:
        for h in headings:
            df = tested_pts_record[smaller_gain][h]
            pts_dict[h].append(df[df["audible"] == 1])
    
    # any inaudible points from a larger gain will still be inaudible
    if larger_gain is not None:
        for h in headings:
            df = tested_pts_record[larger_gain][h]
            pts_dict[h].append(df[df["audible"] == 0])

    for h in pts_dict:
        pts_dict[h] = pd.concat(pts_dict[h], axis=0, ignore_index=True)
    
    return pts_dict


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
    argparse.add_argument('--annotation-file', help="Basename of GEOJSON annotations file to use, if not the default. File should be in the site directory.")
    args = argparse.parse_args()

    ambience_valid = (args.ambience == "nvspl") or (args.ambience == "mennitt") or (args.ambience.endswith(".pkl") and os.path.exists(args.ambience))
    assert ambience_valid, "Ambience argument must be 'nvspl', 'mennitt', or a .pkl file"

    # --------------- INIT --------------- #

    cfg.initialize(f"{DENA_DIR}/config", environment=args.environment)
    site_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"
    logger = get_logger(f"ACTIVE-SPACE: {args.unit}{args.site}{args.year}")

    omni_sources = get_omni_sources(lower=args.omni_min, upper=args.omni_max)

    # --------------- DATA SELECTION --------------- #

    # Load the microphone deployment site metadata and the study area shapefile.
    mic_ = get_deployment(cfg.read('project', 'dir'), args.unit, args.site, args.year, elevation=False)
    study_area = gpd.read_file(glob.glob(f"{site_dir}/*study*.shp")[0])

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
    if args.annotation_file is not None:
        print(f"Using non-default annotation file: {args.annotation_file}")
        annotations = Annotations(f"{site_dir}/{args.annotation_file}", only_valid=True)
    else:
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
        root_dir=site_dir,
        study_area=study_area,
        ambience=ambience,
        dem_src=cfg.read('data', 'dem'),
    )
    logger.info('Setting dem...')
    generator_.set_dem(mic_)

    # Create active space for each omni source.
    logger.info(f"Generating active spaces for: {args.unit}{args.site}{args.year}...")
    results = []
    tested_pts_record = {}

    _run = partial(_run_active_space, generator=generator_, headings=args.headings, microphone=mic_, altitude=altitude_)
    
    with mp.Pool(mp.cpu_count() - 1) as pool:
        with tqdm(desc='Omni Sources', unit='omni source', colour='green', total=len(omni_sources)) as pbar:
            _handle_error = lambda error: print(f'Error: {error}', flush=True)
            _update_pbar = lambda _: pbar.update()

            for group in group_omni_sources(omni_sources):
                tqdm.write(f"GROUP ({len(group)}) -------------------")
                processes = []
                for omni_source_ in group:
                    gain = omni_to_gain(omni_source_)
                    tqdm.write(f"gain {gain}")
                    outfile_ = f'{site_dir}/{args.unit}{args.site}{args.year}_{Path(omni_source_).stem}.geojson'
                    pretested_pts_dict = get_pretested_pts(tested_pts_record, gain, args.headings)

                    processes.append(pool.apply_async(_run, callback=_update_pbar, error_callback=_handle_error,
                                                      kwds={'outfile': outfile_, 'omni_source': omni_source_,
                                                            'pretested_pts_dict': pretested_pts_dict}))
                outputs = [p.get() for p in processes]

                # process outputs
                for output in outputs:
                    if output is None:
                        continue
                    omni, active, tested_pts_dict = output
                    results.append((omni, active))
                    tested_pts_record[omni_to_gain(omni)] = tested_pts_dict


    # Clean up intermediary files if the user requests.
    if args.cleanup:
        for file in glob.glob(f"{site_dir}/control*"):
            os.remove(file)
        for file in glob.glob(f"{site_dir}/batch*"):
            os.remove(file)

    # --------------- ANALYSIS --------------- #

    for beta_ in args.beta:

        plot_savepath = f'{site_dir}/PrecisionRecallPlot_{args.unit}{args.site}{args.year}_{str(beta_).replace(".","p")}.png'
        best_omni, max_fbeta, best_precision, best_recall, detection_results = select_optimal(unit=args.unit,
                                                                                              site=args.site,
                                                                                              year=args.year,
                                                                                              valid_points=valid_points,
                                                                                              active_space_polygons=results,
                                                                                              beta_=beta_,
                                                                                              plot=True,
                                                                                              plot_savepath=plot_savepath)

        logger.info(f"The best performing omni source for F-{beta_} is: {best_omni} (fbeta: {max_fbeta})")
