import logging
import multiprocessing as mp
import os
import re
import signal
import numpy as np
from argparse import ArgumentParser
from copy import deepcopy
from functools import partial
from pathlib import Path
from typing import List, Optional, Tuple, TYPE_CHECKING
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from tqdm import tqdm
import pickle
import sys

# Contour extraction uses matplotlib inside worker processes. TkAgg on Windows
# raises "main thread is not in main loop" / Tcl_AsyncDelete during figure teardown.
os.environ.setdefault("MPLBACKEND", "Agg")

import nps_active_space.utils.config as cfg
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.utils.helpers import get_deployment, get_logger, get_omni_sources, load_annotations, omni_to_gain
from nps_active_space.utils import paths as p
from nps_active_space.utils.models import Annotations
from nps_active_space.utils.computation import (
    select_optimal,
    compute_ambience_from_nvspl_archive,
    ambience_from_raster,
    normalize_point_density,
    load_spectral_ambience_pickle,
)
from nps_active_space.active_space import ActiveSpaceGenerator
from nps_active_space.active_space.active_space_setup import (
    DEFAULT_SRC_PT_DENSITY,
    build_active_space_generator,
    build_batch_run_results,
    cleanup_propagation_artifacts,
    resolve_acoustic_model,
    write_batch_run_results,
)

if TYPE_CHECKING:
    from nps_active_space.utils.models import Microphone

OMNI_GAIN_STEM_RE = re.compile(r"O_[+-]\d{3}$")


def _run_active_space(outfile: str, omni_source: str, generator: ActiveSpaceGenerator, headings: List[int],
                      microphone: 'Microphone', altitude: int, src_pt_density: int,
                      tested_pts_outfile: Optional[str] = None,
                      pretested_pts_dict : Optional[dict] = None) -> Tuple[str, gpd.GeoDataFrame, dict]:
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
    tested_pts_outfile : str, default None
        If provided, name of the .pkl file to save tested points to.
    pretested_pts_dict : dict, default None
        Dictionary storing points we know are audible/inaudible already. This can happen when a quieter source
        is determined to be audible somewhere; a louder source will still be audible there.
        Keys are headings, values are GeoDataFrames of 3D points with an "audible" field = 0 or 1.

    Returns
    -------
    omni_source : str
        The path to the omni source file that was used to create the active space.
    dissolved_active_space : gpd.GeoDataFrame
        The final generated active space for the given parameters.
    tested_pts_dict : dict
        Dictionary containing the points that were tested for audibility.
        Keys are headings, values are GeoDataFrames of 3D points with an "audible" field = 0 or 1.
    """
    assert outfile.endswith(".geojson")
    assert tested_pts_outfile is None or tested_pts_outfile.endswith(".pkl")

    # NOTE: Since the microphone is being used in multiple processes and in those processes is altered, it's safer to
    #  make copies of the microphone with unique names to avoid any issues with shared resources.
    mic_copy = deepcopy(microphone)
    mic_copy.name = f"{microphone.name}{Path(omni_source).stem}"

    active_space_list = []
    tested_pts_dict = {}
    for heading in headings:
        # get the predetermined points corresponding to this heading if they exist
        predetermined_audibility_pts = None if pretested_pts_dict is None else pretested_pts_dict[heading]

        active_space, tested_pts = generator.generate(
            omni_source=omni_source,
            mic=mic_copy,
            heading=heading,
            altitude_m=altitude,
            src_pt_density=src_pt_density,
            predetermined_audibility_pts=predetermined_audibility_pts
        )
        active_space_list.append(active_space)
        tested_pts_dict[heading] = tested_pts

    # Combine the active spaces from each heading into a single active space and write it to a geojson file.
    active_spaces = pd.concat(active_space_list, ignore_index=True)
    dissolved_active_space = active_spaces.dissolve()
    dissolved_active_space.to_file(outfile, driver='GeoJSON', mode='w', index=False)

    # Save the points we sampled and whether they were audible. This is useful for debugging and presentations.
    if tested_pts_outfile is not None:
        with open(tested_pts_outfile, "wb") as f:
            pickle.dump(tested_pts_dict, f)

    return Path(omni_source).stem, dissolved_active_space, tested_pts_dict


def group_omni_sources(omnis: List[str]) -> List[List[str]]:
    """
    Sort and order omni sources into groups, to get maximum benefit from get_pretested_pts().

    Rationale:
    Ideally, we want previous gains to closely bound new gains, so that most of the space can
    be predetermined as audible or inaudible based on previous results (see get_pretested_pts()).
    So instead of calculating gains in increasing order, it makes sense to first calculate the middle gain,
    then the gains at the middle of the upper half and lower half of gains, etc. This quickly ensures
    that any future gain we want to test will have close lower and upper bounds.

    To make this work with multiprocessing, we compute active space gains in groups, where each group
    is multiprocessed. Future groups can make use of past groups' results. So to follow the above rationale,
    the first group is the gain in the middle, the next group is the two gains dividing the upper/lower halves, etc.
    So group size goes 1, 2, 4, 8, etc. However, starting with very small groups doesn't utilize the CPU well
    when multiprocessing, and this is more important than benefitting from previous results. So we collapse the
    first several groups into a single group to improve this.

    Parameters
    ----------
    omnis: List[str]
        List of absolute paths to omni source files, with names like O_+125.src
    
    Returns
    -------
    groups: List[List[str]]
        A list of omni source groups, in the order we should compute them. Each group is a list of
        absolute paths to omni source files.
    """
    # store omni source filename indexed by the corresponding gain
    gains = list(map(omni_to_gain, omnis))
    df = pd.DataFrame({"omni": omnis}, index=gains)
    df = df.sort_index()

    # define recursive function that will look for the halfway gain recursively and assign group numbers
    # the group number will be equivalent to the recursion depth; this matches the rationale in the docstring
    def recurse(idxs, group=0):
        if not idxs:
            return
        mid = len(idxs) // 2
        df.loc[idxs[mid], "group"] = group
        recurse(idxs[:mid], group+1)   # left half
        recurse(idxs[mid+1:], group+1) # right half

    # run recursion on the list of gain indices
    recurse(df.index.tolist())

    # collapse the early groups together so that we can make good use of the CPU
    # in the first round of multiprocessing. 3 is empirically chosen to utilize the CPU well.
    df.loc[df["group"] < 3, "group"] = 3

    # create list of groups
    groups = []
    for _, df_group in df.groupby("group"):
        groups.append(df_group["omni"].tolist())

    return groups


def get_pretested_pts(tested_pts_record: dict, gain: float, headings: List[int]) -> dict:
    """
    Get points that we already know audibility/inaudibility for, given previous results.
    This is based on the principle that if a point is audible and you increase the gain,
    it will remain audible. Similarly, if a point is inaudible and you decrease the gain,
    it will still be inaudible.

    Parameters
    ----------
    tested_pts_record: dict
        Dictionary storing the history of tested points and their audibility.
        Keys are gains, values are the dictionary "tested_pts_dict" returned by _run_active_space(),
        which has keys for each heading, and values are GeoDataFrames of 3D points with an "audible" field = 0 or 1.
    gain: float
        The active space gain we wish to compute next, that we are interested in pretested points for.
    headings: List[int]
        A list of headings that we want pretested points for.
    
    Returns
    -------
    pts_dict: dict
        A dictionary containing points with their expected audibility based on past results.
        Keys are headings, values are GeoDataFrames of 3D points with an "audible" field = 0 or 1.
    """
    # search for previous gain(s) that most closely lower and upper bound this gain
    # these will restrict the space we need to test the most, providing maximum speed benefits

    prev_gains = np.array(list(tested_pts_record.keys()))
    
    smaller_gain = None
    smaller_gains = prev_gains[prev_gains < gain]
    if smaller_gains.shape[0] > 0:
        smaller_gain = smaller_gains.max()
    
    larger_gain = None
    larger_gains = prev_gains[prev_gains > gain]
    if larger_gains.shape[0] > 0:
        larger_gain = larger_gains.min()

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

    # combine audible and inaudible points into a single GeoDataFrame for each heading
    for h in pts_dict:
        pts_dict[h] = pd.concat(pts_dict[h], axis=0, ignore_index=True)
    
    return pts_dict


def init_worker():
    """Worker initializer to allow clean Ctrl+C of multiprocessing.
    This makes workers ignore Ctrl+C so that pool.terminate() can take care
    of cleanly terminating the multiprocess workers."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def resolve_pool_n_workers(model: AcousticModel) -> int:
    """Omni-source pool size. Docker+Wine caps this via ``AAM_PARALLEL_N`` in run_activespace.sh."""
    max_workers = max(1, mp.cpu_count() - 1)
    if AcousticModel.parse(model) is not AcousticModel.AAM:
        return max_workers
    raw = os.environ.get("AAM_PARALLEL_N")
    if raw is None:
        return max_workers
    return max(1, min(int(raw), max_workers))


def _nonempty_active_space_count(results: list[tuple[str, gpd.GeoDataFrame]]) -> int:
    """Count generated active space layers that contain at least one non-empty geometry."""
    count = 0
    for _, active_space in results:
        if active_space.empty:
            continue
        geometries = active_space.geometry
        if geometries.notna().any() and (~geometries.is_empty).any():
            count += 1
    return count


def _fail_active_space_generation(message: str) -> None:
    print(message, flush=True)
    logging.getLogger(__name__).error(message)
    sys.exit(1)


def build_parser() -> ArgumentParser:
    parser = ArgumentParser()

    parser.add_argument('-e', '--environment', required=True,
                        help="The configuration environment to run the script in.")
    parser.add_argument('-u', '--unit', required=True,
                        help="Four letter unit code. E.g. DENA")
    parser.add_argument('-s', '--site', required=True,
                        help="Four letter site code. E.g. TRLA")
    parser.add_argument('-y', '--year', type=int, required=True,
                        help="Four digit year. E.g. 2018")
    parser.add_argument('-a', '--ambience', default='nvspl',
                        help="Ambience for audibility: 'nvspl', 'mennitt', or a path to an ambience .pkl file")
    parser.add_argument('--model', type=AcousticModel, choices=list(AcousticModel),
                        default=None,
                        help="Propagation model (default: nmsim). Docker -m aam also sets this.")
    parser.add_argument('--headings', nargs='+', type=int, default=[0, 120, 240],
                        help="Headings of active spaces to dissolve. Accepts one or more values.")
    parser.add_argument(
        '--source',
        action='append',
        dest='sources',
        metavar='PATH',
        help="Additional source path (.src for both models, .nc for AAM only). Repeatable.",
    )
    parser.add_argument('--omni-min', type=float, default=-10,
                        help="The minimum omni source to run the mesh for.")
    parser.add_argument('--omni-max', type=float, default=40,
                        help="The maximum omni source to run the mesh for.")
    parser.add_argument(
        '--omni-step',
        type=float,
        default=0.5,
        help="Spacing between omni gains in dB (multiple of 0.5). Default 0.5 matches full NMSim ladder.",
    )
    parser.add_argument('-l', '--altitude', type=int, required=False,
                        help="Source altitude in meters (default: mean audible annotation altitude).")
    parser.add_argument('--density', type=int, default=None,
                        help="Source-point mesh density (default: pipeline default).")
    parser.add_argument('-b', '--beta', nargs='+', type=float, default=[1.0],
                        help="Beta value(s) to use when calculating fbeta. Accepts one or more values.")
    parser.add_argument(
        '--cleanup-nmsim-scratch',
        '--cleanup',
        action='store_true',
        dest='cleanup_nmsim_scratch',
        help=(
            "After the run, delete NMSim scratch only: control*/batch* at the site root, "
            ".trj under Input_Data/03_TRAJECTORY, and .tis under Output_Data/nmsim/scratch. "
            "No effect on AAM. Does not remove prediction CSV caches or ACTIVESPACES geojson."
        ),
    )
    parser.add_argument('--annotation-file',
                        help="Basename of GEOJSON annotations file to use, if not the default. File should be in the site directory.")
    parser.add_argument(
        '--results-out',
        help="Path to write structured run results as JSON (used by generate_active_space_batch.py).",
    )
    return parser


def _process_omni_group(
    group: list[str],
    run_fn,
    tested_pts_record: dict,
    headings: list[int],
    usy: str,
    active_savedir: str,
    tested_pts_savedir: str,
    results: list,
) -> None:
    processes = []
    for omni_source_ in group:
        stem = Path(omni_source_).stem
        pretested_pts_dict = None
        if OMNI_GAIN_STEM_RE.fullmatch(stem):
            pretested_pts_dict = get_pretested_pts(
                tested_pts_record,
                omni_to_gain(omni_source_),
                headings,
            )
        name = f"{usy}_{stem}"
        kwds = {
            'omni_source': omni_source_,
            'outfile': f'{active_savedir}/{name}.geojson',
            'tested_pts_outfile': f'{tested_pts_savedir}/{name}.pkl',
            'pretested_pts_dict': pretested_pts_dict,
        }
        processes.append(run_fn(kwds=kwds))

    outputs = [p.get() for p in processes]
    for output in outputs:
        if output is None:
            continue
        omni, active, tested_pts_dict = output
        results.append((omni, active))
        if OMNI_GAIN_STEM_RE.fullmatch(omni):
            tested_pts_record[omni_to_gain(omni)] = tested_pts_dict


if __name__ == '__main__':

    args = build_parser().parse_args()
    model = resolve_acoustic_model(args.model)

    ambience_valid = (
        args.ambience in {"nvspl", "mennitt"}
        or args.ambience.endswith(".pkl")
    )
    assert ambience_valid, "Ambience argument must be 'nvspl', 'mennitt', or a .pkl file"

    # --------------- INIT --------------- #

    cfg.initialize(environment=args.environment)
    layout = p.SiteModelPaths.from_project(
        cfg.read('project', 'dir'),
        args.unit,
        args.site,
        args.year,
        model,
    )
    site_dir = layout.site_dir
    logger = get_logger(f"ACTIVE-SPACE: {layout.usy}")

    ladder_sources = get_omni_sources(
        lower=args.omni_min, upper=args.omni_max, step_db=args.omni_step,
    )
    extra_sources = args.sources or []
    for src in extra_sources:
        if Path(src).suffix.lower() == ".nc" and model is not AcousticModel.AAM:
            _fail_active_space_generation(
                f"--source {src}: .nc sources are AAM-only (use --model aam)",
            )
    omni_sources = ladder_sources + extra_sources

    # --------------- DATA SELECTION --------------- #

    # Load the microphone deployment site metadata and the study area shapefile.
    mic_ = get_deployment(cfg.read('project', 'dir'), args.unit, args.site, args.year, elevation=False)
    study_area = gpd.read_file(p.study_area_shapefile(cfg.read('project', 'dir'), args.unit, args.site))

    # Compute ambience
    # Load NVSPL data or the mennitt raster depending on the user input.
    if args.ambience == 'nvspl':
        ambience_quantile = 90  # L90 = 90% exceedance = 10% quantile sound level
        ambience = compute_ambience_from_nvspl_archive(
            cfg.read('data', 'nvspl_archive'),
            args.unit,
            args.site,
            args.year,
            ambience_quantile,
            broadband=False,
        )
    elif args.ambience == 'mennitt':
        ambience = ambience_from_raster(cfg.read('data', 'mennitt'), mic_)
    else:
        ambience = load_spectral_ambience_pickle(args.ambience)
        if ambience is None:
            _fail_active_space_generation(
                f"Ambience pickle is missing or has no usable spectral bands: {args.ambience}"
            )
        print(f"Read ambience from {p.display_path(args.ambience)}")

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
    points_before_kde = len(valid_points)
    valid_points = normalize_point_density(valid_points, study_area, random_seed=679)
    points_after_kde = len(valid_points)

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

    usy = layout.usy
    src_pt_density = args.density if args.density is not None else DEFAULT_SRC_PT_DENSITY
    logger.info(
        f"Generating active spaces for {usy} using {model} "
        f"at {altitude_}m (density={src_pt_density})..."
    )

    logger.info("Caching project_setup elevation...")
    try:
        generator_ = build_active_space_generator(
            site_dir,
            study_area,
            ambience,
            mic_,
            model,
        )
    except FileNotFoundError as exc:
        _fail_active_space_generation(str(exc))

    active_savedir = layout.layer_dir(altitude_)
    tested_pts_savedir = layout.tested_points_dir(altitude_)
    os.makedirs(active_savedir, exist_ok=True)
    os.makedirs(tested_pts_savedir, exist_ok=True)
    os.makedirs(layout.precision_recall_dir, exist_ok=True)

    results = []
    tested_pts_record = {}
    _run = partial(
        _run_active_space,
        generator=generator_,
        headings=args.headings,
        microphone=mic_,
        altitude=altitude_,
        src_pt_density=src_pt_density,
    )

    if extra_sources:
        ladder_groups = group_omni_sources(ladder_sources)
        extra_groups = [[src] for src in extra_sources]
        omni_groups = ladder_groups + extra_groups
    else:
        omni_groups = group_omni_sources(ladder_sources)
    with tqdm(desc='Omni Sources', unit='omni source', colour='green', total=len(omni_sources)) as pbar:
        try:
            pool = mp.Pool(resolve_pool_n_workers(model), init_worker)
            try:
                with pool:
                    _handle_error = lambda error: print(f'Error: {error}', flush=True)
                    _update_pbar = lambda _: pbar.update()
                    for group in omni_groups:
                        _process_omni_group(
                            group,
                            lambda kwds: pool.apply_async(
                                _run,
                                callback=_update_pbar,
                                error_callback=_handle_error,
                                kwds=kwds,
                            ),
                            tested_pts_record,
                            args.headings,
                            usy,
                            active_savedir,
                            tested_pts_savedir,
                            results,
                        )
            except KeyboardInterrupt:
                pool.terminate()
                pool.join()
                raise
        except KeyboardInterrupt:
            if args.cleanup_nmsim_scratch:
                cleanup_propagation_artifacts(site_dir, model)
            sys.exit(1)

    if args.cleanup_nmsim_scratch:
        cleanup_propagation_artifacts(site_dir, model)

    # --------------- ANALYSIS --------------- #

    if not results:
        _fail_active_space_generation(
            "No active space layers were generated successfully. "
            "Check worker errors above, model configuration, and site inputs under "
            f"{p.display_path(site_dir)}."
        )

    nonempty_layers = _nonempty_active_space_count(results)
    if nonempty_layers == 0:
        _fail_active_space_generation(
            f"Active space generation finished but all {len(results)} geojson layers are empty. "
            "The model likely produced no audible source points. Check DEM elevation files under "
            f"{p.display_path(site_dir)}/Input_Data/01_ELEVATION before continuing the 3D workflow."
        )

    best_omni_for_results: str | None = None
    f1_for_results: float | None = None

    for beta_ in args.beta:
        plot_savepath = layout.precision_recall_plot(altitude_, beta_)
        os.makedirs(os.path.dirname(plot_savepath), exist_ok=True)
        best_omni, max_fbeta, best_precision, best_recall, _ = select_optimal(
            unit=args.unit,
            site=args.site,
            year=args.year,
            valid_points=valid_points,
            active_space_polygons=results,
            beta_=beta_,
            plot=True,
            plot_savepath=plot_savepath,
            verbose=False,
        )

        logger.info(
            f"The best performing omni source for F-{beta_} is: {best_omni} (fbeta: {max_fbeta})"
        )
        logger.info("PR plot -> %s", p.display_path(plot_savepath))

        if beta_ == 1.0:
            best_omni_for_results = best_omni
            f1_for_results = max_fbeta

    if args.results_out is not None:
        run_results = build_batch_run_results(
            len(annotations),
            altitude_,
            points_before_kde,
            points_after_kde,
            best_omni=best_omni_for_results,
            f1=f1_for_results,
        )
        write_batch_run_results(args.results_out, run_results)
