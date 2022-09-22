import os
from argparse import ArgumentParser

import geopandas as gpd
from tqdm import tqdm

import _DENA.resource.config as cfg
from _DENA import _DENA_DIR
from _DENA.resource.helpers import get_logger, get_omni_sources
from nps_active_space.active_space import ActiveSpaceGenerator


if __name__ == '__main__':

    argparse = ArgumentParser()

    argparse.add_argument('-e', '--environment', required=True,
                          help="The configuration environment to run the script in.")
    argparse.add_argument('-n', '--name', required=True,
                          help="What to name this mesh.")
    argparse.add_argument('-s', '--study-area', required=True,
                          help='Absolute path to a study area shapefile.')
    argparse.add_argument('--headings', nargs='+', type=int, default=[0, 120, 240],
                          help='Headings of active spaces to dissolve.')
    argparse.add_argument('--omni-min', type=float, default=-20,
                          help='The minimum omni source to run the mesh for.')
    argparse.add_argument('--omni-max', type=float, default=30,
                          help='The maximum omni source to run the mesh for.')
    argparse.add_argument('--mesh-spacing', type=int, default=1,
                          help='How far apart in km mesh square center points should be.')
    argparse.add_argument('--mesh-size', type=int, default=25,
                          help='How large in km each mesh square should be. mesh-size x mesh-size.')
    argparse.add_argument('-a', '--altitude', type=int, default=3658,
                          help='Altitude to run NSMIM with in meters.')
    args = argparse.parse_args()

    cfg.initialize(f"{_DENA_DIR}/config", environment=args.environment)
    project_dir = f"{cfg.read('project', 'dir')}/{args.name}"
    logger = get_logger(f"ACTIVE-SPACE: {args.name}")

    # Get the requested omni sources and load the study area shapefile and Mennitt ambience raster.
    omni_sources = get_omni_sources(lower=args.omni_min, upper=args.omni_max)
    study_area = gpd.read_file(args.study_area)
    ambience = cfg.read('data', 'mennitt')

    logger.info('Instantiating ActiveSpaceGenerator...')
    active_space_generator = ActiveSpaceGenerator(
        NMSIM=cfg.read('project', 'nmsim'),
        root_dir=project_dir,
        study_area=study_area,
        ambience_src=ambience,
        dem_src=cfg.read('data', 'dem'),
    )

    logger.info(f"Generating active space mesh for: {args.name}...\n")
    for omni_source in tqdm(omni_sources, desc='Omni Source', unit='omni sources', colour='green'):
        active_spaces = gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs='epsg:4269')
        outfile = f'{project_dir}/{args.name}_{os.path.basename(omni_source)}.geojson'
        logger.info(f"Run attributes:\nomni_source: {os.path.basename(omni_source)}\naltitude (m): {args.altitude}\n"
                    f"mesh spacing: {args.mesh_spacing}\nmesh size: {args.mesh_size}x{args.mesh_size}\n"
                    f"headings: {args.headings}\noutfile: {outfile}\n")
        for heading in tqdm(args.headings, desc='Heading', unit='headings', colour='cyan', leave=False):
            active_space = active_space_generator.generate_mesh(
                omni_source=omni_source,
                heading=heading,
                altitude_m=args.altitude,
                mesh_density=(args.mesh_spacing, args.mesh_size),
            )
            active_spaces = active_spaces.append(active_space, ignore_index=True)
        dissolved_active_space = active_spaces.dissolve()
        dissolved_active_space.to_file(outfile, driver='GeoJSON', mode='w', index=False)
