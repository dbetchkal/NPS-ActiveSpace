import glob
import os
from argparse import ArgumentParser

import geopandas as gpd
from tqdm import tqdm

import _DENA.resource.config as cfg
from _DENA import DENA_DIR
from _DENA.resource.helpers import get_logger, get_omni_sources
from nps_active_space.active_space import ActiveSpaceGenerator


if __name__ == '__main__':

    argparse = ArgumentParser()

    argparse.add_argument('-e', '--environment', required=True,
                          help="The configuration environment to run the script in.")
    argparse.add_argument('-n', '--name', required=True,
                          help="The name of this mesh. Will be the name of the directory for all associated files.")
    argparse.add_argument('-s', '--study-area', required=True,
                          help='Absolute path to a study area shapefile.')
    argparse.add_argument('--headings', nargs='+', type=int, default=[0, 120, 240],
                          help='Headings of active spaces to dissolve.')
    argparse.add_argument('--omni-min', type=float, default=0,
                          help='The minimum omni source to run the mesh for.')
    argparse.add_argument('--omni-max', type=float, default=0,
                          help='The maximum omni source to run the mesh for.')
    argparse.add_argument('--mesh-spacing', type=int, default=1,
                          help='How far apart in km mesh square center points should be.')
    argparse.add_argument('--mesh-size', type=int, default=25,
                          help='How large in km each mesh square should be. mesh-size x mesh-size.')
    argparse.add_argument('-a', '--altitude', type=int, default=3658,
                          help='Altitude to run NSMIM with in meters.')
    argparse.add_argument('--cleanup', action='store_true',
                          help="Remove intermediary control and batch files.")
    args = argparse.parse_args()

    cfg.initialize(f"{DENA_DIR}/config", environment=args.environment)
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
                    f"mesh spacing: {args.mesh_spacing}km\nmesh size: {args.mesh_size}kmx{args.mesh_size}km\n"
                    f"headings: {args.headings}\noutfile: {outfile}\n")

        for heading in tqdm(args.headings, desc='Heading', unit='headings', colour='cyan', leave=False):
            heading_outfile = f'{project_dir}/{args.name}_{os.path.basename(omni_source)}_{heading}.geojson'
            active_space = active_space_generator.generate_mesh(
                omni_source=omni_source,
                heading=heading,
                altitude_m=args.altitude,
                mesh_density=(args.mesh_spacing, args.mesh_size),
            )
            active_spaces = active_spaces.append(active_space, ignore_index=True)
            active_spaces.to_file(heading_outfile, driver='GeoJSON', mode='w', index=False)

        # Dissolve the active spaces from each heading into one and output to a geojson file.
        dissolved_active_space = active_spaces.dissolve()
        dissolved_active_space.to_file(outfile, driver='GeoJSON', mode='w', index=False)

        # Clean up intermediary files if the user requests.
        if args.cleanup:
            for file in glob.glob(f"{project_dir}/control*"):
                os.remove(file)
            for file in glob.glob(f"{project_dir}/batch*"):
                os.remove(file)
