import glob
import os
from argparse import ArgumentParser

import geopandas as gpd
import sqlalchemy

import iyore
import nps_active_space.ground_truthing as app
from nps_active_space.utils import Nvspl, Tracks

import _DENA.resource.config as cfg
from _DENA import _DENA_DIR
from _DENA.resource.helpers import get_deployment, get_logger, query_adsb, query_tracks


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
    argparse.add_argument('-t', '--track-source', default='Database', choices=["Database", "ADSB", "AIS"],
                          help="Enter 'Database', 'ADSB', or 'AIS")

    args = argparse.parse_args()

    cfg.initialize(f"{_DENA_DIR}/config", environment=args.environment)
    logger = get_logger('GROUND-TRUTHING')
    engine = sqlalchemy.create_engine(
        'postgresql://{username}:{password}@{host}:{port}/{name}'.format(**cfg.read('database:overflights'))
    )

    logger.info(f"Beginning ground truthing process for {args.unit}{args.site}{args.year}...")

    # Set the various path variables.
    archive = iyore.Dataset(cfg.read('data', 'nvspl_archive'))
    project_dir = f"{cfg.read('project', 'dir')}/{args.unit}{args.site}"

    # Load the microphone deployment site metadata and the study area shapefile.
    microphone = get_deployment(args.unit, args.site, args.year, cfg.read('data', 'site_metadata'))
    study_area = gpd.read_file(glob.glob(f"{project_dir}/*study*.shp")[0])

    # Retrieve the days for which at least some NVSPL data exist.
    nvspl_dates = sorted(set([f"{e.year}-{e.month}-{e.day}" for e in archive.nvspl(unit=args.unit, site=args.site, year=args.year)]))

    # Query flight tracks from days there is NVSPL data for.
    logger.info("Querying tracks...")

    if args.track_source == 'ADSB':
        raw_tracks = query_adsb(
            adsb_files=glob.glob(os.path.join(cfg.read('data', 'adsb'), "*.TSV")),
            start_date=nvspl_dates[0],
            end_date=nvspl_dates[-1],
            mask=study_area
        )

        raw_tracks["local_hourtime"] = raw_tracks["TIME"].apply(lambda t: t.replace(minute=0, second=0, microsecond=0))
        tracks = Tracks(raw_tracks, id_col='flight_id', datetime_col='TIME', z_col='altitude')
        hourtimes = tracks.local_hourtime.astype(object).unique()

    elif args.track_source == 'Database':
        raw_tracks = query_tracks(engine=engine, start_date=nvspl_dates[0], end_date=nvspl_dates[-1], mask=study_area)
        tracks = Tracks(raw_tracks, 'flight_id', 'ak_datetime', 'altitude_m')
        hourtimes = tracks.ak_hourtime.astype(object).unique()

    else:
        raise NotImplementedError('Code for AIS is not ready yet.')

    track_hours = [{'year': hourtime.year,
                    'month': hourtime.month,
                    'day': hourtime.day,
                    'hour': hourtime.hour}
                   for hourtime in hourtimes]

    # Open NVSPL data files during hours in which there is flight data.
    nvspl_files = [e.path for e in archive.nvspl(unit=args.unit, site=args.site, year=str(args.year), items=track_hours)]
    nvspl = Nvspl(nvspl_files)

    logger.info("Launching application...")
    app.launch(
        tracks=tracks,
        nvspl=nvspl,
        mic=microphone,
        crs=microphone.crs,
        study_area=study_area,
        clip=False
    )
