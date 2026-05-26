import glob
from argparse import ArgumentParser
import geopandas as gpd
import os
import sqlalchemy
import iyore
import nps_active_space.ground_truthing as app
from nps_active_space.utils.models import Nvspl, Tracks
import nps_active_space.utils.config as cfg
from nps_active_space.utils.helpers import get_deployment, get_logger, query_adsb, query_tracks, load_DEM
from nps_active_space.utils.computation import coords_to_utm
from nps_active_space.utils.clock_drift import correct_clock_drift


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
    argparse.add_argument('-t', '--track-source', default='GPS', choices=["GPS", "ADSB", "AIS"],
                          help="Enter 'GPS', 'ADSB', or 'AIS")

    args = argparse.parse_args()

    config = cfg.load_config(args.environment)
    logger = get_logger('GROUND-TRUTHING')
    db_conf = config.database
    engine = sqlalchemy.create_engine(
        f'postgresql://{db_conf.username}:{db_conf.password}'
        f'@{db_conf.host}:{db_conf.port}/{db_conf.name}'
    )

    logger.info(f"Beginning ground truthing process for {args.unit}{args.site}{args.year}...")

    # Set the various path variables.
    archive = iyore.Dataset(config.data.nvspl_archive)
    site_dir = f"{config.project.dir}/{args.unit}{args.site}"
    faa_path = None
    faa_corrections_path = None

    # Load the microphone deployment site metadata and the study area shapefile.
    microphone = get_deployment(config.project.dir, args.unit, args.site, args.year)
    study_area = gpd.read_file(glob.glob(f"{site_dir}/*study*.shp")[0])  # In NAD83, epsg:4269

    # Retrieve the days for which at least some NVSPL data exist.
    nvspl_dates = sorted(set([f"{e.year}-{e.month}-{e.day}" for e in archive.nvspl(unit=args.unit, site=args.site, year=args.year)]))
    assert len(nvspl_dates) > 0, f"No NVSPL data found in archive {config.data.nvspl_archive}"

    # Query flight tracks from days there is NVSPL data for.
    logger.info("Querying tracks...")

    if args.track_source == 'ADSB':
        raw_tracks = query_adsb(
            adsb_path=config.data.adsb,
            start_date=nvspl_dates[0],
            end_date=nvspl_dates[-1],
            mask=study_area
        )
        tracks = Tracks(raw_tracks, id_col='flight_id', datetime_col='TIME', z_col='altitude')
        faa_path = config.project.FAA_Releasable_db
        faa_corrections_path = config.project.FAA_type_corrections

    elif args.track_source == 'GPS':
        raw_tracks = query_tracks(engine=engine, start_date=nvspl_dates[0], end_date=nvspl_dates[-1], mask=study_area)
        tracks = Tracks(raw_tracks, 'flight_id', 'ak_datetime', 'altitude_m')
        faa_path = config.project.FAA_Releasable_db
        faa_corrections_path = config.project.FAA_type_corrections

    else:
        raise NotImplementedError('Code for AIS is not ready yet.')

    assert not tracks.empty, "No tracks loaded, is your track source correct?"
    
    # correct for clock drift
    clock_drift_file = os.path.join(site_dir, f"{args.unit}{args.site}{args.year}_clock_drift_{args.track_source}.csv")
    if os.path.exists(clock_drift_file):
        print(f"Found clock drift correction file, using it: {os.path.basename(clock_drift_file)}")
        correct_clock_drift(tracks, clock_drift_file, inplace=True)
        
    # Open NVSPL data files during hours in which there is flight data.
    hourtimes = tracks["point_dt"].dt.floor("h").unique()
    track_hours = [{'year': hourtime.year,
                    'month': hourtime.month,
                    'day': hourtime.day,
                    'hour': hourtime.hour}
                   for hourtime in hourtimes]
    nvspl_files = [e.path for e in archive.nvspl(unit=args.unit, site=args.site, year=str(args.year), items=track_hours)]
    nvspl = Nvspl(nvspl_files)

    # Load DEM
    dem = load_DEM(config.project.dir, args.unit, args.site)

    logger.info("Launching application...")
    app.launch(
        tracks=tracks,
        nvspl=nvspl,
        mic=microphone,
        crs=coords_to_utm(microphone.lat, microphone.lon),
        study_area=study_area,
        database_type=args.track_source,
        dem=dem,
        clip=False,
        faa_path=faa_path,
        faa_corrections_path=faa_corrections_path
    )
