from argparse import ArgumentParser
import os
import iyore
import nps_active_space.ground_truthing as app
from nps_active_space.ground_truthing.load_tracks import load_tracks
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.models import Nvspl
import nps_active_space.utils.config as cfg
from nps_active_space.utils.helpers import (
    get_deployment,
    get_logger,
    load_DEM,
    load_studyarea,
)
from nps_active_space.utils import paths as p
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
    argparse.add_argument(
        '-t', '--track-source',
        default=TrackSource.GPS,
        type=TrackSource,
        choices=list(TrackSource),
        help="Enter 'GPS', 'ADSB', or 'AIS'",
    )

    args = argparse.parse_args()

    cfg.initialize(environment=args.environment)
    logger = get_logger('GROUND-TRUTHING')

    logger.info(f"Beginning ground truthing process for {args.unit}{args.site}{args.year}...")

    # Set the various path variables.
    archive = iyore.Dataset(cfg.read('data', 'nvspl_archive'))
    site_dir = p.site_dir(cfg.read('project', 'dir'), args.unit, args.site)

    # Load the microphone deployment site metadata and the study area shapefile.
    project_dir = cfg.read("project", "dir")
    microphone = get_deployment(project_dir, args.unit, args.site, args.year)
    study_area = load_studyarea(project_dir, args.unit, args.site, args.year)

    # Retrieve the days for which at least some NVSPL data exist.
    nvspl_dates = sorted(set([f"{e.year}-{e.month}-{e.day}" for e in archive.nvspl(unit=args.unit, site=args.site, year=args.year)]))
    assert len(nvspl_dates) > 0, f"No NVSPL data found in archive {cfg.read('data', 'nvspl_archive')}"

    # Query flight tracks from days there is NVSPL data for.
    logger.info("Querying tracks...")

    loaded = load_tracks(
        args.track_source,
        start_date=nvspl_dates[0],
        end_date=nvspl_dates[-1],
        study_area=study_area,
        microphone=microphone,
    )
    tracks = loaded.tracks
    faa_path = loaded.faa_path
    faa_corrections_path = loaded.faa_corrections_path

    assert not tracks.empty, "No tracks loaded, is your track source correct?"
    logger.info(
        f"Loaded {tracks['track_id'].nunique()} tracks "
        f"({len(tracks)} points) from {nvspl_dates[0]} to {nvspl_dates[-1]}"
    )

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
    dem = load_DEM(cfg.read('project', 'dir'), args.unit, args.site)

    logger.info("Launching application...")
    app.launch(
        tracks=tracks,
        nvspl=nvspl,
        mic=microphone,
        crs=microphone.crs,
        study_area=study_area,
        track_source=args.track_source,
        dem=dem,
        clip=False,
        faa_path=faa_path,
        faa_corrections_path=faa_corrections_path
    )
