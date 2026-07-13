import os
import pandas as pd
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import pickle
from nps_active_space.scripts.run_audible_transits import init_audible_transits, AudibleTransits
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.metrics import get_obs_periods, get_all_geo_stats
import nps_active_space.utils.config as cfg


def get_optimal_3d_gain(project_dir, unit, site, year):
    """Get optimal 3D gain from the fits.csv file located in the project directory."""
    csv_3d_fits = os.path.join(project_dir, "fits.csv")
    if not os.path.exists(csv_3d_fits):
        raise FileNotFoundError(f"Couldn't find {csv_3d_fits} for getting the optimal gain.")
    fit_results = pd.read_csv(csv_3d_fits, index_col="Designator")
    if unit+site+year in fit_results.index:
        gain_3d = float(fit_results.loc[unit+site+year, "1/3rd Octave Gain (F1)"])
    else:
        raise RuntimeError(f"No fitted gain for {unit}{site}{year}")
    
    return gain_3d


def get_geographic_metrics(unit, site, year, env, track_source: TrackSource, transits_pkl=None):
    """
    Gets geographic metrics for the period(s) of time with overlapping acoustic and causal data.
    
    Parameters
    ----------
    unit: str
        Four letter unit code. E.g. DENA
    site: str
        Four letter site code. E.g. TRLA
    year: str
        Four digit year. E.g. "2018"
    env: str
        The configuration environment to run the script in. E.g. "DENA_streamline"
    track_source: TrackSource
        GPS, ADSB, or AIS. Metrics will only be calculated for
        time periods with overlapping acoustic and track data.
    transits_pkl: str (Optional)
        If provided, path to a .pkl file to load audible transits from instead of the default.

    Returns
    -------
    Tuple (stats, CIs, data, obs_periods)
        What is stored in the output .pkl file. See get_all_geo_stats() and get_obs_periods() for documentation.
    """
    year = str(year)
    print(unit + site + year)

    cfg.initialize(env)
    project_dir = cfg.read("project", "dir")
    nvspl_archive = cfg.read("data", "nvspl_archive")

    # process track source
    adsb_dir = None
    match track_source:
        case TrackSource.ADSB:
            adsb_dir = cfg.read("data", "adsb")
        case TrackSource.AIS:
            raise NotImplementedError('Code for AIS is not ready yet.')
    
    obs_periods = get_obs_periods(unit, site, year, nvspl_archive, adsb_dir)
    print(f"Time periods with acoustic and causal data:\n{obs_periods}")

    
    print("Loading audible transits")

    # determine audible transits pkl file path if not provided
    if transits_pkl is None:

        # determine optimal gain from CSV file
        gain_3d = get_optimal_3d_gain(project_dir, unit, site, year)
    
        # make a dummy listener to determine default save path
        metadata = {"unit": unit, "site": site, "activespace year": year, "gain": gain_3d,
                    "study start": f"{year}-01-01", "study end": f"{year}-12-31",
                    "database type": track_source, "env": env}
        listener = init_audible_transits(metadata, paths={})
        transits_pkl = os.path.join(listener.default_output_dir(), listener.pkl_filename)
    
    if not os.path.exists(transits_pkl):
        raise FileNotFoundError(f"Couldn't find audible transits output at {transits_pkl}")
    
    # load audible transits
    listener = AudibleTransits.from_pickle(transits_pkl)
    study_year = str(pd.Timestamp(listener.study_start).year)
    assert study_year == year, f"Audible transits study year ({study_year}) doesn't match 'year' argument ({year})."
    tracks = listener.tracks
    tracks = tracks[tracks["aircraft_type"] == "Fixed-wing"]
    # close excess figures
    plt.close(listener.overflights_fig)
    plt.close(listener.transits_fig)
    
    # get metrics
    print("Computing metrics")
    geo_stats, geo_CIs, geo_data = get_all_geo_stats(tracks, obs_periods)

    # save output
    savedir = os.path.join(project_dir, unit+site, "Output_Data", "METRICS")
    os.makedirs(savedir, exist_ok=True)
    # include active space year in filename in case we are experimenting with
    # other years' active spaces to predict metrics for this year
    pkl_file = os.path.join(savedir, f"geographic_metrics_{unit}{site}{year}" \
                            f"_active_space_{listener.activespace_year}.pkl")
    with open(pkl_file, "wb") as f:
        pickle.dump((geo_stats, geo_CIs, geo_data, obs_periods), f)
    print(f"Saved tuple (stats, conf intervals, data, observation periods) to {pkl_file}")

    return geo_stats, geo_CIs, geo_data, obs_periods


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-e", "--environment", required=True,
                        help="The configuration environment to run the script in.")
    parser.add_argument("-u", "--unit", help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", help="Four letter site code. E.g. TRLA")
    parser.add_argument("-y", "--year",  help="Four digit year. E.g. 2018")
    parser.add_argument(
        '-t', '--track-source',
        default=TrackSource.GPS,
        type=TrackSource,
        choices=list(TrackSource),
        help="Enter 'GPS', 'ADSB', or 'AIS'. Metrics will only be calculated for "
        "time periods with overlapping acoustic and track data.",
    )
    parser.add_argument("--transits-pkl", help="Path to .pkl file to load audible transits from, if not the default.")

    args = parser.parse_args()
    
    get_geographic_metrics(args.unit, args.site, args.year, args.environment,
                           args.track_source, args.transits_pkl)