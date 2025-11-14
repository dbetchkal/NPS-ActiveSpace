from argparse import ArgumentParser
import os
import iyore
import pickle
from nps_active_space.scripts.generate_metrics import get_obs_periods, get_all_srcid_stats
import nps_active_space.utils.config as cfg
from nps_active_space.utils.models import Srcid


def get_acoustic_metrics(unit, site, year, env, track_source):
    """
    Gets acoustic metrics for the period(s) of time with overlapping acoustic and causal data.
    
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
    track_source: str
        'GPS', 'ADSB', or 'AIS'. Metrics will only be calculated for
        time periods with overlapping acoustic and track data.
    """
    assert track_source in ["GPS", "ADSB", "AIS"], "Invalid track source"
    print(unit + site + year)

    cfg.initialize(env)
    project_dir = cfg.read("project", "dir")
    nvspl_archive = cfg.read("data", "nvspl_archive")

    # process track source
    adsb_dir = None
    if track_source == "ADSB":
        adsb_dir = cfg.read("data", "adsb")
    elif track_source == "AIS":
        raise NotImplementedError('Code for AIS is not ready yet.')
    
    obs_periods = get_obs_periods(unit, site, year, nvspl_archive, adsb_dir)
    print(f"Time periods with acoustic and causal data:\n{obs_periods}")

    # load SRCID files (true acoustic record)
    print("Loading Srcid")
    ds = iyore.Dataset(nvspl_archive)
    paths = [e.path for e in ds.srcid(unit=unit, site=site, year=year)]
    assert len(paths) > 0, "No SPLAT annotation data"
    src_obj = Srcid(paths[0])
    src = src_obj.data

    # get metrics
    print("Computing metrics")
    src_stats, src_CIs, src_data = get_all_srcid_stats(src, obs_periods)

    # save output
    savedir = os.path.join(project_dir, unit+site, "Output_Data", "METRICS")
    pkl_file = os.path.join(savedir, f"acoustic_metrics_{unit}{site}{year}.pkl")
    os.makedirs(savedir, exist_ok=True)
    with open(pkl_file, "wb") as f:
        pickle.dump((src_stats, src_CIs, src_data, obs_periods), f)
    print(f"Saved tuple (stats, conf intervals, data, observation periods) to {pkl_file}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-e", "--environment", required=True,
                        help="The configuration environment to run the script in.")
    parser.add_argument("-u", "--unit", help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", help="Four letter site code. E.g. TRLA")
    parser.add_argument("-y", "--year",  help="Four digit year. E.g. 2018")
    parser.add_argument('-t', '--track-source', default='GPS', choices=["GPS", "ADSB", "AIS"],
                          help="Enter 'GPS', 'ADSB', or 'AIS'. Metrics will only be calculated for " \
                          "time periods with overlapping acoustic and track data.")
    args = parser.parse_args()
    
    get_acoustic_metrics(args.unit, args.site, args.year, args.environment, args.track_source)