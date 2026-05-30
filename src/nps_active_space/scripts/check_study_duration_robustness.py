from argparse import ArgumentParser
from nps_active_space.validation.study_duration_stability import fit_varying_n_tracks, plot_stability
import nps_active_space.config as cfg

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("deployments", nargs="*",
                        help="List of deployments to plot, e.g. DENATRLA2023 DENATRLA2024 DENATRLA2025."\
                            "If only one deployment is passed, the plot will be saved in that site's outputs directory instead of the project directory.")
    parser.add_argument("-e", "--environment", required=True,
                        help="The configuration environment to run the script in.")
    parser.add_argument("--col", default="area",
                        help="Which column from the study duration csv to plot. Common values are 'area', 'gain', 'f1'. Default is 'area'.")
    parser.add_argument("-k", type=int, default=10, help="Top k best fits will be plotted in a range around the best fit. Default k=10.")
    parser.add_argument("-n", "--max-n-tracks", type=int,
                        help="Rightmost xlim for plotting; maximum number of tracks to plot.")

    args = parser.parse_args()
    assert len(args.deployments) > 0, "Must pass at least one deployment as a positional argument."

    config = cfg.load_config(args.environment)
    project_dir = config.project.dir

    # for each deployment, fit active space with increasing # of tracks
    # and save results to a CSV
    print("Fitting increasing number of tracks for each deployment...")
    for deployment in args.deployments:
        unit, site, year = deployment[:4], deployment[4:-4], deployment[-4:]
        fit_varying_n_tracks(project_dir, unit, site, year)
        print("")
    
    # make plot
    plot_stability(project_dir, args.deployments, args.col, args.k, args.max_n_tracks)