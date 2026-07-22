"""
This script creates the overarching NMSIM project that is used by every module of ``NPS-ActiveSpace``.
This requires two conceptual inputs: (1) a study area polygon, and (2) a listening location point.
We use the study area to prepare a portion of a Digital Elevation Model for use in computation and visualization.
The elevation data is also converted to an antiquated ESRI GridFloat format required by ``NMSIM``.
We represent the listening location coordinates as a canonical NMSIM site (.sit) file.

Library functions live in :mod:`nps_active_space.setup`.
"""

import argparse
from pathlib import Path

import nps_active_space.utils.config as cfg
from nps_active_space.setup import NMSIM_DST_CRS, setup_site

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-e",
        "--environment",
        required=True,
        help="The configuration environment to run the script in.",
    )
    parser.add_argument("-u", "--unit", required=True, help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", required=True, help="Four letter site code. E.g. TRLA")
    parser.add_argument("-y", "--year", type=int, required=True, help="Four digit year. E.g. 2018")
    parser.add_argument(
        "--mic-coord",
        nargs=2,
        type=float,
        required=True,
        help="Mic coordinates x y (WGS84). Example: -136.088360 58.569310",
    )
    parser.add_argument(
        "--studyarea-sw",
        nargs=2,
        type=float,
        required=True,
        help="Study area SW corner x y (WGS84). Example: -136.088360 58.569310",
    )
    parser.add_argument(
        "--studyarea-ne",
        nargs=2,
        type=float,
        required=True,
        help="Study area NE corner x y (WGS84). Example: -135.818994 58.706095",
    )
    args = parser.parse_args()

    cfg.initialize(environment=args.environment)
    project_dir = Path(cfg.read("project", "dir"))
    site_dir = project_dir / f"{args.unit}{args.site}"

    result = setup_site(
        site_dir=site_dir,
        unit=args.unit,
        site=args.site,
        year=args.year,
        mic_coord=tuple(args.mic_coord),
        studyarea_sw=tuple(args.studyarea_sw),
        studyarea_ne=tuple(args.studyarea_ne),
        source_dem=cfg.read("data", "dem"),
        dst_crs=NMSIM_DST_CRS,
        source_elevation_units=cfg.read_dem_elevation_units(),
    )

    print(f"Created a new NMSIM site directory {args.unit}{args.site} in {project_dir}.")
    print(f"Saved study area to {result.study_area_path}")
    print(f"Elevation data have been written to {result.elevation_tif_path}.")
    print("Grid float elevation file components (.flt, .hdr) have been written.")
    print("NMSIM site file has been written. Finished!")

