#!/usr/bin/env python
"""Verify AAM terrain cross-run cache (expect <1s load, not ~10min write_terrain).

Run: docker/run_activespace.sh -m aam docker/validate_aam_terrain_cache.py
"""
import glob
import time

import geopandas as gpd

import nps_active_space.utils.config as cfg
from nps_active_space.propagation_model.aam.model import AamPropagationModel
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.utils.computation import study_area_utm_crs
from nps_active_space.utils.helpers import get_deployment


def main() -> None:
    cfg.initialize("container")
    proj_dir = cfg.read("project", "dir")
    site_dir = f"{proj_dir}/DENATRLA"
    dem, _ = get_project_setup_elevation(site_dir)
    study_area = gpd.read_file(glob.glob(f"{site_dir}/*study*.shp")[0])
    mic = get_deployment(proj_dir, "DENA", "TRLA", 2025, elevation=False)
    mic = mic.to_crs(study_area_utm_crs(study_area.iloc[[0]]))

    model = AamPropagationModel(site_dir)
    start = time.perf_counter()
    model.prepare_site(
        str(dem),
        study_area.iloc[[0]],
        mic,
        project_dem=False,
        suffix=f"_{mic.name}",
    )
    elapsed_s = time.perf_counter() - start
    print(f"[terrain-cache] prepare_site done in {elapsed_s:.1f}s", flush=True)
    if elapsed_s > 30:
        raise SystemExit("cache load slower than expected (>30s); check [aam-terrain] logs")


if __name__ == "__main__":
    main()
