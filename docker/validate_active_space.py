#!/usr/bin/env python
"""
TODO split into some functions and add type hints
Integration check for running NMSim through Wine inside the Linux/Mac container.

This script uses ActiveSpaceGenerator.generate() directly for a single site, gain, altitude
and heading, so it exercises the real NMSim-via-wine path (control/batch file creation,
the `project.nmsim` shim, .tis parsing, audibility) WITHOUT requiring ground-truth
annotations (those are only used by generate_active_space.py for gain auto-selection and
precision/recall scoring, not to run the basic active space computation step).

Run inside the container, e.g.:
  docker/run_activespace.sh docker/validate_active_space.py -u DENA -s TRLA -y 2025 \
      --gain 0 --altitude 1000 --density 12 --heading 0
"""
import argparse
import glob
import multiprocessing as mp
import os
from copy import deepcopy
from pathlib import Path

import geopandas as gpd
import iyore

import nps_active_space.utils.config as cfg
from nps_active_space.utils.helpers import get_deployment, get_omni_sources
from nps_active_space.utils.models import Nvspl
from nps_active_space.utils.computation import ambience_from_nvspl
from nps_active_space.active_space import ActiveSpaceGenerator


def _run_one(gain, gen, mic, heading, altitude, density, out_dir, unit, site, year):
    """Run a single gain end-to-end (NMSim via wine) and write the active-space geojson.
    Used both directly and as an mp.Pool worker (Linux fork inherits `gen`)."""
    omni = get_omni_sources(lower=gain, upper=gain)[0]
    # Unique mic name per source so concurrent workers don't collide on trajectory/job
    # filenames (mirrors generate_active_space.py._run_active_space).
    mic = deepcopy(mic)
    mic.name = f"{mic.name}{Path(omni).stem}"
    active_space, tested_pts = gen.generate(
        omni_source=omni, mic=mic, heading=heading, altitude_m=altitude, src_pt_density=density,
    )
    n_aud = int((tested_pts['audible'] == 1).sum()) if 'audible' in tested_pts else -1
    out = os.path.join(out_dir, f"{unit}{site}{year}_{altitude}m_gain{gain}_VALIDATE.geojson")
    active_space.to_file(out, driver='GeoJSON')
    area_km2 = active_space.to_crs(active_space.estimate_utm_crs()).area.sum() / 1e6
    return gain, len(tested_pts), n_aud, round(area_km2, 2), os.path.basename(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-e', '--environment', default='container')
    ap.add_argument('-u', '--unit', default='DENA')
    ap.add_argument('-s', '--site', default='TRLA')
    ap.add_argument('-y', '--year', type=int, default=2025)
    ap.add_argument('--gains', type=float, nargs='+', default=[0.0],
                    help='One or more omni source gains (dB). >1 gain runs concurrently via mp.Pool '
                         '(mirrors generate_active_space.py) to validate Wine under concurrency.')
    ap.add_argument('--altitude', type=int, default=1000, help='Source altitude in meters.')
    ap.add_argument('--heading', type=int, default=0)
    ap.add_argument('--density', type=int, default=12,
                    help='src_pt_density (NxN mesh). Keep small for a quick check; 48 is the pipeline default.')
    ap.add_argument('--cpus', type=int, default=0, help='Worker processes for multi-gain runs (0 = len(gains)).')
    args = ap.parse_args()

    cfg.initialize(environment=args.environment)
    proj_dir = cfg.read('project', 'dir')
    site_dir = f"{proj_dir}/{args.unit}{args.site}"

    nmsim = cfg.read('project', 'nmsim')
    dem = cfg.read('data', 'dem')
    print(f"[validate] env={args.environment} site={args.unit}{args.site}{args.year}")
    print(f"[validate] nmsim={nmsim}")
    print(f"[validate] dem={dem}")

    mic = get_deployment(proj_dir, args.unit, args.site, args.year, elevation=False)
    print(f"[validate] mic={mic.name} lat={mic.lat:.5f} lon={mic.lon:.5f} z={mic.z}")

    study_area = gpd.read_file(glob.glob(f"{site_dir}/*study*.shp")[0])

    archive = iyore.Dataset(cfg.read('data', 'nvspl_archive'))
    nvspl_files = [e.path for e in archive.nvspl(unit=args.unit, site=args.site, year=str(args.year))]
    assert nvspl_files, "No NVSPL files found for ambience"
    print(f"[validate] {len(nvspl_files)} NVSPL file(s) -> ambience (L90)")
    nvspl = Nvspl(nvspl_files)
    ambience = ambience_from_nvspl(nvspl, quantile=90, broadband=False)

    gen = ActiveSpaceGenerator(
        NMSIM=nmsim,
        root_dir=site_dir,
        study_area=study_area,
        ambience=ambience,
        dem_src=dem,
    )
    print("[validate] masking/projecting DEM (set_dem)...")
    gen.set_dem(mic)  # done ONCE, before any concurrency (mirrors generate_active_space.py)

    out_dir = os.path.join(site_dir, "Output_Data", "ACTIVESPACES")
    os.makedirs(out_dir, exist_ok=True)
    common = dict(gen=gen, mic=mic, heading=args.heading, altitude=args.altitude,
                  density=args.density, out_dir=out_dir, unit=args.unit, site=args.site, year=args.year)

    if len(args.gains) == 1:
        print(f"[validate] running NMSim via wine: gain={args.gains[0]} "
              f"alt={args.altitude}m heading={args.heading} density={args.density}")
        results = [_run_one(args.gains[0], **common)]
    else:
        ncpu = args.cpus or len(args.gains)
        print(f"[validate] CONCURRENCY: {len(args.gains)} gains {args.gains} across {ncpu} workers "
              f"(shared Wine prefix/server) -> testing wine under multiprocessing")
        with mp.Pool(ncpu) as pool:
            results = pool.starmap(_run_one, [(g,) + tuple(common.values()) for g in args.gains])

    print("[validate] === results ===")
    for gain, n_tested, n_aud, area, name in sorted(results):
        print(f"[validate]   gain={gain:>5}  tested={n_tested:>4}  audible={n_aud:>4}  area~={area} km^2  -> {name}")
    print(f"[validate] OK: {len(results)}/{len(args.gains)} active spaces generated")


if __name__ == '__main__':
    main()
