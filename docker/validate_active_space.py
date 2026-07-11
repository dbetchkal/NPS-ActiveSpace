#!/usr/bin/env python
"""
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
import threading
import time
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

# Before any package import that might pull matplotlib.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import geopandas as gpd
import iyore
import pandas as pd
from tqdm import tqdm

import nps_active_space.utils.config as cfg
from nps_active_space.utils.computation import ambience_from_nvspl
from nps_active_space.utils.helpers import get_deployment, get_omni_sources
from nps_active_space.utils.models import Microphone, Nvspl

from nps_active_space.active_space import ActiveSpaceGenerator

HEARTBEAT_INTERVAL_S = 15.0


@dataclass(frozen=True)
class ValidateArgs:
    environment: str
    unit: str
    site: str
    year: int
    gains: list[float]
    altitude: int
    heading: int
    density: int
    cpus: int


@dataclass(frozen=True)
class ValidationResult:
    gain: float
    n_tested: int
    n_audible: int
    area_km2: float
    outfile: str
    elapsed_s: float


@dataclass(frozen=True)
class RunContext:
    gen: ActiveSpaceGenerator
    mic: Microphone
    heading: int
    altitude: int
    density: int
    out_dir: str
    unit: str
    site: str
    year: int


def log(msg: str) -> None:
    print(f"[validate] {msg}", flush=True)


@contextmanager
def timed_step(label: str) -> Iterator[None]:
    log(f"{label}...")
    start = time.perf_counter()
    try:
        yield
    finally:
        log(f"{label} done ({time.perf_counter() - start:.1f}s)")


@contextmanager
def heartbeat(label: str, interval_s: float = HEARTBEAT_INTERVAL_S) -> Iterator[None]:
    """Emit periodic messages while a long-running step (NMSim via Wine) is in progress."""
    stop = threading.Event()
    start = time.perf_counter()

    def _pulse() -> None:
        while not stop.wait(interval_s):
            elapsed = time.perf_counter() - start
            log(f"{label} ... still running ({elapsed:.0f}s elapsed)")

    thread = threading.Thread(target=_pulse, daemon=True)
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join(timeout=0.1)


def parse_args() -> ValidateArgs:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-e", "--environment", default="container")
    ap.add_argument("-u", "--unit", default="DENA")
    ap.add_argument("-s", "--site", default="TRLA")
    ap.add_argument("-y", "--year", type=int, default=2025)
    ap.add_argument(
        "--gains",
        type=float,
        nargs="+",
        default=[0.0],
        help="One or more omni source gains (dB). >1 gain runs concurrently via mp.Pool "
        "(mirrors generate_active_space.py) to validate Wine under concurrency.",
    )
    ap.add_argument("--altitude", type=int, default=1000, help="Source altitude in meters.")
    ap.add_argument("--heading", type=int, default=0)
    ap.add_argument(
        "--density",
        type=int,
        default=12,
        help="src_pt_density (NxN mesh). Keep small for a quick check; 48 is the pipeline default.",
    )
    ap.add_argument("--cpus", type=int, default=0, help="Worker processes for multi-gain runs (0 = len(gains)).")
    ns = ap.parse_args()
    return ValidateArgs(
        environment=ns.environment,
        unit=ns.unit,
        site=ns.site,
        year=ns.year,
        gains=ns.gains,
        altitude=ns.altitude,
        heading=ns.heading,
        density=ns.density,
        cpus=ns.cpus,
    )


def initialize_site(args: ValidateArgs) -> tuple[str, str, Microphone, gpd.GeoDataFrame]:
    cfg.initialize(environment=args.environment)
    proj_dir = cfg.read("project", "dir")
    site_dir = f"{proj_dir}/{args.unit}{args.site}"

    nmsim = cfg.read("project", "nmsim")
    dem = cfg.read("data", "dem")
    log(f"env={args.environment} site={args.unit}{args.site}{args.year}")
    log(f"nmsim={nmsim}")
    log(f"dem={dem}")

    with timed_step("loading microphone deployment"):
        mic = get_deployment(proj_dir, args.unit, args.site, args.year, elevation=False)
    log(f"mic={mic.name} lat={mic.lat:.5f} lon={mic.lon:.5f} z={mic.z}")

    with timed_step("loading study area"):
        study_area = gpd.read_file(glob.glob(f"{site_dir}/*study*.shp")[0])

    return site_dir, dem, mic, study_area


def load_ambience(unit: str, site: str, year: int) -> pd.Series:
    archive = iyore.Dataset(cfg.read("data", "nvspl_archive"))
    nvspl_files = [e.path for e in archive.nvspl(unit=unit, site=site, year=str(year))]
    if not nvspl_files:
        raise RuntimeError("No NVSPL files found for ambience")
    log(f"{len(nvspl_files)} NVSPL file(s) -> ambience (L90)")
    with timed_step("computing ambience from NVSPL"):
        nvspl = Nvspl(nvspl_files)
        return ambience_from_nvspl(nvspl, quantile=90, broadband=False)


def build_generator(
    site_dir: str,
    dem: str,
    study_area: gpd.GeoDataFrame,
    ambience: pd.Series,
    mic: Microphone,
) -> ActiveSpaceGenerator:
    nmsim = cfg.read("project", "nmsim")
    gen = ActiveSpaceGenerator(
        NMSIM=nmsim,
        root_dir=site_dir,
        study_area=study_area,
        ambience=ambience,
        dem_src=dem,
    )
    with timed_step("masking/projecting DEM (set_dem)"):
        gen.set_dem(mic)
    return gen


def _output_path(ctx: RunContext, gain: float) -> str:
    return os.path.join(
        ctx.out_dir,
        f"{ctx.unit}{ctx.site}{ctx.year}_{ctx.altitude}m_gain{gain}_VALIDATE.geojson",
    )


def _area_km2(active_space: gpd.GeoDataFrame) -> float:
    return round(active_space.to_crs(active_space.estimate_utm_crs()).area.sum() / 1e6, 2)


def run_one_gain(gain: float, ctx: RunContext) -> ValidationResult:
    """Run a single gain end-to-end (NMSim via wine) and write the active-space geojson.

    Used both directly and as an mp.Pool worker (Linux fork inherits ``ctx.gen``).
    """
    omni = get_omni_sources(lower=gain, upper=gain)[0]
    mic = deepcopy(ctx.mic)
    mic.name = f"{mic.name}{Path(omni).stem}"

    label = (
        f"gain={gain} alt={ctx.altitude}m heading={ctx.heading} "
        f"density={ctx.density} (NMSim via Wine)"
    )
    start = time.perf_counter()
    with heartbeat(label):
        active_space, tested_pts = ctx.gen.generate(
            omni_source=omni,
            mic=mic,
            heading=ctx.heading,
            altitude_m=ctx.altitude,
            src_pt_density=ctx.density,
        )

    n_audible = int((tested_pts["audible"] == 1).sum()) if "audible" in tested_pts else -1
    out = _output_path(ctx, gain)
    active_space.to_file(out, driver="GeoJSON")
    elapsed_s = time.perf_counter() - start
    return ValidationResult(
        gain=gain,
        n_tested=len(tested_pts),
        n_audible=n_audible,
        area_km2=_area_km2(active_space),
        outfile=os.path.basename(out),
        elapsed_s=round(elapsed_s, 1),
    )


def _run_one_worker(gain: float, ctx: RunContext) -> ValidationResult:
    """Pool entry point: unpack is avoided so starmap can pass a single context object."""
    return run_one_gain(gain, ctx)


def run_all_gains(gains: Sequence[float], ctx: RunContext, ncpu: int) -> list[ValidationResult]:
    if len(gains) == 1:
        gain = gains[0]
        log(
            f"running NMSim via wine: gain={gain} alt={ctx.altitude}m "
            f"heading={ctx.heading} density={ctx.density}"
        )
        return [run_one_gain(gain, ctx)]

    log(
        f"CONCURRENCY: {len(gains)} gains {list(gains)} across {ncpu} workers "
        f"(shared Wine prefix/server) -> testing wine under multiprocessing"
    )
    results: list[ValidationResult] = []
    with mp.Pool(ncpu) as pool:
        with tqdm(desc="Gains", unit="gain", colour="green", total=len(gains)) as pbar:
            processes = [
                pool.apply_async(_run_one_worker, args=(gain, ctx), callback=lambda _: pbar.update())
                for gain in gains
            ]
            for proc in processes:
                results.append(proc.get())
    return results


def print_results(results: Sequence[ValidationResult], total_elapsed_s: float) -> None:
    log("=== results ===")
    for result in sorted(results, key=lambda r: r.gain):
        log(
            f"  gain={result.gain:>5}  tested={result.n_tested:>4}  "
            f"audible={result.n_audible:>4}  area~={result.area_km2} km^2  "
            f"elapsed={result.elapsed_s:>5}s  -> {result.outfile}"
        )
    log(f"OK: {len(results)} active space(s) generated in {total_elapsed_s:.1f}s wall time")


def main() -> None:
    args = parse_args()
    wall_start = time.perf_counter()

    site_dir, dem, mic, study_area = initialize_site(args)
    ambience = load_ambience(args.unit, args.site, args.year)
    gen = build_generator(site_dir, dem, study_area, ambience, mic)

    out_dir = os.path.join(site_dir, "Output_Data", "ACTIVESPACES")
    os.makedirs(out_dir, exist_ok=True)
    ctx = RunContext(
        gen=gen,
        mic=mic,
        heading=args.heading,
        altitude=args.altitude,
        density=args.density,
        out_dir=out_dir,
        unit=args.unit,
        site=args.site,
        year=args.year,
    )

    ncpu = args.cpus or len(args.gains)
    results = run_all_gains(args.gains, ctx, ncpu)
    print_results(results, time.perf_counter() - wall_start)


if __name__ == "__main__":
    main()
