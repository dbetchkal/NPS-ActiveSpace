#!/usr/bin/env python
"""
Integration check for ActiveSpaceGenerator with NMSim or AAM propagation (Docker+Wine).

Exercises the real propagation path (control/batch or AAM .inp, shim, audibility, polygon)
for a single site, altitude, and heading. Optional ``--fit`` loads ground-truthing
annotations (audible/inaudible track segments), generates active spaces across a gain
range, and selects the optimal gain via precision/recall (same logic as
``generate_active_space.py``).

Run inside the container, e.g.:
  docker/run_activespace.sh docker/validate_active_space.py -u DENA -s TRLA -y 2025 \
      --gains 0 --altitude 1000 --density 10 --heading 0

  docker/run_activespace.sh -m aam docker/validate_active_space.py --model aam \
      -u DENA -s TRLA -y 2025 --fit --omni-min 0 --omni-max 2 --altitude 1000 --density 10
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
from shapely.geometry import Point
from tqdm import tqdm

import nps_active_space.utils.config as cfg
import nps_active_space.utils.paths as p
from nps_active_space.utils.computation import (
    ambience_from_nvspl,
    normalize_point_density,
    select_optimal,
)
from nps_active_space.utils.helpers import (
    get_deployment,
    get_omni_sources,
    load_annotations,
    omni_to_gain,
)
from nps_active_space.utils.models import Microphone, Nvspl
from nps_active_space.utils.enums import AcousticModel

from nps_active_space.active_space import ActiveSpaceGenerator
from nps_active_space.active_space.active_space_setup import (
    build_active_space_generator,
    precision_recall_plot_path,
    resolve_acoustic_model,
)

HEARTBEAT_INTERVAL_S = 15.0
AAM_SHIM = Path("/usr/local/bin/aam")


@dataclass(frozen=True)
class ValidateArgs:
    environment: str
    unit: str
    site: str
    year: int
    model: AcousticModel
    gains: list[float]
    omni_min: float
    omni_max: float
    omni_step: float
    fit: bool
    beta: float
    altitude: int
    heading: int
    density: int
    cpus: int


@dataclass(frozen=True)
class ValidationResult:
    gain: float
    omni_stem: str
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
    model: AcousticModel


@dataclass(frozen=True)
class FitResult:
    best_omni: str
    max_fbeta: float
    best_precision: float
    best_recall: float
    plot_path: str
    csv_path: str


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
    """Emit periodic messages while a long-running propagation step is in progress."""
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
        "--model",
        type=AcousticModel,
        choices=list(AcousticModel),
        default=None,
        help="Propagation model (default: nmsim, or aam when ACOUSTIC_MODEL=aam).",
    )
    ap.add_argument(
        "--gains",
        type=float,
        nargs="+",
        default=None,
        help="Omni source gains (dB). Ignored when --fit is set (uses --omni-min/--omni-max).",
    )
    ap.add_argument(
        "--fit",
        action="store_true",
        help="Fit optimal gain using ground-truthing annotations (saved_annotations geojson).",
    )
    ap.add_argument(
        "--omni-min",
        type=float,
        default=0.0,
        help="Minimum gain (dB) for --fit gain sweep (matches generate_active_space.py default).",
    )
    ap.add_argument(
        "--omni-max",
        type=float,
        default=2.0,
        help="Maximum gain (dB) for --fit gain sweep.",
    )
    ap.add_argument(
        "--omni-step",
        type=float,
        default=0.5,
        help="Spacing between omni gains in dB for --fit (multiple of 0.5).",
    )
    ap.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="F-beta for --fit (1.0 = F1).",
    )
    ap.add_argument("--altitude", type=int, default=1000, help="Source altitude in meters.")
    ap.add_argument("--heading", type=int, default=0)
    ap.add_argument(
        "--density",
        type=int,
        default=12,
        help="src_pt_density (NxN mesh). Keep small for a quick check; 48 is the pipeline default.",
    )
    ap.add_argument(
        "--cpus",
        type=int,
        default=0,
        help="Worker processes for multi-gain runs (0 = len(gains)).",
    )
    ns = ap.parse_args()

    model = resolve_acoustic_model(ns.model)

    if ns.fit:
        gains = [
            omni_to_gain(src)
            for src in get_omni_sources(ns.omni_min, ns.omni_max, ns.omni_step)
        ]
    elif ns.gains is not None:
        gains = ns.gains
    else:
        gains = [0.0]

    return ValidateArgs(
        environment=ns.environment,
        unit=ns.unit,
        site=ns.site,
        year=ns.year,
        model=model,
        gains=gains,
        omni_min=ns.omni_min,
        omni_max=ns.omni_max,
        omni_step=ns.omni_step,
        fit=ns.fit,
        beta=ns.beta,
        altitude=ns.altitude,
        heading=ns.heading,
        density=ns.density,
        cpus=ns.cpus,
    )


def initialize_site(args: ValidateArgs) -> tuple[str, Microphone, gpd.GeoDataFrame]:
    cfg.initialize(environment=args.environment)
    proj_dir = cfg.read("project", "dir")
    site_dir = f"{proj_dir}/{args.unit}{args.site}"

    log(f"env={args.environment} site={args.unit}{args.site}{args.year} model={args.model}")
    if args.model is AcousticModel.NMSIM:
        log(f"nmsim={cfg.read('project', 'nmsim')}")
    log("elevation from project_setup artifacts under Input_Data/01_ELEVATION")

    with timed_step("loading microphone deployment"):
        mic = get_deployment(proj_dir, args.unit, args.site, args.year, elevation=False)
    log(f"mic={mic.name} lat={mic.lat:.5f} lon={mic.lon:.5f} z={mic.z}")

    with timed_step("loading study area"):
        study_area = gpd.read_file(glob.glob(f"{site_dir}/*study*.shp")[0])

    return site_dir, mic, study_area


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
    study_area: gpd.GeoDataFrame,
    ambience: pd.Series,
    mic: Microphone,
    model: AcousticModel,
) -> ActiveSpaceGenerator:
    with timed_step("caching project_setup elevation (set_dem)"):
        return build_active_space_generator(
            site_dir, study_area, ambience, mic, model, aam_shim=str(AAM_SHIM),
        )


def output_dir(site_dir: str, args: ValidateArgs) -> str:
    return p.activespace_layer_dir(
        site_dir, args.unit, args.site, args.year, args.altitude, args.model,
    )


def _output_path(ctx: RunContext, gain: float, omni_stem: str) -> str:
    usy = p.deployment_id(ctx.unit, ctx.site, ctx.year)
    return os.path.join(ctx.out_dir, f"{usy}_{omni_stem}.geojson")


def _area_km2(active_space: gpd.GeoDataFrame) -> float:
    return round(active_space.to_crs(active_space.estimate_utm_crs()).area.sum() / 1e6, 2)


def _model_label(model: AcousticModel) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            return "AAM"
        case AcousticModel.NMSIM:
            return "NMSim via Wine"


def run_one_gain(gain: float, ctx: RunContext) -> ValidationResult:
    """Run a single gain end-to-end and write the active-space geojson."""
    omni_sources = get_omni_sources(gain, gain)
    omni = omni_sources[0]
    omni_stem = Path(omni).stem
    mic = deepcopy(ctx.mic)
    mic.name = f"{mic.name}{omni_stem}"

    label = (
        f"gain={gain} alt={ctx.altitude}m heading={ctx.heading} "
        f"density={ctx.density} ({_model_label(ctx.model)})"
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
    out = _output_path(ctx, gain, omni_stem)
    active_space.to_file(out, driver="GeoJSON")
    elapsed_s = time.perf_counter() - start
    return ValidationResult(
        gain=gain,
        omni_stem=omni_stem,
        n_tested=len(tested_pts),
        n_audible=n_audible,
        area_km2=_area_km2(active_space),
        outfile=os.path.basename(out),
        elapsed_s=round(elapsed_s, 1),
    )


def _run_one_worker(gain: float, ctx: RunContext) -> ValidationResult:
    return run_one_gain(gain, ctx)


def run_all_gains(gains: Sequence[float], ctx: RunContext, ncpu: int) -> list[ValidationResult]:
    if len(gains) == 1:
        gain = gains[0]
        log(
            f"running {_model_label(ctx.model)}: gain={gain} alt={ctx.altitude}m "
            f"heading={ctx.heading} density={ctx.density}"
        )
        return [run_one_gain(gain, ctx)]

    log(
        f"CONCURRENCY: {len(gains)} gains {list(gains)} across {ncpu} workers "
        f"({_model_label(ctx.model)})"
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


def load_valid_points(
    site_dir: str,
    study_area: gpd.GeoDataFrame,
    args: ValidateArgs,
) -> gpd.GeoDataFrame:
    proj_dir = cfg.read("project", "dir")
    with timed_step("loading ground-truthing annotations"):
        annotations = load_annotations(proj_dir, args.unit, args.site, args.year)
    if annotations.empty:
        raise RuntimeError(
            f"No saved_annotations geojson found for {args.unit}{args.site}{args.year}",
        )
    log(f"{len(annotations)} valid annotated segments")

    valid_points_lst = []
    for idx, row in tqdm(
        annotations.iterrows(),
        total=len(annotations),
        desc="Extracting valid points",
        unit="segment",
    ):
        valid_points_lst.extend([
            {
                "annotation_idx": idx,
                "audible": row.audible,
                "geometry": Point(coords),
            }
            for coords in row.geometry.coords
        ])
    valid_points = gpd.GeoDataFrame(
        data=valid_points_lst,
        geometry="geometry",
        crs=annotations.crs,
    )
    valid_points = normalize_point_density(valid_points, study_area, random_seed=679)
    log(f"{len(valid_points)} annotation points after density normalization")
    return valid_points


def run_fit(
    results: Sequence[ValidationResult],
    valid_points: gpd.GeoDataFrame,
    site_dir: str,
    args: ValidateArgs,
) -> FitResult:
    active_space_polygons: list[tuple[str, gpd.GeoDataFrame]] = []
    for result in results:
        path = os.path.join(output_dir(site_dir, args), result.outfile)
        active_space_polygons.append((result.omni_stem, gpd.read_file(path)))

    plot_path = precision_recall_plot_path(
        site_dir, args.unit, args.site, args.year, args.altitude, args.beta, args.model,
    )
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)

    log(
        f"fitting F-{args.beta} across {len(active_space_polygons)} active spaces "
        f"vs {len(valid_points)} annotation points"
    )
    best_omni, max_fbeta, best_precision, best_recall, _ = select_optimal(
        unit=args.unit,
        site=args.site,
        year=args.year,
        valid_points=valid_points,
        active_space_polygons=active_space_polygons,
        beta_=args.beta,
        plot=True,
        plot_savepath=plot_path,
        verbose=True,
    )

    project_dir = cfg.read("project", "dir")
    csv_path = p.fits_csv(project_dir)

    return FitResult(
        best_omni=best_omni,
        max_fbeta=max_fbeta,
        best_precision=best_precision,
        best_recall=best_recall,
        plot_path=plot_path,
        csv_path=csv_path,
    )


def print_results(results: Sequence[ValidationResult], total_elapsed_s: float) -> None:
    log("=== results ===")
    for result in sorted(results, key=lambda r: r.gain):
        log(
            f"  gain={result.gain:>5}  omni={result.omni_stem}  tested={result.n_tested:>4}  "
            f"audible={result.n_audible:>4}  area~={result.area_km2} km^2  "
            f"elapsed={result.elapsed_s:>5}s  -> {result.outfile}"
        )
    log(f"OK: {len(results)} active space(s) generated in {total_elapsed_s:.1f}s wall time")


def append_aam_site_log(
    site_dir: str,
    results: Sequence[ValidationResult],
    fit: FitResult | None = None,
) -> None:
    from nps_active_space.propagation_model.aam.run_log import append_aam_run_summary

    summary_lines = [
        (
            f"gain={result.gain} omni={result.omni_stem} tested={result.n_tested} "
            f"audible={result.n_audible} area_km2={result.area_km2} "
            f"elapsed_s={result.elapsed_s} -> {result.outfile}"
        )
        for result in sorted(results, key=lambda r: r.gain)
    ]
    if fit is not None:
        summary_lines.append(
            f"fit best_omni={fit.best_omni} F1={fit.max_fbeta:.4f} "
            f"precision={fit.best_precision:.4f} recall={fit.best_recall:.4f}"
        )
        summary_lines.append(f"fit plot={fit.plot_path}")
        summary_lines.append(
            f"canonical fits -> {fit.csv_path} (populate via fit_3d_active_space.py)"
        )
    append_aam_run_summary(summary_lines)


def print_fit(fit: FitResult) -> None:
    log("=== fit (annotations) ===")
    log(f"  best_omni={fit.best_omni}  F1={fit.max_fbeta:.4f}  "
        f"precision={fit.best_precision:.4f}  recall={fit.best_recall:.4f}")
    log(f"  plot -> {fit.plot_path}")
    log(f"  canonical fits -> {fit.csv_path} (use fit_3d_active_space.py to populate)")


def main() -> None:
    args = parse_args()
    wall_start = time.perf_counter()

    site_dir, mic, study_area = initialize_site(args)
    ambience = load_ambience(args.unit, args.site, args.year)
    gen = build_generator(site_dir, study_area, ambience, mic, args.model)

    out_dir = output_dir(site_dir, args)
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
        model=args.model,
    )

    ncpu = args.cpus or len(args.gains)
    if args.fit:
        log(
            f"fit mode: gains {args.gains[0]}–{args.gains[-1]} dB "
            f"({len(args.gains)} values, step {args.omni_step:g} dB), beta={args.beta}"
        )

    results = run_all_gains(args.gains, ctx, ncpu)
    print_results(results, time.perf_counter() - wall_start)

    fit: FitResult | None = None
    if args.fit:
        valid_points = load_valid_points(site_dir, study_area, args)
        fit = run_fit(results, valid_points, site_dir, args)
        print_fit(fit)

    if args.model is AcousticModel.AAM:
        append_aam_site_log(site_dir, results, fit)


if __name__ == "__main__":
    main()
