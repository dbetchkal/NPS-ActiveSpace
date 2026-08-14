"""CLI for estimating and fitting clock drift between NVSPL and causal track data."""

from __future__ import annotations

import logging
from argparse import ArgumentParser

import iyore
import numpy as np
import pandas as pd

import nps_active_space.utils.config as cfg
from nps_active_space.ground_truthing.load_tracks import load_tracks
from nps_active_space.utils.clock_drift import (
    DEFAULT_DRIFT_SEC_PER_DAY,
    DEFAULT_POST_RESET_DRIFT_SEC,
    DEFAULT_SEASON_START_DRIFT_SEC,
    ClockDriftFixer,
    default_clock_drift_file_path,
    infer_correction_period,
    write_constant_drift_csv,
)
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.helpers import get_deployment, get_logger, load_studyarea
from nps_active_space.utils.models import Nvspl

CLOCK_DRIFT_MODULE_LOGGER = "nps_active_space.utils.clock_drift"


def attach_module_logging(script_logger: logging.Logger, module_name: str) -> None:
    """Route a library module's log records to the script logger's console handlers."""
    module_logger = logging.getLogger(module_name)
    module_logger.handlers = script_logger.handlers
    module_logger.setLevel(script_logger.level)
    module_logger.propagate = False


def parse_indices(indices_str: str) -> list[int]:
    """Parse a comma-separated list of integer indices."""
    return [int(i.strip()) for i in indices_str.split(",")]


def parse_maintenance_times(times_str: str | None) -> list[pd.Timestamp]:
    """Parse optional comma-separated maintenance visit datetimes."""
    if not times_str:
        return []
    return [pd.Timestamp(t.strip()) for t in times_str.split(",")]


def format_drift_status(drift: float) -> str:
    """Describe a single drift estimate for summary output."""
    if drift is None or (isinstance(drift, float) and np.isnan(drift)):
        return "invalid"
    return "valid"


def print_drift_summary(times: pd.DatetimeIndex, drifts: np.ndarray) -> None:
    """Print a table of daily drift estimates with plot index labels."""
    label = 0
    rows: list[tuple[str, str, str, str]] = []
    for anchor, drift in zip(times, drifts, strict=True):
        status = format_drift_status(drift)
        if status == "valid":
            index_label = str(label)
            drift_str = f"{drift:.2f}"
            label += 1
        else:
            index_label = "-"
            drift_str = "NaN" if drift is None or np.isnan(drift) else str(drift)
        rows.append((index_label, str(anchor), drift_str, status))

    if not rows:
        print("No daily drift estimates were produced.")
        return

    def column_width(header: str, values: list[str]) -> int:
        return max([len(header), *(len(value) for value in values)])

    col_widths = (
        column_width("index", [r[0] for r in rows]),
        column_width("anchor", [r[1] for r in rows]),
        column_width("drift_sec", [r[2] for r in rows]),
        column_width("status", [r[3] for r in rows]),
    )
    header = (
        f"{'index':>{col_widths[0]}}  "
        f"{'anchor':<{col_widths[1]}}  "
        f"{'drift_sec':>{col_widths[2]}}  "
        f"{'status':<{col_widths[3]}}"
    )
    print(header)
    print("-" * len(header))
    for index_label, anchor, drift_str, status in rows:
        print(
            f"{index_label:>{col_widths[0]}}  "
            f"{anchor:<{col_widths[1]}}  "
            f"{drift_str:>{col_widths[2]}}  "
            f"{status:<{col_widths[3]}}"
        )


def build_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description="Estimate and fit clock drift between NVSPL and causal track data.",
    )
    parser.add_argument(
        "-e", "--environment", required=True,
        help="The configuration environment to run the script in.",
    )
    parser.add_argument("-u", "--unit", required=True, help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", required=True, help="Four letter site code. E.g. TRLA")
    parser.add_argument("-y", "--year", type=int, required=True, help="Four digit year. E.g. 2025")
    parser.add_argument(
        "-t", "--track-source",
        default=TrackSource.ADSB,
        type=TrackSource,
        choices=list(TrackSource),
        help="Track source (default ADSB; primary use case for clock drift correction).",
    )
    parser.add_argument(
        "--plot-dir",
        help="Directory for QC plots (default: <site_dir>/clock_drift_qc/).",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not call plt.show() (useful for CI/automation).",
    )
    parser.add_argument(
        "--max-drift-minutes",
        type=float,
        default=5.0,
        help="Maximum expected clock drift magnitude in minutes (default: 5).",
    )
    parser.add_argument(
        "--method",
        default="correlation",
        choices=["correlation", "constant_drift"],
        help=(
            "Workflow to use (default: correlation). 'correlation' estimates daily drift from "
            "NVSPL/track cross-correlation and requires --fit --indices to write the CSV. "
            "'constant_drift' writes a piecewise-linear CSV from assumed drift rates and optional "
            "maintenance visits."
        ),
    )
    parser.add_argument(
        "--drift-sec-per-day",
        type=float,
        default=DEFAULT_DRIFT_SEC_PER_DAY,
        help=(
            "Only used by --method constant_drift. Assumed ADS-B clock drift rate in sec/day "
            f"(default: {DEFAULT_DRIFT_SEC_PER_DAY})."
        ),
    )
    parser.add_argument(
        "--season-start-drift-sec",
        type=float,
        default=DEFAULT_SEASON_START_DRIFT_SEC,
        help=(
            "Only used by --method constant_drift. Drift at season start before the first "
            f"maintenance visit (default: {DEFAULT_SEASON_START_DRIFT_SEC})."
        ),
    )
    parser.add_argument(
        "--post-reset-drift-sec",
        type=float,
        default=DEFAULT_POST_RESET_DRIFT_SEC,
        help=(
            "Only used by --method constant_drift. Drift immediately after each maintenance "
            f"visit (default: {DEFAULT_POST_RESET_DRIFT_SEC}, based on DENATRLA history)."
        ),
    )
    parser.add_argument(
        "--fit",
        action="store_true",
        help="Fit linear drift lines and write the clock drift CSV (correlation method only).",
    )
    parser.add_argument(
        "--indices",
        help="Comma-separated credible point indices from the estimate plot (required with --fit).",
    )
    parser.add_argument(
        "--maintenance-times",
        help=(
            "Comma-separated maintenance visit datetimes (ISO format, e.g. "
            "2025-07-18T11:13:00). Used by constant_drift and optional piecewise correlation fits."
        ),
    )
    return parser


def run_constant_drift(
    project_dir: str,
    unit: str,
    site: str,
    year_str: str,
    track_source: TrackSource,
    nvspl_dates: list[str],
    tracks: pd.DataFrame,
    maintenance_times: list[pd.Timestamp],
    drift_sec_per_day: float,
    season_start_drift_sec: float,
    post_reset_drift_sec: float,
) -> str:
    """Build and write a constant-drift correction CSV."""
    start_dt, end_dt = infer_correction_period(nvspl_dates, tracks)
    clock_drift_file = default_clock_drift_file_path(
        project_dir, unit, site, year_str, track_source,
    )
    write_constant_drift_csv(
        clock_drift_file=clock_drift_file,
        start_dt=start_dt,
        end_dt=end_dt,
        maintenance_times=maintenance_times,
        drift_sec_per_day=drift_sec_per_day,
        season_start_drift_sec=season_start_drift_sec,
        post_reset_drift_sec=post_reset_drift_sec,
    )
    return clock_drift_file


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.method == "correlation":
        if args.fit and not args.indices:
            parser.error("--indices is required when using --fit")
    elif args.method == "constant_drift":
        if args.fit or args.indices:
            parser.error("--fit and --indices are only used with --method correlation")

    cfg.initialize(environment=args.environment)
    logger = get_logger("CLOCK-DRIFT")
    attach_module_logging(logger, CLOCK_DRIFT_MODULE_LOGGER)

    project_dir = cfg.read("project", "dir")
    site_dir = f"{project_dir}/{args.unit}{args.site}"
    plot_dir = args.plot_dir or f"{site_dir}/clock_drift_qc"
    show_plots = not args.no_show
    max_clock_drift = pd.Timedelta(minutes=args.max_drift_minutes)
    year_str = str(args.year)
    maintenance_times = parse_maintenance_times(args.maintenance_times)

    logger.info(
        f"Beginning clock drift workflow for {args.unit}{args.site}{year_str} "
        f"({args.track_source}), method={args.method}..."
    )

    archive = iyore.Dataset(cfg.read("data", "nvspl_archive"))
    microphone = get_deployment(project_dir, args.unit, args.site, year_str)
    study_area = load_studyarea(project_dir, args.unit, args.site, year_str)

    nvspl_dates = sorted({
        f"{e.year}-{e.month}-{e.day}"
        for e in archive.nvspl(unit=args.unit, site=args.site, year=year_str)
    })
    if not nvspl_dates:
        raise AssertionError(
            f"No NVSPL data found in archive {cfg.read('data', 'nvspl_archive')}"
        )

    logger.info("Loading tracks...")
    loaded = load_tracks(
        args.track_source,
        start_date=nvspl_dates[0],
        end_date=nvspl_dates[-1],
        study_area=study_area,
        microphone=microphone,
    )
    tracks = loaded.tracks
    if tracks.empty:
        raise AssertionError("No tracks loaded; is your track source correct?")

    logger.info(
        f"Loaded {tracks['track_id'].nunique()} tracks "
        f"({len(tracks)} points) from {nvspl_dates[0]} to {nvspl_dates[-1]}"
    )

    if args.method == "constant_drift":
        clock_drift_file = run_constant_drift(
            project_dir=project_dir,
            unit=args.unit,
            site=args.site,
            year_str=year_str,
            track_source=args.track_source,
            nvspl_dates=nvspl_dates,
            tracks=tracks,
            maintenance_times=maintenance_times,
            drift_sec_per_day=args.drift_sec_per_day,
            season_start_drift_sec=args.season_start_drift_sec,
            post_reset_drift_sec=args.post_reset_drift_sec,
        )
        print(f"\nWrote clock drift file: {clock_drift_file}")
        print(
            "Verify alignment in the ground truthing GUI before relying on this correction."
        )
        return

    logger.info("Loading NVSPL for full deployment year...")
    nvspl_files = [
        e.path
        for e in archive.nvspl(unit=args.unit, site=args.site, year=year_str)
    ]
    nvspl = Nvspl(nvspl_files)
    logger.info(f"Loaded {len(nvspl_files)} NVSPL files")

    if loaded.faa_path is None or loaded.faa_corrections_path is None:
        raise ValueError(
            f"FAA paths are required for clock drift with track source {args.track_source}"
        )

    fixer = ClockDriftFixer(
        project_dir=project_dir,
        unit=args.unit,
        site=args.site,
        year=year_str,
        pts=tracks,
        nvspl=nvspl,
        database_type=args.track_source,
        faa_path=loaded.faa_path,
        faa_corrections_path=loaded.faa_corrections_path,
        plot_dir=plot_dir,
    )

    times, drifts = fixer.drift_time_series(
        max_clock_drift=max_clock_drift,
        show_plots=show_plots,
    )

    print("\nDaily clock drift estimates:")
    print_drift_summary(times, drifts)

    if args.fit:
        indices = parse_indices(args.indices)
        fixer.fit_drift_lines(
            indices_to_use=indices,
            maintenance_times=maintenance_times,
            show_plots=show_plots,
        )
        print(f"\nWrote clock drift file: {fixer.default_clock_drift_file}")
    else:
        print(
            "\nReview the QC plots, note credible index labels, then re-run with "
            f"--fit --indices <comma-separated indices>."
        )
        if plot_dir:
            print(f"Plots saved to: {plot_dir}")


if __name__ == "__main__":
    main()
