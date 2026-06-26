from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

from nps_active_space.utils.enums import TrackSource
from nps_active_space.viz.visualizer import Visualizer


def parse_deployment(name: str) -> tuple[str, str, int]:
    """Parse deployment name like DENATRLA2024 into unit, site, and year."""
    if len(name) < 9:
        raise argparse.ArgumentTypeError(
            f"deployment must be at least 9 characters "
            f"(4-char unit + site code + 4-digit year), got {len(name)}: {name!r}"
        )
    unit, site, year_str = name[:4], name[4:-4], name[-4:]
    if not site:
        raise argparse.ArgumentTypeError(
            f"deployment must include a site code between unit and year, got {name!r}"
        )
    if len(year_str) != 4 or not year_str.isdigit():
        raise argparse.ArgumentTypeError(
            f"deployment year must be 4 digits, got {year_str!r} in {name!r}"
        )
    return unit, site, int(year_str)


def parse_iso_date(value: str, *, arg_name: str) -> str:
    """Validate YYYY-MM-DD date string."""
    try:
        datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{arg_name}: expected YYYY-MM-DD date, got {value!r}"
        ) from None
    return value


def parse_existing_file(path: str, *, arg_name: str) -> str:
    """Validate that a CLI path refers to an existing file."""
    file_path = Path(path)
    if not file_path.is_file():
        raise argparse.ArgumentTypeError(f"{arg_name}: file not found: {path}")
    return path


def parse_max_tracks(value: str) -> int:
    """Validate positive integer for --max-tracks."""
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--max-tracks must be a positive integer, got {value!r}"
        ) from None
    if n <= 0:
        raise argparse.ArgumentTypeError(
            f"--max-tracks must be a positive integer, got {n}"
        )
    return n


def resolve_viz_plot_flags(
    *,
    active_space: bool = False,
    annotations: bool = False,
    audible_transits: bool = False,
    track_source: TrackSource | None = None,
    plot_all: bool = False,
    annotation_file: str | None = None,
    transits_pkl: str | None = None,
) -> tuple[bool, bool, bool, TrackSource | None]:
    """Map CLI flags to Visualizer layer toggles."""
    do_active = active_space or plot_all
    do_annotations = annotations or plot_all or annotation_file is not None
    do_transits = audible_transits or plot_all or transits_pkl is not None
    return do_active, do_annotations, do_transits, track_source


def resolve_track_source_args(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> TrackSource | None:
    """Apply deprecated vessel flags and validate track date options."""
    track_source = args.track_source
    if args.vessels:
        if track_source is not None and track_source is not TrackSource.AIS:
            parser.error(
                "--vessels is an alias for --track-source AIS; "
                "cannot combine with another --track-source"
            )
        if track_source is None:
            print(
                "Warning: -v/--vessels is deprecated; use --track-source AIS",
                file=sys.stderr,
            )
            track_source = TrackSource.AIS

    if (args.start_date or args.end_date) and track_source is None:
        parser.error("--start-date and --end-date require --track-source")

    if args.start_date and args.end_date and args.start_date > args.end_date:
        parser.error("--start-date must be on or before --end-date")

    return track_source


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "deployment",
        type=parse_deployment,
        help="Deployment name, e.g. DENATRLA2024",
    )
    parser.add_argument(
        "-e",
        "--environment",
        default="DENA_streamline",
        help="Config environment name, e.g. DENA_streamline. Default is 'DENA_streamline'",
    )
    parser.add_argument(
        "-s",
        "--active-space",
        action="store_true",
        help="If included, load and plot the active space.",
    )
    parser.add_argument("-g", "--gain", type=float, help="Active space gain, if not the default.")
    parser.add_argument(
        "-a",
        "--annotations",
        action="store_true",
        help="If included, load and plot annotations",
    )
    parser.add_argument(
        "-t",
        "--audible-transits",
        action="store_true",
        help="If included, load and plot audible transits",
    )
    parser.add_argument(
        "--track-source",
        type=TrackSource,
        choices=list(TrackSource),
        help="Load and plot causal tracks (GPS, ADSB, or AIS). Not included in --all.",
    )
    parser.add_argument(
        "-v",
        "--vessels",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Load everything, shorthand for --active-space --annotations --audible-transits",
    )
    parser.add_argument(
        "-m",
        "--max-tracks",
        type=parse_max_tracks,
        default=500,
        help="Maximum tracks to show (annotations, audible transits, or --track-source).",
    )
    parser.add_argument(
        "--annotation-file",
        type=lambda p: parse_existing_file(p, arg_name="--annotation-file"),
        help="Path to .geojson annotations (implies -a).",
    )
    parser.add_argument(
        "--transits-pkl",
        type=lambda p: parse_existing_file(p, arg_name="--transits-pkl"),
        help="Path to .pkl audible transits (implies -t).",
    )
    parser.add_argument(
        "--start-date",
        type=lambda d: parse_iso_date(d, arg_name="--start-date"),
        help="Track query start date (YYYY-MM-DD). Requires --track-source. Default: Jan 1 of deployment year.",
    )
    parser.add_argument(
        "--end-date",
        type=lambda d: parse_iso_date(d, arg_name="--end-date"),
        help="Track query end date (YYYY-MM-DD). Requires --track-source. Default: Dec 31 of deployment year.",
    )
    parser.add_argument(
        "--terraced",
        action="store_true",
        help="If included, render the active space as the terraced surface instead of contours.",
    )
    parser.add_argument(
        "--fill-layers",
        action="store_true",
        help="If included, fill the interior of each active space contour polygon.",
    )

    args = parser.parse_args()
    unit, site, year = args.deployment

    track_source = resolve_track_source_args(args, parser)
    do_active, do_annotations, do_transits, _ = resolve_viz_plot_flags(
        active_space=args.active_space,
        annotations=args.annotations,
        audible_transits=args.audible_transits,
        track_source=track_source,
        plot_all=args.all,
        annotation_file=args.annotation_file,
        transits_pkl=args.transits_pkl,
    )

    Visualizer(
        unit,
        site,
        year,
        args.environment,
        do_active,
        args.gain,
        do_annotations,
        do_transits,
        track_source,
        args.annotation_file,
        args.transits_pkl,
        args.start_date,
        args.end_date,
        args.terraced,
        args.fill_layers,
        args.max_tracks,
    )
