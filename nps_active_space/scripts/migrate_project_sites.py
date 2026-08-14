"""
Rewrite legacy ``.sit`` microphone coordinates to the NMSIM project UTM zone.

Use this after upgrading to the project-zone ``.sit`` writer fix when on-disk deployments
were created with mic-local UTM easting/northing. This updates only the ``.sit`` file; it
does not re-clip elevation data or require ``data.dem`` in the config.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import nps_active_space.utils.config as cfg
from nps_active_space.setup.site_migrate import (
    SitMigrationResult,
    migrate_deployment_sit,
    migrate_project_sites,
)


def _print_result(result: SitMigrationResult) -> None:
    label = f"{result.unit} {result.site} {result.year}"
    match result.action:
        case "skipped_ok":
            print(f"ok       {label}  already project zone ({result.project_utm})")
        case "migrated":
            print(
                f"migrate  {label}  {result.decoded_utm} -> {result.project_utm}  "
                f"{result.sit_path}"
            )
        case "dry_run_would_migrate":
            print(
                f"dry-run  {label}  would migrate {result.decoded_utm} -> "
                f"{result.project_utm}"
            )
        case "failed":
            print(f"failed   {label}  {result.message or result.sit_path}")


def _summarize(results: list[SitMigrationResult]) -> int:
    migrated = sum(result.action == "migrated" for result in results)
    ok = sum(result.action == "skipped_ok" for result in results)
    dry_run = sum(result.action == "dry_run_would_migrate" for result in results)
    failed = sum(result.action == "failed" for result in results)
    print("---")
    print(f"Summary: {migrated} migrated, {ok} ok, {dry_run} dry-run, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rewrite legacy .sit files to use the NMSIM project UTM zone."
    )
    parser.add_argument(
        "-e",
        "--environment",
        required=True,
        help="The configuration environment to run the script in.",
    )
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--all", action="store_true", help="Migrate every deployment found.")
    target.add_argument("-u", "--unit", help="Four letter unit code. E.g. DENA")
    parser.add_argument("-s", "--site", help="Four letter site code. E.g. BULL")
    parser.add_argument("-y", "--year", type=int, help="Four digit year. E.g. 2019")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report legacy .sit files without writing changes.",
    )
    args = parser.parse_args()

    if args.unit and not (args.site and args.year):
        parser.error("-u requires -s and -y")
    if (args.site or args.year) and not args.unit:
        parser.error("-s and -y require -u")

    cfg.initialize(environment=args.environment)
    project_dir = Path(cfg.read("project", "dir"))

    if args.all:
        results = migrate_project_sites(project_dir, dry_run=args.dry_run)
    else:
        results = [
            migrate_deployment_sit(
                project_dir,
                args.unit,
                args.site,
                args.year,
                dry_run=args.dry_run,
            )
        ]

    for result in results:
        _print_result(result)

    sys.exit(_summarize(results))
