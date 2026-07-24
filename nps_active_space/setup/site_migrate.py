"""Rewrite legacy ``.sit`` files to use the NMSIM project UTM zone."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import geopandas as gpd

from nps_active_space.setup.site_decoder import (
    diagnose_sit_coords,
    format_project_setup_command,
    project_zone_utm_coords,
    read_sit_file,
    study_area_corner_coords_wgs84,
)
from nps_active_space.setup.site_writer import (
    deployment_sit_name,
    sit_file_path,
    write_site_file,
)
from nps_active_space.utils.helpers import load_studyarea

logger = logging.getLogger(__name__)

SitMigrationAction = Literal["skipped_ok", "migrated", "dry_run_would_migrate", "failed"]


@dataclass(frozen=True)
class SitMigrationResult:
    """Outcome of migrating one deployment ``.sit`` file."""

    unit: str
    site: str
    year: int
    sit_path: Path
    action: SitMigrationAction
    project_utm: str
    decoded_utm: str | None = None
    lon: float | None = None
    lat: float | None = None
    message: str | None = None


def discover_deployment_sits(project_dir: str | Path) -> list[tuple[str, str, int, Path]]:
    """
    Find canonical deployment ``.sit`` files under ``project_dir``.

    Returns sorted ``(unit, site, year, sit_path)`` tuples.
    """
    project_dir = Path(project_dir)
    deployments: list[tuple[str, str, int, Path]] = []

    for sit_path in sorted(project_dir.glob("*/Input_Data/05_SITES/*.sit")):
        site_folder = sit_path.parents[2].name
        stem = sit_path.stem
        if len(stem) < 5 or not stem[-4:].isdigit():
            logger.warning("Skipping non-canonical .sit file: %s", sit_path)
            continue

        year = int(stem[-4:])
        unit_site = stem[:-4]
        if unit_site != site_folder or len(unit_site) < 5:
            logger.warning("Skipping non-canonical .sit file: %s", sit_path)
            continue

        unit, site = unit_site[:4], unit_site[4:]
        deployments.append((unit, site, year, sit_path))

    return deployments


def migrate_sit_file(
    sit_path: str | Path,
    study_area: gpd.GeoDataFrame,
    *,
    unit: str,
    site: str,
    year: int,
    dry_run: bool = False,
) -> SitMigrationResult:
    """Rewrite a legacy ``.sit`` file to use the NMSIM project UTM zone when needed."""
    sit_path = Path(sit_path)
    sit_contents = read_sit_file(sit_path)
    status, lon, lat, project_utm, decoded_utm = diagnose_sit_coords(
        sit_contents.easting_m,
        sit_contents.northing_m,
        study_area,
    )

    if status == "ok":
        return SitMigrationResult(
            unit=unit,
            site=site,
            year=year,
            sit_path=sit_path,
            action="skipped_ok",
            project_utm=project_utm,
            decoded_utm=decoded_utm,
            lon=lon,
            lat=lat,
        )

    if status == "failed":
        studyarea_sw, studyarea_ne = study_area_corner_coords_wgs84(study_area)
        message = (
            "Could not decode .sit coordinates inside the study area. "
            "Re-run project_setup manually:\n"
            f"  {format_project_setup_command(unit, site, year, (lon, lat), studyarea_sw, studyarea_ne)}"
        )
        return SitMigrationResult(
            unit=unit,
            site=site,
            year=year,
            sit_path=sit_path,
            action="failed",
            project_utm=project_utm,
            decoded_utm=decoded_utm,
            lon=lon,
            lat=lat,
            message=message,
        )

    easting_m, northing_m = project_zone_utm_coords(lon, lat, study_area)
    if not dry_run:
        write_site_file(
            sit_path,
            easting_m,
            northing_m,
            sit_contents.height_agl_m,
            sit_contents.label,
            sit_contents.dem_flt_path,
        )

    return SitMigrationResult(
        unit=unit,
        site=site,
        year=year,
        sit_path=sit_path,
        action="dry_run_would_migrate" if dry_run else "migrated",
        project_utm=project_utm,
        decoded_utm=decoded_utm,
        lon=lon,
        lat=lat,
        message=f"{decoded_utm} -> {project_utm}",
    )


def migrate_deployment_sit(
    project_dir: str | Path,
    unit: str,
    site: str,
    year: int,
    *,
    dry_run: bool = False,
) -> SitMigrationResult:
    """Migrate the canonical ``.sit`` file for one deployment."""
    project_dir = Path(project_dir)
    site_dir = project_dir / f"{unit}{site}"
    sit_path = sit_file_path(site_dir, deployment_sit_name(unit, site, year))
    if not sit_path.is_file():
        return SitMigrationResult(
            unit=unit,
            site=site,
            year=year,
            sit_path=sit_path,
            action="failed",
            project_utm="",
            message=f"Microphone site file not found: {sit_path}",
        )

    study_area = load_studyarea(str(project_dir), unit, site, year)
    return migrate_sit_file(
        sit_path,
        study_area,
        unit=unit,
        site=site,
        year=year,
        dry_run=dry_run,
    )


def migrate_project_sites(
    project_dir: str | Path,
    *,
    deployments: list[tuple[str, str, int]] | None = None,
    dry_run: bool = False,
) -> list[SitMigrationResult]:
    """Migrate one or all deployments under ``project_dir``."""
    project_dir = Path(project_dir)
    discovered = discover_deployment_sits(project_dir)
    if deployments is None:
        selected = discovered
    else:
        requested = set(deployments)
        selected = [
            (unit, site, year, sit_path)
            for unit, site, year, sit_path in discovered
            if (unit, site, year) in requested
        ]
        missing = requested - {(unit, site, year) for unit, site, year, _ in selected}
        for unit, site, year in sorted(missing):
            logger.warning("No canonical .sit found for %s%s%s", unit, site, year)

    results: list[SitMigrationResult] = []
    for unit, site, year, sit_path in selected:
        study_area = load_studyarea(str(project_dir), unit, site, year)
        results.append(
            migrate_sit_file(
                sit_path,
                study_area,
                unit=unit,
                site=site,
                year=year,
                dry_run=dry_run,
            )
        )
    return results
