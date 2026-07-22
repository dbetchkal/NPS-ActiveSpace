"""Tests for legacy ``.sit`` migration."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import Point, box

from nps_active_space.setup import create_site_file, setup_site
from nps_active_space.setup.site_decoder import (
    diagnose_sit_coords,
    read_sit_file,
)
from nps_active_space.setup.site_migrate import (
    discover_deployment_sits,
    migrate_deployment_sit,
    migrate_project_sites,
    migrate_sit_file,
)
from nps_active_space.utils.computation import coords_to_utm
from nps_active_space.utils.helpers import get_deployment, load_studyarea
from fixtures import (
    DENABULL_BOUNDS_4269,
    DENABULL_MIC_COORD,
    DENABULL_STUDYAREA_NE,
    DENABULL_STUDYAREA_SW,
    STUDY_BOUNDS_4269,
    write_source_dem,
)

def _create_denabull_deployment(project_dir: Path, *, legacy_sit: bool) -> tuple[str, str, int]:
    unit, site, year = "DENA", "BULL", 2019
    site_dir = project_dir / f"{unit}{site}"
    dem_path = project_dir / "source_dem.tif"
    write_source_dem(dem_path, bounds_4269=DENABULL_BOUNDS_4269, elevation=5000.0)

    setup_site(
        site_dir=site_dir,
        unit=unit,
        site=site,
        year=year,
        mic_coord=DENABULL_MIC_COORD,
        studyarea_sw=DENABULL_STUDYAREA_SW,
        studyarea_ne=DENABULL_STUDYAREA_NE,
        source_dem=dem_path,
    )

    if legacy_sit:
        mic_zone, _ = coords_to_utm(lat=DENABULL_MIC_COORD[1], lon=DENABULL_MIC_COORD[0])
        mic_utm = gpd.GeoSeries([Point(DENABULL_MIC_COORD)], crs="EPSG:4326").to_crs(mic_zone)
        create_site_file(
            site_dir=site_dir,
            unit=unit,
            site=site,
            year=year,
            easting_m=mic_utm.x.item(),
            northing_m=mic_utm.y.item(),
        )

    return unit, site, year


class TestDiagnoseSitCoords:
    def test_legacy_and_ok_statuses(self, tmp_path):
        unit, site, year = _create_denabull_deployment(tmp_path, legacy_sit=True)
        study_area = load_studyarea(str(tmp_path), unit, site, year)
        sit_path = tmp_path / f"{unit}{site}" / "Input_Data" / "05_SITES" / f"{unit}{site}{year}.sit"
        sit_contents = read_sit_file(sit_path)

        status, lon, lat, project_utm, decoded_utm = diagnose_sit_coords(
            sit_contents.easting_m,
            sit_contents.northing_m,
            study_area,
        )
        assert status == "legacy"
        assert decoded_utm != project_utm
        assert lon == pytest.approx(DENABULL_MIC_COORD[0], abs=1e-4)
        assert lat == pytest.approx(DENABULL_MIC_COORD[1], abs=1e-4)

        migrate_sit_file(
            sit_path,
            study_area,
            unit=unit,
            site=site,
            year=year,
        )
        migrated = read_sit_file(sit_path)
        status, _, _, _, decoded_utm = diagnose_sit_coords(
            migrated.easting_m,
            migrated.northing_m,
            study_area,
        )
        assert status == "ok"
        assert decoded_utm == project_utm


class TestMigrateSitFile:
    def test_migrates_legacy_sit_and_preserves_metadata(self, tmp_path):
        unit, site, year = _create_denabull_deployment(tmp_path, legacy_sit=True)
        study_area = load_studyarea(str(tmp_path), unit, site, year)
        sit_path = tmp_path / f"{unit}{site}" / "Input_Data" / "05_SITES" / f"{unit}{site}{year}.sit"
        before = read_sit_file(sit_path)

        result = migrate_sit_file(
            sit_path,
            study_area,
            unit=unit,
            site=site,
            year=year,
        )

        assert result.action == "migrated"
        assert result.decoded_utm != result.project_utm
        after = read_sit_file(sit_path)
        assert after.label == before.label
        assert after.height_agl_m == before.height_agl_m
        assert after.dem_flt_path == before.dem_flt_path
        assert (after.easting_m, after.northing_m) != (before.easting_m, before.northing_m)

        mic = get_deployment(str(tmp_path), unit, site, year, elevation=False)
        assert mic.lon == pytest.approx(DENABULL_MIC_COORD[0], abs=1e-4)
        assert mic.lat == pytest.approx(DENABULL_MIC_COORD[1], abs=1e-4)

    def test_skips_project_zone_sit(self, tmp_path):
        unit, site, year = _create_denabull_deployment(tmp_path, legacy_sit=False)
        study_area = load_studyarea(str(tmp_path), unit, site, year)
        sit_path = tmp_path / f"{unit}{site}" / "Input_Data" / "05_SITES" / f"{unit}{site}{year}.sit"
        before = read_sit_file(sit_path)

        result = migrate_sit_file(
            sit_path,
            study_area,
            unit=unit,
            site=site,
            year=year,
        )

        assert result.action == "skipped_ok"
        after = read_sit_file(sit_path)
        assert after.easting_m == before.easting_m
        assert after.northing_m == before.northing_m

    def test_dry_run_does_not_write(self, tmp_path):
        unit, site, year = _create_denabull_deployment(tmp_path, legacy_sit=True)
        study_area = load_studyarea(str(tmp_path), unit, site, year)
        sit_path = tmp_path / f"{unit}{site}" / "Input_Data" / "05_SITES" / f"{unit}{site}{year}.sit"
        before = read_sit_file(sit_path)

        result = migrate_sit_file(
            sit_path,
            study_area,
            unit=unit,
            site=site,
            year=year,
            dry_run=True,
        )

        assert result.action == "dry_run_would_migrate"
        after = read_sit_file(sit_path)
        assert after.easting_m == before.easting_m
        assert after.northing_m == before.northing_m


class TestMigrateProjectSites:
    def test_discover_and_migrate_all(self, tmp_path):
        unit, site, year = _create_denabull_deployment(tmp_path, legacy_sit=True)

        discovered = discover_deployment_sits(tmp_path)
        assert discovered == [(unit, site, year, discovered[0][3])]

        results = migrate_project_sites(tmp_path)
        assert len(results) == 1
        assert results[0].action == "migrated"

        get_deployment(str(tmp_path), unit, site, year, elevation=False)

    def test_migrate_single_deployment(self, tmp_path):
        unit, site, year = _create_denabull_deployment(tmp_path, legacy_sit=True)

        result = migrate_deployment_sit(tmp_path, unit, site, year)
        assert result.action == "migrated"

    def test_missing_sit_returns_failed(self, tmp_path):
        result = migrate_deployment_sit(tmp_path, "DENA", "BULL", 2019)
        assert result.action == "failed"
        assert "not found" in (result.message or "").lower()

    def test_failed_coords_outside_study_area(self, tmp_path):
        unit, site, year = "TEST", "0001", 2025
        site_dir = tmp_path / f"{unit}{site}"
        dem_path = tmp_path / "source_dem.tif"
        write_source_dem(dem_path, bounds_4269=(-136.2, 58.4, -135.8, 58.8), elevation=1000.0)

        setup_site(
            site_dir=site_dir,
            unit=unit,
            site=site,
            year=year,
            mic_coord=(-136.0, 58.6),
            studyarea_sw=(-136.1, 58.5),
            studyarea_ne=(-135.9, 58.7),
            source_dem=dem_path,
        )

        study_area = load_studyarea(str(tmp_path), unit, site, year)
        sit_path = site_dir / "Input_Data" / "05_SITES" / f"{unit}{site}{year}.sit"
        create_site_file(
            site_dir=site_dir,
            unit=unit,
            site=site,
            year=year,
            easting_m=1.0,
            northing_m=1.0,
        )

        result = migrate_sit_file(
            sit_path,
            study_area,
            unit=unit,
            site=site,
            year=year,
        )
        assert result.action == "failed"
        assert "project_setup" in (result.message or "")
