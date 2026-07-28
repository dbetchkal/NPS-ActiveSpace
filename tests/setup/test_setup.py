"""Integration and contract tests for ``nps_active_space.setup``."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from affine import Affine
from shapely.geometry import Point, box

from nps_active_space.setup import NMSIM_DST_CRS, setup_site
from nps_active_space.setup.elevation import NODATA_INT16, create_elevation_tif, write_gridfloat
from nps_active_space.setup.site_decoder import parse_sit_coords_line
from nps_active_space.setup.site_writer import create_site_file
from nps_active_space.utils.computation import NMSIM_bbox_utm, coords_to_utm
from nps_active_space.utils.enums import SourceElevationUnits
from nps_active_space.utils.helpers import get_deployment, load_studyarea
from fixtures import (
    DENABULL_BOUNDS_4269,
    DENABULL_MIC_COORD,
    DENABULL_STUDYAREA_NE,
    DENABULL_STUDYAREA_SW,
    EXAMPLE_PROJECT_DIR,
    STUDY_BOUNDS_4269,
    write_source_dem,
)

class TestExampleDenatrlaArtifacts:
    """Contract tests against checked-in example deployment data."""

    def test_study_area_and_mic_crs_align(self):
        study_area = load_studyarea(str(EXAMPLE_PROJECT_DIR), "DENA", "TRLA", 2025)
        mic = get_deployment(str(EXAMPLE_PROJECT_DIR), "DENA", "TRLA", 2025, elevation=False)

        assert study_area.crs.to_epsg() == 4269
        assert mic.crs == NMSIM_bbox_utm(study_area)

    def test_elevation_tif_is_nad83_geographic_not_utm(self):
        tif_path = (
            EXAMPLE_PROJECT_DIR
            / "DENATRLA"
            / "Input_Data"
            / "01_ELEVATION"
            / "elevation_m_nad83_utm6.tif"
        )
        with rasterio.open(tif_path) as ds:
            assert ds.crs.to_epsg() == 4269
            # Geographic CRS uses degree cell sizes, not meter-based UTM.
            assert abs(ds.transform.a) < 1.0


class TestCreateElevationUnits:
    @pytest.mark.parametrize(
        ("source_units", "source_elevation", "expected_mean"),
        [
            (SourceElevationUnits.FEET, 5000.0, 5000.0 * 0.3048),
            (SourceElevationUnits.METERS, 1000.0, 1000.0),
        ],
    )
    def test_converts_source_elevation_to_meters(
        self,
        tmp_path,
        source_units,
        source_elevation,
        expected_mean,
    ):
        dem_path = tmp_path / "source_dem.tif"
        write_source_dem(dem_path, bounds_4269=STUDY_BOUNDS_4269, elevation=source_elevation)
        study_area = gpd.GeoDataFrame(
            {"Unit": ["T"], "Site": ["0001"], "Year": [2025]},
            geometry=[box(-136.1, 58.5, -135.9, 58.7)],
            crs="EPSG:4326",
        )
        output_tif = tmp_path / "elevation_m.tif"

        create_elevation_tif(
            source_dem=dem_path,
            study_area_wgs84=study_area,
            dst_crs=NMSIM_DST_CRS,
            output_tif=output_tif,
            nodata_int16=NODATA_INT16,
            source_elevation_units=source_units,
        )

        with rasterio.open(output_tif) as ds:
            values = ds.read(1, masked=True)
            assert float(values.mean()) == pytest.approx(expected_mean, abs=1.0)


class TestSetupSiteIntegration:
    def test_mic_coords_round_trip_through_sit_file(self, tmp_path):
        project_dir = tmp_path
        unit, site, year = "TEST", "0001", 2025
        site_dir = project_dir / f"{unit}{site}"

        mic_coord = (-136.0, 58.6)
        studyarea_sw = (-136.1, 58.5)
        studyarea_ne = (-135.9, 58.7)

        dem_path = tmp_path / "source_dem.tif"
        write_source_dem(dem_path, bounds_4269=STUDY_BOUNDS_4269, elevation=5000.0)

        result = setup_site(
            site_dir=site_dir,
            unit=unit,
            site=site,
            year=year,
            mic_coord=mic_coord,
            studyarea_sw=studyarea_sw,
            studyarea_ne=studyarea_ne,
            source_dem=dem_path,
        )

        mic = get_deployment(str(project_dir), unit, site, year, elevation=False)
        assert mic.lon == pytest.approx(mic_coord[0], abs=1e-4)
        assert mic.lat == pytest.approx(mic_coord[1], abs=1e-4)

        sit_lines = result.sit_path.read_text().splitlines()
        _, _, height = parse_sit_coords_line(sit_lines[2])
        assert height == pytest.approx(1.5)

        flt_path = result.elevation_tif_path.with_suffix(".flt")
        assert flt_path.exists()
        assert flt_path.with_suffix(".hdr").exists()

        with rasterio.open(result.elevation_tif_path) as ds:
            assert ds.crs.to_epsg() == 4269

        saved_study_area = gpd.read_file(result.study_area_path)
        assert saved_study_area.crs.to_epsg() == 4269

    def test_mic_coords_round_trip_when_mic_and_project_utm_zones_differ(self, tmp_path):
        """Regression: .sit easting/northing must use the project UTM zone, not the mic's local zone."""
        project_dir = tmp_path
        unit, site, year = "DENA", "BULL", 2019
        site_dir = project_dir / f"{unit}{site}"

        mic_zone, _ = coords_to_utm(lat=DENABULL_MIC_COORD[1], lon=DENABULL_MIC_COORD[0])
        project_zone, _ = coords_to_utm(
            lat=DENABULL_STUDYAREA_SW[1], lon=DENABULL_STUDYAREA_SW[0]
        )
        assert mic_zone != project_zone

        dem_path = tmp_path / "source_dem.tif"
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

        mic = get_deployment(str(project_dir), unit, site, year, elevation=False)
        assert mic.lon == pytest.approx(DENABULL_MIC_COORD[0], abs=1e-4)
        assert mic.lat == pytest.approx(DENABULL_MIC_COORD[1], abs=1e-4)

    def test_legacy_mic_local_zone_sit_raises_with_project_setup_command(self, tmp_path):
        project_dir = tmp_path
        unit, site, year = "DENA", "BULL", 2019
        site_dir = project_dir / f"{unit}{site}"

        dem_path = tmp_path / "source_dem.tif"
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

        with pytest.raises(ValueError, match="legacy mic-local UTM") as exc_info:
            get_deployment(str(project_dir), unit, site, year, elevation=False)

        message = str(exc_info.value)
        assert "project_setup" in message
        assert "epsg:26906" in message
        assert "epsg:26905" in message
        assert f"-u {unit} -s {site} -y {year}" in message
        assert (
            f"--studyarea-sw {DENABULL_STUDYAREA_SW[0]:.8f} {DENABULL_STUDYAREA_SW[1]:.8f}" in message
        )
        assert (
            f"--studyarea-ne {DENABULL_STUDYAREA_NE[0]:.8f} {DENABULL_STUDYAREA_NE[1]:.8f}" in message
        )

    def test_write_gridfloat_encodes_geometry_and_nodata(self, tmp_path):
        dst_int16 = np.array([[100, 200], [300, 32767]], dtype=np.int16)
        dst_transform = Affine(0.0004, 0.0, -136.1, 0.0, -0.0004, 58.7)
        output_base = tmp_path / "elevation_m_nad83_utm6"
        nodata_int16 = np.int16(32767)
        gridfloat_nodata = np.float32(-32768.0)

        write_gridfloat(
            output_base=output_base,
            dst_int16=dst_int16,
            dst_transform=dst_transform,
            dst_width=2,
            dst_height=2,
            nodata_int16=nodata_int16,
            gridfloat_nodata=gridfloat_nodata,
        )

        hdr_lines = {
            line.split()[0]: line.split()[1]
            for line in output_base.with_suffix(".hdr").read_text().splitlines()
        }
        assert int(hdr_lines["ncols"]) == 2
        assert int(hdr_lines["nrows"]) == 2
        assert float(hdr_lines["NODATA_value"]) == pytest.approx(float(gridfloat_nodata))
        assert float(hdr_lines["xllcorner"]) == pytest.approx(dst_transform.c)
        assert float(hdr_lines["yllcorner"]) == pytest.approx(dst_transform.f + 2 * dst_transform.e)
        assert float(hdr_lines["cellsize"]) == pytest.approx(dst_transform.a)

        flt = np.fromfile(output_base.with_suffix(".flt"), dtype=np.float32).reshape(2, 2)
        assert flt[0, 0] == pytest.approx(100.0)
        assert flt[0, 1] == pytest.approx(200.0)
        assert flt[1, 0] == pytest.approx(300.0)
        assert flt[1, 1] == pytest.approx(float(gridfloat_nodata))

