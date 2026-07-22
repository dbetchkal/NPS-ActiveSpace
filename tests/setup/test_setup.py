"""Integration and contract tests for ``nps_active_space.setup``."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from affine import Affine
from shapely.geometry import box

from nps_active_space.setup import (
    NMSIM_DST_CRS,
    parse_sit_coords_line,
    setup_site,
)
from nps_active_space.setup.elevation import NODATA_INT16, create_elevation_tif, write_gridfloat
from nps_active_space.utils.computation import NMSIM_bbox_utm
from nps_active_space.utils.enums import SourceElevationUnits
from nps_active_space.utils.helpers import get_deployment, load_studyarea

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_PROJECT_DIR = REPO_ROOT / "example_data" / "site_projects"
STUDY_BOUNDS_4269 = (-136.2, 58.4, -135.8, 58.8)


def _write_source_dem(
    path: Path,
    bounds_4269: tuple[float, float, float, float],
    elevation: float,
    cellsize_deg: float = 0.02,
) -> None:
    """Write a minimal single-band DEM in NAD83 geographic (EPSG:4269)."""
    minx, miny, maxx, maxy = bounds_4269
    width = max(2, int(np.ceil((maxx - minx) / cellsize_deg)))
    height = max(2, int(np.ceil((maxy - miny) / cellsize_deg)))
    transform = Affine(cellsize_deg, 0.0, minx, 0.0, -cellsize_deg, maxy)
    data = np.full((height, width), elevation, dtype=np.float32)

    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs=NMSIM_DST_CRS,
        transform=transform,
        nodata=-9999.0,
    ) as ds:
        ds.write(data, 1)


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
    def test_converts_feet_source_to_meters(self, tmp_path):
        dem_path = tmp_path / "source_dem.tif"
        _write_source_dem(dem_path, bounds_4269=STUDY_BOUNDS_4269, elevation=5000.0)
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
            source_elevation_units=SourceElevationUnits.FEET,
        )

        with rasterio.open(output_tif) as ds:
            values = ds.read(1, masked=True)
            assert float(values.mean()) == pytest.approx(5000.0 * 0.3048, abs=1.0)

    def test_passthrough_meters_source(self, tmp_path):
        dem_path = tmp_path / "source_dem.tif"
        _write_source_dem(dem_path, bounds_4269=STUDY_BOUNDS_4269, elevation=1000.0)
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
            source_elevation_units=SourceElevationUnits.METERS,
        )

        with rasterio.open(output_tif) as ds:
            values = ds.read(1, masked=True)
            assert float(values.mean()) == pytest.approx(1000.0, abs=1.0)


class TestSetupSiteIntegration:
    def test_mic_coords_round_trip_through_sit_file(self, tmp_path):
        project_dir = tmp_path
        unit, site, year = "TEST", "0001", 2025
        site_dir = project_dir / f"{unit}{site}"

        mic_coord = (-136.0, 58.6)
        studyarea_sw = (-136.1, 58.5)
        studyarea_ne = (-135.9, 58.7)

        dem_path = tmp_path / "source_dem.tif"
        _write_source_dem(dem_path, bounds_4269=STUDY_BOUNDS_4269, elevation=5000.0)

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

    def test_gridfloat_dimensions_match_header(self, tmp_path):
        dst_int16 = np.array([[100, 200], [300, 32767]], dtype=np.int16)
        dst_transform = Affine(0.0004, 0.0, -136.1, 0.0, -0.0004, 58.7)
        output_base = tmp_path / "elevation_m_nad83_utm6"

        write_gridfloat(
            output_base=output_base,
            dst_int16=dst_int16,
            dst_transform=dst_transform,
            dst_width=2,
            dst_height=2,
            nodata_int16=np.int16(32767),
            gridfloat_nodata=np.float32(-32768.0),
        )

        flt = np.fromfile(output_base.with_suffix(".flt"), dtype=np.float32)
        hdr_lines = {
            line.split()[0]: line.split()[1]
            for line in output_base.with_suffix(".hdr").read_text().splitlines()
        }
        assert flt.size == int(hdr_lines["ncols"]) * int(hdr_lines["nrows"])
