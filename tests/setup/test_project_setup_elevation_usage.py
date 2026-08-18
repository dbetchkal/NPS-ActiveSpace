"""Tests for requiring project_setup elevation artifacts in active space generation."""

from __future__ import annotations

import shutil
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
from nps_active_space.setup import setup_site
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.setup.site_writer import create_site_dir
from nps_active_space.utils.computation import NMSIM_bbox_utm
from nps_active_space.utils.helpers import get_deployment, load_studyarea
from nps_active_space.utils.models import Microphone
from fixtures import (
    EXAMPLE_PROJECT_DIR,
    STUDY_BOUNDS_4269,
    write_project_setup_elevation_artifacts,
    write_source_dem,
)


class TestGetProjectSetupElevation:
    def test_returns_tif_and_flt_from_site_artifacts(self, tmp_path: Path):
        site_dir = tmp_path / "TESTSITE"
        expected_tif, expected_flt = write_project_setup_elevation_artifacts(site_dir)

        tif_path, flt_path = get_project_setup_elevation(site_dir)
        assert tif_path == expected_tif
        assert flt_path == expected_flt
        assert flt_path.with_suffix(".hdr").is_file()

    def test_raises_when_missing(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="project_setup"):
            get_project_setup_elevation(tmp_path)

    def test_example_denatrla_has_complete_meter_artifacts(self):
        site_dir = EXAMPLE_PROJECT_DIR / "DENATRLA"
        tif_path, flt_path = get_project_setup_elevation(site_dir)
        assert tif_path.name == "elevation_m_nad83_utm6.tif"
        assert flt_path.name == "elevation_m_nad83_utm6.flt"

    def test_example_denatrla_set_dem_uses_project_setup_artifacts(self, tmp_path: Path):
        site_dir = tmp_path / "DENATRLA"
        shutil.copytree(EXAMPLE_PROJECT_DIR / "DENATRLA", site_dir)
        nmsim_exe = tmp_path / "nmsim.exe"
        nmsim_exe.write_text("stub")

        study_area = load_studyarea(str(EXAMPLE_PROJECT_DIR), "DENA", "TRLA", 2025)
        mic = get_deployment(str(EXAMPLE_PROJECT_DIR), "DENA", "TRLA", 2025, elevation=False)
        generator = ActiveSpaceGenerator(
            NMSIM=str(nmsim_exe),
            study_area=study_area,
            root_dir=str(site_dir),
            ambience=pd.Series({"1000": 40.0}),
        )
        generator.set_dem(mic)

        assert generator._dem_file.endswith("elevation_m_nad83_utm6.tif")
        assert generator._flt_file.endswith("elevation_m_nad83_utm6.flt")
        assert Path(generator._flt_file).is_file()
        sit_path = site_dir / "Input_Data/05_SITES/DENATRLA2025.sit"
        assert "elevation_m_nad83_utm6.flt" in sit_path.read_text()


class TestActiveSpaceGeneratorSetDem:
    def test_uses_project_setup_elevation_from_setup_site(self, tmp_path: Path):
        project_dir = tmp_path / "project"
        site_dir = project_dir / "TESTSITE"
        source_dem = tmp_path / "global_dem.tif"
        write_source_dem(source_dem, bounds_4269=STUDY_BOUNDS_4269, elevation=5000.0)
        nmsim_exe = tmp_path / "nmsim.exe"
        nmsim_exe.write_text("stub")

        setup_site(
            site_dir=site_dir,
            unit="TEST",
            site="SITE",
            year=2025,
            mic_coord=(-136.05, 58.55),
            studyarea_sw=(-136.15, 58.45),
            studyarea_ne=(-135.85, 58.75),
            source_dem=source_dem,
        )

        study_area = load_studyarea(str(project_dir), "TEST", "SITE", 2025)
        mic = get_deployment(str(project_dir), "TEST", "SITE", 2025, elevation=False)
        generator = ActiveSpaceGenerator(
            NMSIM=str(nmsim_exe),
            study_area=study_area,
            root_dir=str(site_dir),
            ambience=pd.Series({"1000": 40.0}),
        )
        generator.set_dem(mic)

        assert "elevation_m_nad83_utm" in generator._dem_file
        assert Path(generator._flt_file).suffix == ".flt"
        assert Path(generator._dem_file).name != "elevation.tif"

    def test_raises_when_setup_elevation_missing(self, tmp_path: Path):
        site_dir = tmp_path / "TESTSITE"
        create_site_dir(site_dir)
        nmsim_exe = tmp_path / "nmsim.exe"
        nmsim_exe.write_text("stub")

        study_area = gpd.GeoDataFrame(
            {"Unit": ["TEST"], "Site": ["SITE"], "Year": [2025]},
            geometry=[box(*STUDY_BOUNDS_4269)],
            crs="EPSG:4269",
        )
        study_area.to_file(site_dir / "TESTSITE_study_area.shp")

        mic_crs = NMSIM_bbox_utm(study_area)
        mic = Microphone(name="TESTSITE2025", lat=58.55, lon=-136.05, z=100.0, crs=mic_crs)

        generator = ActiveSpaceGenerator(
            NMSIM=str(nmsim_exe),
            study_area=study_area,
            root_dir=str(site_dir),
            ambience=pd.Series({"1000": 40.0}),
        )

        with pytest.raises(FileNotFoundError, match="project_setup"):
            generator.set_dem(mic)
