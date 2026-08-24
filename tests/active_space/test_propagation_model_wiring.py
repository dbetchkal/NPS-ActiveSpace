"""Tests for PropagationModel wiring in ActiveSpaceGenerator."""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
from nps_active_space.active_space.propagation_model import NMSIM_BAND_COLUMNS


class _StubPropagationModel:
    max_points_per_run = 7

    def prepare_site(self, dem_src, study_area, mic, *, project_dem=True, suffix=""):
        return {"dem_file": dem_src}

    def predict(self, site, source_pts, omni_source, altitude_m, job_name, heading=None):
        return pd.DataFrame({
            "Xpos": source_pts.geometry.x.values,
            "Ypos": source_pts.geometry.y.values,
            "Zpos": source_pts.geometry.z.values,
            "A": [50.0] * len(source_pts),
            **{col: [40.0] * len(source_pts) for col in NMSIM_BAND_COLUMNS},
        })


class TestPropagationModelWiring:
    def test_generator_accepts_custom_propagation_model(self, tmp_path) -> None:
        dem = tmp_path / "dem.tif"
        dem.write_bytes(b"")
        study_area = gpd.GeoDataFrame(
            geometry=[Point(0, 0).buffer(1)],
            crs="EPSG:4326",
        )
        stub = _StubPropagationModel()
        gen = ActiveSpaceGenerator(
            NMSIM=None,
            study_area=study_area,
            root_dir=str(tmp_path),
            dem_src=str(dem),
            ambience=30.0,
            propagation_model=stub,
        )
        assert gen.propagation_model.max_points_per_run == 7

    def test_preprocess_uses_model_batch_cap(self, tmp_path) -> None:
        dem = tmp_path / "dem.tif"
        dem.write_bytes(b"")
        crs = "EPSG:32606"
        study_area = gpd.GeoDataFrame(
            geometry=[Point(500000, 6000000).buffer(5000)],
            crs=crs,
        )
        stub = _StubPropagationModel()
        gen = ActiveSpaceGenerator(
            NMSIM=None,
            study_area=study_area,
            root_dir=str(tmp_path),
            dem_src=str(dem),
            ambience=30.0,
            propagation_model=stub,
        )
        pts = gpd.GeoDataFrame(
            geometry=[Point(500000 + i * 10, 6000000) for i in range(20)],
            crs=crs,
        )
        region = Point(500000, 6000000).buffer(10000)
        tested = gpd.GeoDataFrame(columns=["audible", "geometry"], geometry="geometry", crs=crs)
        filtered = gen._preprocess_source_points(
            pts, region, tested, max_pts=gen.propagation_model.max_points_per_run,
        )
        assert len(filtered) == 7
