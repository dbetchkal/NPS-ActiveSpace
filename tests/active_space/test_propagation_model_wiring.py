"""Tests for PropagationModel wiring in ActiveSpaceGenerator."""

from __future__ import annotations

from unittest.mock import patch

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
from nps_active_space.active_space.propagation_model import NMSIM_BAND_COLUMNS


class _StubPropagationModel:
    max_points_per_run = 7
    predictions_subdir = "Output_Data/aam/predictions"

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
        study_area = gpd.GeoDataFrame(
            geometry=[Point(0, 0).buffer(1)],
            crs="EPSG:4326",
        )
        stub = _StubPropagationModel()
        gen = ActiveSpaceGenerator(
            NMSIM=None,
            study_area=study_area,
            root_dir=str(tmp_path),
            ambience=30.0,
            propagation_model=stub,
        )
        assert gen.propagation_model.max_points_per_run == 7

    def test_preprocess_uses_model_batch_cap(self, tmp_path) -> None:
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


class _EmptyPredictionPropagationModel:
    max_points_per_run = 7
    predictions_subdir = "Output_Data/aam/predictions"

    def prepare_site(self, dem_src, study_area, mic, *, project_dem=True, suffix=""):
        return {"dem_file": dem_src}

    def predict(self, site, source_pts, omni_source, altitude_m, job_name, heading=None):
        return pd.DataFrame()


class _PartialFailurePropagationModel:
    max_points_per_run = 7
    predictions_subdir = "Output_Data/aam/predictions"

    def prepare_site(self, dem_src, study_area, mic, *, project_dem=True, suffix=""):
        return {"dem_file": dem_src}

    def predict(self, site, source_pts, omni_source, altitude_m, job_name, heading=None):
        return pd.DataFrame({
            "Xpos": [source_pts.geometry.x.iloc[0]],
            "Ypos": [source_pts.geometry.y.iloc[0]],
            "Zpos": [source_pts.geometry.z.iloc[0]],
            "A": [50.0],
            **{col: [40.0] for col in NMSIM_BAND_COLUMNS},
        })


class TestMissingPredictionHandling:
    def test_run_propagation_model_raises_on_total_predict_failure(self, tmp_path) -> None:
        crs = "EPSG:32606"
        study_area = gpd.GeoDataFrame(
            geometry=[Point(500000, 6000000).buffer(5000)],
            crs=crs,
        )
        gen = ActiveSpaceGenerator(
            NMSIM=None,
            study_area=study_area,
            root_dir=str(tmp_path),
            ambience=pd.Series({"1000": 40.0, "12500": 40.0}),
            propagation_model=_EmptyPredictionPropagationModel(),
        )
        source_pts = gpd.GeoDataFrame(
            geometry=[
                Point(500000, 6000000, 1000),
                Point(500010, 6000000, 1000),
            ],
            crs=crs,
        )

        with patch.object(
            ActiveSpaceGenerator,
            "_determine_underground_pts",
            return_value=(source_pts, source_pts.iloc[0:0]),
        ):
            with pytest.raises(RuntimeError, match="no predictions"):
                gen._run_propagation_model(
                    "test_job",
                    source_pts,
                    "/fake/omni.omni",
                    altitude_m=1000,
                )

    def test_run_propagation_model_marks_partial_failures_inaudible(self, tmp_path) -> None:
        crs = "EPSG:32606"
        study_area = gpd.GeoDataFrame(
            geometry=[Point(500000, 6000000).buffer(5000)],
            crs=crs,
        )
        gen = ActiveSpaceGenerator(
            NMSIM=None,
            study_area=study_area,
            root_dir=str(tmp_path),
            ambience=40.0,
            propagation_model=_PartialFailurePropagationModel(),
        )
        source_pts = gpd.GeoDataFrame(
            geometry=[
                Point(500000, 6000000, 1000),
                Point(500010, 6000000, 1000),
            ],
            crs=crs,
        )

        with patch.object(
            ActiveSpaceGenerator,
            "_determine_underground_pts",
            return_value=(source_pts, source_pts.iloc[0:0]),
        ):
            audibility_pts = gen._run_propagation_model(
                "test_job",
                source_pts,
                "/fake/omni.omni",
                altitude_m=1000,
            )

        assert len(audibility_pts) == 2
        by_x = audibility_pts.assign(x=audibility_pts.geometry.x).sort_values("x")
        assert by_x["audible"].tolist() == [1, 0]
        crs = "EPSG:32606"
        study_area = gpd.GeoDataFrame(
            geometry=[Point(500000, 6000000).buffer(5000)],
            crs=crs,
        )
        gen = ActiveSpaceGenerator(
            NMSIM=None,
            study_area=study_area,
            root_dir=str(tmp_path),
            ambience=pd.Series({"1000": 40.0, "12500": 40.0}),
            propagation_model=_EmptyPredictionPropagationModel(),
        )

        audible_pts = gen._find_audible_points(pd.DataFrame(), crs)

        assert len(audible_pts) == 0
        assert list(audible_pts.columns) == ["audible", "geometry"]
        assert audible_pts.crs.to_string().lower() == crs.lower()

    def test_source_pts_missing_predictions(self) -> None:
        crs = "EPSG:32606"
        pts = gpd.GeoDataFrame(
            geometry=[Point(500000, 6000000, 1000), Point(500010, 6000000, 1000)],
            crs=crs,
        )
        preds = pd.DataFrame({
            "Xpos": [500000.0],
            "Ypos": [6000000.0],
            "Zpos": [1000.0],
            "A": [50.0],
        })
        missing = ActiveSpaceGenerator._source_pts_missing_predictions(pts, preds)
        assert len(missing) == 1
        assert missing.geometry.x.iloc[0] == 500010.0
