"""Tests for PropagationModel wiring in ActiveSpaceGenerator."""

from __future__ import annotations

from typing import Callable
from unittest.mock import patch

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
from nps_active_space.propagation_model.protocol import THIRD_OCTAVE_BANDS

_WIRING_CRS = "EPSG:32606"


def _study_area_utm() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        geometry=[Point(500000, 6000000).buffer(5000)],
        crs=_WIRING_CRS,
    )


def _two_source_pts() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        geometry=[
            Point(500000, 6000000, 1000),
            Point(500010, 6000000, 1000),
        ],
        crs=_WIRING_CRS,
    )


class _ConfigurableStubPropagationModel:
    max_points_per_run = 7
    predictions_subdir = "Output_Data/stub/predictions"

    def __init__(
        self,
        predict_fn: Callable[..., pd.DataFrame] | None = None,
    ) -> None:
        self._predict_fn = predict_fn or self._predict_all

    def prepare_site(self, dem_src, study_area, mic, *, project_dem=True, suffix=""):
        return {"dem_file": dem_src}

    def _predict_all(self, site, source_pts, omni_source, altitude_m, job_name, heading=None):
        return pd.DataFrame({
            "Xpos": source_pts.geometry.x.values,
            "Ypos": source_pts.geometry.y.values,
            "Zpos": source_pts.geometry.z.values,
            "A": [50.0] * len(source_pts),
            **{col: [40.0] * len(source_pts) for col in THIRD_OCTAVE_BANDS},
        })

    def predict(self, site, source_pts, omni_source, altitude_m, job_name, heading=None):
        return self._predict_fn(site, source_pts, omni_source, altitude_m, job_name, heading)

    def filter_below_terrain(self, site, source_pts, *, job_name=""):
        return source_pts, source_pts.iloc[0:0]


class TestPropagationModelWiring:
    def test_generator_accepts_custom_propagation_model(self, tmp_path) -> None:
        study_area = gpd.GeoDataFrame(
            geometry=[Point(0, 0).buffer(1)],
            crs="EPSG:4326",
        )
        stub = _ConfigurableStubPropagationModel()
        gen = ActiveSpaceGenerator(
            study_area=study_area,
            root_dir=str(tmp_path),
            ambience=30.0,
            propagation_model=stub,
        )
        assert gen.propagation_model.max_points_per_run == 7

    def test_preprocess_uses_model_batch_cap(self, tmp_path) -> None:
        stub = _ConfigurableStubPropagationModel()
        gen = ActiveSpaceGenerator(
            study_area=_study_area_utm(),
            root_dir=str(tmp_path),
            ambience=30.0,
            propagation_model=stub,
        )
        pts = gpd.GeoDataFrame(
            geometry=[Point(500000 + i * 10, 6000000) for i in range(20)],
            crs=_WIRING_CRS,
        )
        region = Point(500000, 6000000).buffer(10000)
        tested = gpd.GeoDataFrame(columns=["audible", "geometry"], geometry="geometry", crs=_WIRING_CRS)
        filtered = gen._preprocess_source_points(
            pts, region, tested, max_pts=gen.propagation_model.max_points_per_run,
        )
        assert len(filtered) == 7


class TestMissingPredictionHandling:
    def test_run_propagation_model_marks_total_predict_failure_inaudible(
        self, tmp_path,
    ) -> None:
        gen = ActiveSpaceGenerator(
            study_area=_study_area_utm(),
            root_dir=str(tmp_path),
            ambience=pd.Series({"1000": 40.0, "12500": 40.0}),
            propagation_model=_ConfigurableStubPropagationModel(
                predict_fn=lambda *args, **kwargs: pd.DataFrame(),
            ),
        )
        source_pts = _two_source_pts()

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
        assert audibility_pts["audible"].tolist() == [0, 0]

    def test_run_propagation_model_marks_partial_failures_inaudible(self, tmp_path) -> None:
        def predict_first_only(site, source_pts, omni_source, altitude_m, job_name, heading=None):
            return pd.DataFrame({
                "Xpos": [source_pts.geometry.x.iloc[0]],
                "Ypos": [source_pts.geometry.y.iloc[0]],
                "Zpos": [source_pts.geometry.z.iloc[0]],
                "A": [50.0],
                **{col: [40.0] for col in THIRD_OCTAVE_BANDS},
            })

        gen = ActiveSpaceGenerator(
            study_area=_study_area_utm(),
            root_dir=str(tmp_path),
            ambience=40.0,
            propagation_model=_ConfigurableStubPropagationModel(predict_fn=predict_first_only),
        )
        source_pts = _two_source_pts()

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
