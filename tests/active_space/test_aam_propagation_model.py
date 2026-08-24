"""Tests for AAM propagation model output mapping."""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from aam_translator import read_poi, read_run_log

from nps_active_space.active_space.aam_propagation_model import (
    aam_source_id_from_omni,
    poi_history_to_predictions_df,
)
from nps_active_space.active_space.propagation_model import NMSIM_BAND_COLUMNS

FIXTURES = Path(__file__).parent / "fixtures" / "tier4_two_point"


class TestAamSourceMapping:
    def test_omni_o_plus_200_maps_to_flato200(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+200.avg") == "FLATO200"

    def test_nmsim_omni_stem_maps_to_flato200(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+000.src") == "FLATO200"
        assert aam_source_id_from_omni("/data/tuning/O_+005.src") == "FLATO200"


class TestPoiHistoryMapping:
    @pytest.fixture
    def case_meta(self) -> dict:
        return json.loads((FIXTURES / "case_meta.json").read_text())

    @pytest.fixture
    def poi_history(self):
        histories = read_poi(FIXTURES / "scenario.POI")
        assert len(histories) == 1
        return histories[0]

    @pytest.fixture
    def source_pts(self, case_meta: dict) -> gpd.GeoDataFrame:
        rows = case_meta["source_points_utm"]
        crs = "EPSG:32606"
        geoms = [Point(r["x"], r["y"], r["z"]) for r in rows]
        return gpd.GeoDataFrame({"label": [r["label"] for r in rows]}, geometry=geoms, crs=crs)

    def test_poi_maps_to_nmsim_columns(self, poi_history, source_pts: gpd.GeoDataFrame) -> None:
        frame = poi_history_to_predictions_df(poi_history, source_pts)
        expected_cols = {"Xpos", "Ypos", "Zpos", "A", *NMSIM_BAND_COLUMNS}
        assert expected_cols == set(frame.columns)
        assert len(frame) == 2
        assert frame["A"].notna().all()

    def test_run_log_ok(self) -> None:
        log = read_run_log(FIXTURES / "scenario.txt")
        assert log.ok
        assert not log.read_error
