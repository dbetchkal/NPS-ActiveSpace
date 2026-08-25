"""Tests for AAM propagation model output mapping."""

from __future__ import annotations

import json
import os
import pickle
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

pytest.importorskip("aam_translator")

from aam_translator import read_poi, read_run_log

from nps_active_space.active_space.aam_output import (
    aam_source_id_from_omni,
    poi_history_to_predictions_df,
)
from nps_active_space.active_space.aam_propagation_model import (
    AAM_PREDICTIONS_SUBDIR,
    AamPropagationModel,
    _aam_subprocess_env,
    _order_source_pts_for_track,
    resolve_aam_chunk_size,
)
from nps_active_space.active_space.propagation_model import (
    NMSIM_BAND_COLUMNS,
    prediction_cache_csv_path,
)

FIXTURES = Path(__file__).parent / "fixtures" / "two_point_ridge"
CRS = "EPSG:32606"


def _make_source_pts(xs: list[float], y: float = 6000000.0, z: float = 1000.0) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"id": list(range(len(xs)))},
        geometry=[Point(x, y, z) for x in xs],
        crs=CRS,
    )


def _predictions_for(source_pts: gpd.GeoDataFrame, level_db: float = 50.0) -> pd.DataFrame:
    return pd.DataFrame({
        "Xpos": source_pts.geometry.x.values,
        "Ypos": source_pts.geometry.y.values,
        "Zpos": source_pts.geometry.z.values,
        "A": [level_db] * len(source_pts),
        **{col: [level_db - 10.0] * len(source_pts) for col in NMSIM_BAND_COLUMNS},
    })


class TestAamSourceMapping:
    def test_omni_o_plus_200_maps_to_flato200(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+200.avg") == "FLATO200"

    def test_nmsim_omni_stem_maps_to_flato200(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+000.src") == "FLATO200"
        assert aam_source_id_from_omni("/data/tuning/O_+005.src") == "FLATO200"


class TestAamPredictionsLayout:
    def test_predictions_subdir(self) -> None:
        assert AamPropagationModel("/tmp/site").predictions_subdir == AAM_PREDICTIONS_SUBDIR

    def test_prediction_cache_path(self, tmp_path: Path) -> None:
        path = prediction_cache_csv_path(
            str(tmp_path),
            AAM_PREDICTIONS_SUBDIR,
            1000,
            "O_+000",
            0,
        )
        assert path.endswith("Output_Data/aam/predictions/1000m_O_+000_0deg.csv")
        assert (tmp_path / "Output_Data" / "aam" / "predictions").is_dir()


class TestAamBatching:
    def test_resolve_aam_chunk_size_default(self) -> None:
        assert resolve_aam_chunk_size() == 50

    def test_order_source_pts_sorts_by_xy(self) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2]},
            geometry=[Point(2, 0, 100), Point(0, 0, 100), Point(1, 0, 100)],
            crs="EPSG:32606",
        )
        ordered = _order_source_pts_for_track(pts)
        assert ordered["id"].tolist() == [1, 2, 0]


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


class TestAamMultiprocessPickle:
    def test_unpickle_reconfigures_site_log(self, tmp_path: Path) -> None:
        from nps_active_space.active_space import aam_run_log

        site_root = tmp_path / "site"
        site_root.mkdir()
        model = AamPropagationModel(str(site_root))
        aam_run_log._configured_root = None
        aam_run_log._log_path = None

        restored = pickle.loads(pickle.dumps(model))
        assert restored.root_dir == str(site_root)
        assert restored._runs_dir == model._runs_dir
        assert aam_run_log._log_path == aam_run_log.aam_run_log_path(restored._root)

    def test_site_context_pickles_after_prepare(self, tmp_path: Path) -> None:
        from rasterio import open as rio_open
        from shapely.geometry import box

        from nps_active_space.utils.models import Microphone

        meta = json.loads((FIXTURES / "case_meta.json").read_text())
        dem_path = FIXTURES / "parent_dem_utm.tif"
        root = tmp_path / "site"
        (root / "Input_Data").mkdir(parents=True)
        rx_lon, rx_lat = meta["receiver_lonlat"]
        mic = Microphone(name="Receiver", lat=rx_lat, lon=rx_lon, z=4.92)
        with rio_open(dem_path) as ds:
            bounds = ds.bounds
        aoi = gpd.GeoDataFrame(
            geometry=[box(bounds.left, bounds.bottom, bounds.right, bounds.top)],
            crs="EPSG:32606",
        )

        model = AamPropagationModel(str(root), aam_shim="/usr/local/bin/aam")
        site = model.prepare_site(str(dem_path), aoi, mic, project_dem=False)
        restored_model, restored_site = pickle.loads(pickle.dumps((model, site)))
        assert restored_site.terrain.elv_path == site.terrain.elv_path
        assert restored_model._root == model._root


class TestAamResilientPredict:
    def test_chunk_failure_continues_other_chunks(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        model = AamPropagationModel(str(tmp_path))
        site = object()
        xs = [float(i) for i in range(75)]
        source_pts = _make_source_pts(xs)
        monkeypatch.setenv("AAM_CHUNK_SIZE", "50")

        def fake_batch(
            self,
            site_ctx,
            batch_pts,
            omni_source,
            altitude_m,
            job_name,
            heading=None,
        ) -> pd.DataFrame:
            if "_c001" in job_name:
                raise RuntimeError("simulated AAM abort")
            return _predictions_for(batch_pts)

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "mesh_job")

        assert len(result) == 50
        assert set(result["Xpos"]) == set(xs[:50])

    def test_binary_split_isolates_bad_point(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        model = AamPropagationModel(str(tmp_path))
        site = object()
        bad_x = 500020.0
        xs = [500000.0, 500010.0, bad_x, 500030.0]
        source_pts = _make_source_pts(xs)

        def fake_batch(
            self,
            site_ctx,
            batch_pts,
            omni_source,
            altitude_m,
            job_name,
            heading=None,
        ) -> pd.DataFrame:
            batch_xs = set(batch_pts.geometry.x)
            if len(batch_pts) > 1 and bad_x in batch_xs:
                raise RuntimeError("bad point poisons batch")
            if bad_x in batch_xs:
                raise RuntimeError("bad point alone")
            return _predictions_for(batch_pts)

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "split_job")

        assert len(result) == 3
        assert bad_x not in set(result["Xpos"])
        assert set(result["Xpos"]) == {500000.0, 500010.0, 500030.0}

    def test_single_point_failure_returns_empty(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        model = AamPropagationModel(str(tmp_path))
        site = object()
        source_pts = _make_source_pts([500000.0])

        def fake_batch(self, *args, **kwargs) -> pd.DataFrame:
            raise RuntimeError("single point below ground")

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "solo_job")

        assert result.empty


class TestAamSubprocessEnv:
    def test_shim_does_not_inject_ncfiles(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.delenv("ROTOR_NOISE", raising=False)
        shim = tmp_path / "aam"
        shim.write_text("#!/bin/sh\n")
        env = _aam_subprocess_env(shim)
        assert "ROTOR_NOISE" not in env

    def test_exe_sets_noise_paths_from_sibling_ncfiles(self, tmp_path: Path) -> None:
        exe = tmp_path / "AAM_3.0.0.exe"
        exe.write_bytes(b"")
        nc = tmp_path / "NCfiles"
        nc.mkdir()
        env = _aam_subprocess_env(exe)
        expected = str(nc.resolve()) + os.sep
        assert env["ROTOR_NOISE"] == expected
        assert env["FWING_NOISE"] == expected
        assert env["QUARRY_NOISE"] == expected
