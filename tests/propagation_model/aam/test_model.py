"""Tests for AAM propagation model output mapping and predict orchestration."""

from __future__ import annotations

import pickle
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

pytest.importorskip("aam_translator")

from aam_translator import read_poi

from nps_active_space.propagation_model.aam.model import (
    AamPropagationModel,
    resolve_aam_chunk_size,
)
from nps_active_space.propagation_model.aam.output import poi_history_to_predictions_df
from nps_active_space.propagation_model.nmsim.model import NmsimPropagationModel
from nps_active_space.propagation_model.protocol import (
    DEFAULT_MAX_POINTS_PER_PREDICT,
    THIRD_OCTAVE_BANDS,
)
from nps_active_space.utils.paths import AAM_PREDICTIONS_SUBDIR

from nps_active_space.utils.paths import AAM_PREDICTIONS_SUBDIR

TWO_POINT_RIDGE_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "two_point_ridge"
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
        **{col: [level_db - 10.0] * len(source_pts) for col in THIRD_OCTAVE_BANDS},
    })


class TestAamPredictionsLayout:
    def test_predictions_subdir(self) -> None:
        assert AamPropagationModel("/tmp/site").predictions_subdir == AAM_PREDICTIONS_SUBDIR


class TestAamBatching:
    def test_resolve_aam_chunk_size_default(self) -> None:
        assert resolve_aam_chunk_size() == 400

    def test_resolve_aam_chunk_size_env_override(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.setenv("AAM_CHUNK_SIZE", "123")
        assert resolve_aam_chunk_size() == 123

    def test_max_points_per_run_matches_nmsim(self) -> None:
        assert AamPropagationModel.max_points_per_run == DEFAULT_MAX_POINTS_PER_PREDICT
        assert (
            AamPropagationModel.max_points_per_run
            == NmsimPropagationModel.max_points_per_run
            == 4000
        )
        assert AamPropagationModel.max_points_per_run != resolve_aam_chunk_size()


class TestPoiHistoryMapping:
    @pytest.fixture
    def poi_history(self):
        histories = read_poi(TWO_POINT_RIDGE_FIXTURES / "scenario.POI")
        assert len(histories) == 1
        return histories[0]

    def test_poi_maps_to_nmsim_columns(
        self,
        poi_history,
        ridge_source_pts: gpd.GeoDataFrame,
    ) -> None:
        frame = poi_history_to_predictions_df(poi_history, ridge_source_pts)
        expected_cols = {"Xpos", "Ypos", "Zpos", "A", *THIRD_OCTAVE_BANDS}
        assert expected_cols == set(frame.columns)
        assert len(frame) == 2
        assert frame["A"].notna().all()


class TestAamMultiprocessPickle:
    def test_unpickle_reconfigures_site_log(self, tmp_path: Path) -> None:
        from nps_active_space.propagation_model.aam import run_log as aam_run_log

        site_root = tmp_path / "site"
        site_root.mkdir()
        model = AamPropagationModel(str(site_root))
        aam_run_log._configured_root = None
        aam_run_log._log_path = None

        restored = pickle.loads(pickle.dumps(model))
        assert restored.root_dir == str(site_root)
        assert restored._runs_dir == model._runs_dir
        assert aam_run_log._log_path == aam_run_log.aam_run_log_path(restored._root)

    def test_site_context_pickles_after_prepare(self, tmp_path: Path, case_meta: dict) -> None:
        from rasterio import open as rio_open
        from shapely.geometry import box

        from nps_active_space.utils.models import Microphone

        dem_path = TWO_POINT_RIDGE_FIXTURES / case_meta["dem_utm"]
        root = tmp_path / "site"
        (root / "Input_Data").mkdir(parents=True)
        rx_lon, rx_lat = case_meta["receiver_lonlat"]
        mic = Microphone(name="Receiver", lat=rx_lat, lon=rx_lon, z=4.92)
        with rio_open(dem_path) as ds:
            bounds = ds.bounds
        aoi = gpd.GeoDataFrame(
            geometry=[box(bounds.left, bounds.bottom, bounds.right, bounds.top)],
            crs=CRS,
        )

        model = AamPropagationModel(str(root), aam_shim="/usr/local/bin/aam")
        site = model.prepare_site(str(dem_path), aoi, mic, project_dem=False)
        restored_model, restored_site = pickle.loads(pickle.dumps((model, site)))
        assert restored_site.terrain.elv_path == site.terrain.elv_path
        assert restored_model._root == model._root


class TestAamPredictSkipOnFailure:
    def test_chunk_failure_continues_other_chunks(
        self,
        monkeypatch: pytest.MonkeyPatch,
        aam_predict_harness,
    ) -> None:
        model, site = aam_predict_harness
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
            if "_r001" in job_name:
                raise RuntimeError("simulated AAM abort")
            return _predictions_for(batch_pts)

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "mesh_job")

        assert len(result) == 50
        assert set(result["Xpos"]) == set(xs[:50])

    def test_non_fpa_failure_is_skipped_not_bisected(
        self,
        monkeypatch: pytest.MonkeyPatch,
        aam_predict_harness,
    ) -> None:
        model, site = aam_predict_harness
        source_pts = _make_source_pts([500000.0, 500010.0, 500020.0])

        def fake_batch(self, *args, **kwargs) -> pd.DataFrame:
            raise RuntimeError("AAM abort")

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "skip_job")
        assert result.empty

    def test_fpa_bounds_splits_batch_and_retries(
        self,
        monkeypatch: pytest.MonkeyPatch,
        aam_predict_harness,
    ) -> None:
        model, site = aam_predict_harness
        xs = [float(i) for i in range(4)]
        source_pts = _make_source_pts(xs)
        fpa_error = (
            "forrtl: severe (408): fort: (11): Subscript #2 of the array FPA "
            "has value 0 which is less than the lower bound of 1"
        )

        def fake_batch(
            self,
            site_ctx,
            batch_pts,
            omni_source,
            altitude_m,
            job_name,
            heading=None,
        ) -> pd.DataFrame:
            if len(batch_pts) >= 3:
                raise RuntimeError(fpa_error)
            return _predictions_for(batch_pts)

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1800, "mesh_job")

        assert len(result) == 4
        assert set(result["Xpos"]) == set(xs)

    def test_all_batches_fail_returns_empty(
        self,
        monkeypatch: pytest.MonkeyPatch,
        aam_predict_harness,
    ) -> None:
        model, site = aam_predict_harness
        source_pts = _make_source_pts([500000.0])

        def fake_batch(self, *args, **kwargs) -> pd.DataFrame:
            raise RuntimeError("single point below ground")

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "solo_job")
        assert result.empty

    def test_predict_issues_one_batch_per_hop_run(
        self,
        monkeypatch: pytest.MonkeyPatch,
        aam_predict_harness,
    ) -> None:
        model, site = aam_predict_harness
        source_pts = _make_source_pts([500000.0, 500010.0, 500020.0, 500030.0])
        monkeypatch.setattr(
            "nps_active_space.propagation_model.aam.model.split_safe_aam_track_runs",
            lambda terrain, pts, job_name="": [pts.iloc[:2], pts.iloc[2:]],
        )
        jobs: list[str] = []

        def fake_batch(
            self,
            site_ctx,
            batch_pts,
            omni_source,
            altitude_m,
            job_name,
            heading=None,
        ) -> pd.DataFrame:
            jobs.append(job_name)
            return _predictions_for(batch_pts)

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "split_job")

        assert jobs == ["split_job_r000", "split_job_r001"]
        assert len(result) == 4
