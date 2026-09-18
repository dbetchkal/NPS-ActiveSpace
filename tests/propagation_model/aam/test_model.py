"""Tests for AAM propagation model output mapping."""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

pytest.importorskip("aam_translator")

from aam_translator import read_poi, read_run_log

from nps_active_space.active_space.prediction_cache import prediction_cache_csv_path
from nps_active_space.propagation_model.aam.model import (
    AamPropagationModel,
    resolve_aam_chunk_size,
)
from nps_active_space.propagation_model.nmsim.model import NmsimPropagationModel
from nps_active_space.propagation_model.aam.output import poi_history_to_predictions_df
from nps_active_space.propagation_model.aam.source import aam_source_id_from_omni
from nps_active_space.propagation_model.protocol import (
    DEFAULT_MAX_POINTS_PER_PREDICT,
    THIRD_OCTAVE_BANDS,
)
from nps_active_space.utils.paths import AAM_PREDICTIONS_SUBDIR

FIXTURES = Path(__file__).resolve().parents[2] / "active_space" / "fixtures" / "two_point_ridge"
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


class TestAamSourceMapping:
    def test_omni_o_plus_200_maps_to_omni_200(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+200.avg") == "OMNI_200"

    def test_nmsim_omni_stem_maps_to_omni_tokens(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+000.src") == "OMNI_000"
        assert aam_source_id_from_omni("/data/tuning/O_+005.src") == "OMNI_005"


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
        assert resolve_aam_chunk_size() == 400

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
        expected_cols = {"Xpos", "Ypos", "Zpos", "A", *THIRD_OCTAVE_BANDS}
        assert expected_cols == set(frame.columns)
        assert len(frame) == 2
        assert frame["A"].notna().all()

    def test_run_log_ok(self) -> None:
        log = read_run_log(FIXTURES / "scenario.txt")
        assert log.ok
        assert not log.read_error


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


class TestAamPredictSkipOnFailure:
    def _dummy_site(self):
        return SimpleNamespace(terrain=None)

    def _passthrough_filter(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def passthrough(self, site, source_pts, job_name=""):
            return source_pts, source_pts.iloc[0:0]

        monkeypatch.setattr(AamPropagationModel, "filter_below_terrain", passthrough)
        monkeypatch.setattr(
            "nps_active_space.propagation_model.aam.model.split_safe_aam_track_runs",
            lambda terrain, pts, job_name="": [pts],
        )

    def test_chunk_failure_continues_other_chunks(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        self._passthrough_filter(monkeypatch)
        model = AamPropagationModel(str(tmp_path))
        site = self._dummy_site()
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
        tmp_path: Path,
    ) -> None:
        self._passthrough_filter(monkeypatch)
        model = AamPropagationModel(str(tmp_path))
        site = self._dummy_site()
        source_pts = _make_source_pts([500000.0, 500010.0, 500020.0])

        def fake_batch(self, *args, **kwargs) -> pd.DataFrame:
            raise RuntimeError("AAM abort")

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "skip_job")
        assert result.empty

    def test_fpa_bounds_splits_batch_and_retries(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        self._passthrough_filter(monkeypatch)
        model = AamPropagationModel(str(tmp_path))
        site = self._dummy_site()
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
        tmp_path: Path,
    ) -> None:
        self._passthrough_filter(monkeypatch)
        model = AamPropagationModel(str(tmp_path))
        site = self._dummy_site()
        source_pts = _make_source_pts([500000.0])

        def fake_batch(self, *args, **kwargs) -> pd.DataFrame:
            raise RuntimeError("single point below ground")

        monkeypatch.setattr(AamPropagationModel, "_predict_batch", fake_batch)
        result = model.predict(site, source_pts, "O_+000.src", 1000, "solo_job")
        assert result.empty

    def test_predict_issues_one_batch_per_hop_run(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        source_pts = _make_source_pts([500000.0, 500010.0, 500020.0, 500030.0])

        def passthrough(self, site, pts, job_name=""):
            return pts, pts.iloc[0:0]

        monkeypatch.setattr(AamPropagationModel, "filter_below_terrain", passthrough)
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
        model = AamPropagationModel(str(tmp_path))
        result = model.predict(self._dummy_site(), source_pts, "O_+000.src", 1000, "split_job")

        assert jobs == ["split_job_r000", "split_job_r001"]
        assert len(result) == 4
