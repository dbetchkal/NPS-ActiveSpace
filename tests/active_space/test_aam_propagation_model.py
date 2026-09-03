"""Tests for AAM propagation model output mapping."""

from __future__ import annotations

import json
import math
import os
import pickle
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

pytest.importorskip("aam_translator")

from aam_translator import read_poi, read_run_log
from aam_translator.write_inp import TrackPoint

from nps_active_space.active_space.aam_output import poi_history_to_predictions_df
from nps_active_space.active_space.aam_source import aam_source_id_from_omni, site_ncfiles_dir
from nps_active_space.active_space.aam_source import AAM_TEMPLATE_NC_FILENAME
from nps_active_space.active_space.aam_propagation_model import (
    AAM_PREDICTIONS_SUBDIR,
    SINGLE_TRACK_PAD_M,
    AamPropagationModel,
    _aam_subprocess_env,
    _order_source_pts_for_track,
    _pad_single_point_track,
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

    def test_order_source_pts_sorts_by_xy(self) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2]},
            geometry=[Point(2, 0, 100), Point(0, 0, 100), Point(1, 0, 100)],
            crs="EPSG:32606",
        )
        ordered = _order_source_pts_for_track(pts)
        assert ordered["id"].tolist() == [1, 2, 0]

    def test_order_source_pts_snakes_grid_columns(self) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2, 3]},
            geometry=[
                Point(0, 0, 100),
                Point(0, 1, 100),
                Point(1, 0, 100),
                Point(1, 1, 100),
            ],
            crs="EPSG:32606",
        )
        ordered = _order_source_pts_for_track(pts)
        assert ordered["id"].tolist() == [0, 1, 3, 2]


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


class TestAamPredictSkipOnFailure:
    def _dummy_site(self):
        return SimpleNamespace(terrain=None)

    def _passthrough_filter(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def passthrough(self, site, source_pts, job_name=""):
            return source_pts, source_pts.iloc[0:0]

        monkeypatch.setattr(AamPropagationModel, "filter_below_terrain", passthrough)
        monkeypatch.setattr(
            "nps_active_space.active_space.aam_propagation_model.split_safe_aam_track_runs",
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
        with pytest.raises(RuntimeError, match="no predictions"):
            model.predict(site, source_pts, "O_+000.src", 1000, "skip_job")

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

    def test_all_batches_fail_raises(
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
        with pytest.raises(RuntimeError, match="no predictions"):
            model.predict(site, source_pts, "O_+000.src", 1000, "solo_job")

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
            "nps_active_space.active_space.aam_propagation_model.split_safe_aam_track_runs",
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


class TestAamSubprocessEnv:
    def test_sets_site_ncfiles_for_shim(self, tmp_path: Path) -> None:
        site_root = tmp_path / "site"
        nc_dir = site_ncfiles_dir(site_root)
        nc_dir.mkdir(parents=True)
        shim = tmp_path / "aam"
        shim.write_text("#!/bin/sh\n")
        env = _aam_subprocess_env(shim, nc_dir)
        expected = str(nc_dir.resolve()) + os.sep
        assert env["ROTOR_NOISE"] == expected
        assert env["AAM_NC"] == str(nc_dir.resolve())

    def test_exe_sets_noise_paths_from_site_ncfiles(self, tmp_path: Path) -> None:
        site_root = tmp_path / "site"
        nc_dir = site_ncfiles_dir(site_root)
        nc_dir.mkdir(parents=True)
        exe = tmp_path / "AAM_3.0.0.exe"
        exe.write_bytes(b"")
        env = _aam_subprocess_env(exe, nc_dir)
        expected = str(nc_dir.resolve()) + os.sep
        assert env["ROTOR_NOISE"] == expected
        assert env["FWING_NOISE"] == expected
        assert env["QUARRY_NOISE"] == expected
        assert env["AAM_NC"] == str(nc_dir.resolve())

    def test_template_resolution_prefers_parent_when_bin_stub_lacks_template(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        from nps_active_space.active_space.aam_propagation_model import (
            _resolve_aam_template_ncfiles_dir,
        )

        monkeypatch.delenv("AAM_NC", raising=False)
        bin_dir = tmp_path / "Bin"
        bin_dir.mkdir()
        exe = bin_dir / "AAM_3.0.0.exe"
        exe.write_bytes(b"")
        stub = bin_dir / "NCfiles"
        stub.mkdir()
        nc = tmp_path / "NCfiles"
        nc.mkdir()
        (nc / AAM_TEMPLATE_NC_FILENAME).write_bytes(b"")
        resolved = _resolve_aam_template_ncfiles_dir(exe)
        assert resolved == nc

    def test_aam_home_template_resolution(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        from nps_active_space.active_space.aam_propagation_model import (
            _resolve_aam_template_ncfiles_dir,
        )

        monkeypatch.delenv("AAM_NC", raising=False)
        aam_home = tmp_path / "opt" / "aam"
        nc = aam_home / "NCfiles"
        nc.mkdir(parents=True)
        (nc / AAM_TEMPLATE_NC_FILENAME).write_bytes(b"")
        shim = tmp_path / "usr" / "local" / "bin" / "aam"
        shim.parent.mkdir(parents=True)
        shim.write_text("#!/bin/sh\n")
        monkeypatch.setenv("AAM_HOME", str(aam_home))
        resolved = _resolve_aam_template_ncfiles_dir(shim)
        assert resolved == nc

    def test_aam_nc_override_for_template(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        from nps_active_space.active_space.aam_propagation_model import (
            _resolve_aam_template_ncfiles_dir,
        )

        exe = tmp_path / "Bin" / "AAM_3.0.0.exe"
        exe.parent.mkdir(parents=True)
        exe.write_bytes(b"")
        nc = tmp_path / "custom" / "NCfiles"
        nc.mkdir(parents=True)
        (nc / AAM_TEMPLATE_NC_FILENAME).write_bytes(b"")
        monkeypatch.setenv("AAM_NC", str(nc))
        resolved = _resolve_aam_template_ncfiles_dir(exe)
        assert resolved == nc


class TestPadSinglePointTrack:
    def test_leaves_multi_point_track_unchanged(self) -> None:
        track = [
            TrackPoint(lon=-148.87, lat=63.66, alt_m=1500.0),
            TrackPoint(lon=-148.86, lat=63.66, alt_m=1500.0),
        ]
        assert _pad_single_point_track(track) == track

    def test_pads_one_vertex_about_one_meter_east(self) -> None:
        point = TrackPoint(lon=-148.87, lat=63.66, alt_m=1500.0)
        padded = _pad_single_point_track([point])
        assert len(padded) == 2
        assert padded[0] == point
        assert padded[1].lat == point.lat
        assert padded[1].alt_m == point.alt_m
        meters_per_deg_lon = 111_320.0 * abs(math.cos(math.radians(point.lat)))
        east_m = (padded[1].lon - padded[0].lon) * meters_per_deg_lon
        assert east_m == pytest.approx(SINGLE_TRACK_PAD_M, rel=1e-4)
