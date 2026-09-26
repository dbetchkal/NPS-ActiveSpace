"""Tests for prediction CSV cache read/write."""

from __future__ import annotations

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from nps_active_space.active_space.prediction_cache import (
    load_prediction_cache,
    predict_with_cache,
    prediction_cache_csv_path,
    save_prediction_cache,
    source_pts_missing_predictions,
)
from nps_active_space.utils.paths import AAM_PREDICTIONS_SUBDIR

_CACHE_X = 407202.0
_CACHE_Y = 7060771.0
_CACHE_Z = 600.0
_CACHE_CRS = "epsg:26906"


@pytest.fixture
def cache_source_pt() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        geometry=[Point(_CACHE_X, _CACHE_Y, _CACHE_Z)],
        crs=_CACHE_CRS,
    )


@pytest.fixture
def cache_prediction_row() -> pd.DataFrame:
    return pd.DataFrame({
        "Xpos": [_CACHE_X],
        "Ypos": [_CACHE_Y],
        "Zpos": [_CACHE_Z],
        "A": [45.0],
        "1000": [40.0],
        "12500": [-99.9],
    })


class TestPredictionCache:
    @pytest.mark.parametrize(
        ("csv_content", "log_substrings"),
        [
            ("", ["0 bytes", "Removing file"]),
            ("Xpos,Ypos\n407202,7060771\n", ["missing # units=dB header", "Removing file"]),
            (
                "Xpos,Ypos,A,1000\n407202,7060771,450,400\n",
                ["missing # units=dB header", "Removing file"],
            ),
            (
                "# units=dB\nXpos,Ypos\n407202,7060771\n",
                ["missing required columns", "'A'"],
            ),
        ],
    )
    def test_load_treats_invalid_cache_as_missing(
        self,
        tmp_path: Path,
        caplog,
        csv_content: str,
        log_substrings: list[str],
    ) -> None:
        caplog.set_level(logging.WARNING)
        csv_path = tmp_path / "600m_O_10deg.csv"
        csv_path.write_text(csv_content)

        source_pts = gpd.GeoDataFrame(
            geometry=[Point(_CACHE_X, _CACHE_Y, _CACHE_Z)],
            crs=_CACHE_CRS,
        )
        cached_all, cached_hits, new_pts = load_prediction_cache(
            source_pts, str(csv_path), altitude_m=600,
        )

        assert cached_all.empty
        assert cached_hits.empty
        assert len(new_pts) == 1
        assert not csv_path.exists()
        for substring in log_substrings:
            assert substring in caplog.text

    def test_save_skips_empty_dataframe(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "cache.csv"
        save_prediction_cache(pd.DataFrame(), str(csv_path))
        assert not csv_path.exists()

    def test_save_and_load_roundtrip_db(self, tmp_path: Path, cache_source_pt: gpd.GeoDataFrame) -> None:
        csv_path = tmp_path / "600m_O_0deg.csv"
        predictions = pd.DataFrame({
            "Xpos": [_CACHE_X],
            "Ypos": [_CACHE_Y],
            "Zpos": [_CACHE_Z],
            "A": [45.0],
            "12.5": [21.3],
            "1000": [40.0],
            "12500": [-99.9],
        })
        save_prediction_cache(predictions, str(csv_path))
        written = csv_path.read_text()
        assert written.startswith("# units=dB")
        assert ",45.0," in written or ",45," in written
        assert ",450," not in written

        loaded, hits, new_pts = load_prediction_cache(
            cache_source_pt, str(csv_path), altitude_m=600,
        )
        assert new_pts.empty
        assert float(hits.iloc[0]["A"]) == 45.0
        assert float(hits.iloc[0]["12.5"]) == 21.3
        assert float(loaded.iloc[0]["12500"]) == -99.9

    def test_predict_with_cache_skips_predict_on_full_hit(
        self,
        tmp_path: Path,
        cache_source_pt: gpd.GeoDataFrame,
        cache_prediction_row: pd.DataFrame,
    ) -> None:
        csv_path = tmp_path / "600m_O_0deg.csv"
        save_prediction_cache(cache_prediction_row, str(csv_path))
        calls: list[int] = []

        def predict_fn(pts):
            calls.append(len(pts))
            raise AssertionError("predict should not run on a full cache hit")

        pred_df, failed_pts = predict_with_cache(
            predict_fn, cache_source_pt, str(csv_path), altitude_m=600, job_name="job",
        )
        assert calls == []
        assert failed_pts.empty
        assert float(pred_df.iloc[0]["A"]) == 45.0

    def test_predict_with_cache_rebuilds_when_units_header_missing(
        self,
        tmp_path: Path,
        cache_source_pt: gpd.GeoDataFrame,
    ) -> None:
        csv_path = tmp_path / "600m_O_0deg.csv"
        csv_path.write_text("Xpos,Ypos,A,1000,12500\n407202,7060771,450,400,-999\n")

        def predict_fn(pts):
            assert len(pts) == 1
            return pd.DataFrame({
                "Xpos": [_CACHE_X],
                "Ypos": [_CACHE_Y],
                "Zpos": [_CACHE_Z],
                "A": [45.0],
                "1000": [40.0],
                "12500": [-99.9],
            })

        pred_df, failed_pts = predict_with_cache(
            predict_fn, cache_source_pt, str(csv_path), altitude_m=600, job_name="job",
        )
        assert failed_pts.empty
        assert float(pred_df.iloc[0]["A"]) == 45.0
        written = csv_path.read_text()
        assert written.startswith("# units=dB")
        assert ",450," not in written

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
        missing = source_pts_missing_predictions(pts, preds)
        assert len(missing) == 1
        assert missing.geometry.x.iloc[0] == 500010.0

    def test_aam_prediction_cache_path(self, tmp_path: Path) -> None:
        path = prediction_cache_csv_path(
            str(tmp_path),
            AAM_PREDICTIONS_SUBDIR,
            1000,
            "O_+000",
            0,
        )
        assert path.endswith("Output_Data/aam/predictions/1000m_O_+000_0deg.csv")
        assert (tmp_path / "Output_Data" / "aam" / "predictions").is_dir()
