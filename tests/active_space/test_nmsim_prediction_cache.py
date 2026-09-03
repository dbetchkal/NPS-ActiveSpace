"""Tests for NMSIM prediction CSV cache read/write."""

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
    save_prediction_cache,
)


class TestNmsimPredictionCache:
    @pytest.mark.parametrize(
        ("csv_content", "log_substrings"),
        [
            ("", ["0 bytes", "Removing file"]),
            ("Xpos,Ypos\n407202,7060771\n", ["missing required columns", "'A'"]),
        ],
    )
    def test_load_treats_invalid_cache_as_missing(
        self,
        tmp_path: Path,
        caplog,
        csv_content: str,
        log_substrings: list[str],
    ):
        caplog.set_level(logging.WARNING)
        csv_path = tmp_path / "600m_O_10deg.csv"
        csv_path.write_text(csv_content)

        source_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0)],
            crs="epsg:26906",
        )
        nmsim_df_all, nmsim_df, new_pts = load_prediction_cache(
            source_pts, str(csv_path), altitude_m=600,
        )

        assert nmsim_df_all.empty
        assert nmsim_df.empty
        assert len(new_pts) == 1
        assert not csv_path.exists()
        for substring in log_substrings:
            assert substring in caplog.text

    def test_save_skips_empty_dataframe(self, tmp_path: Path):
        csv_path = tmp_path / "cache.csv"
        save_prediction_cache(pd.DataFrame(), str(csv_path))
        assert not csv_path.exists()

    def test_load_legacy_centibel_csv_converts_to_db(self, tmp_path: Path):
        roundtrip_csv = tmp_path / "600m_O_0deg.csv"
        roundtrip_csv.write_text(
            "Xpos,Ypos,A,1000\n"
            "407202,7060771,450,400\n"
            "408000,7061000,300,250\n"
        )
        roundtrip_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0), Point(408000.0, 7061000.0, 600.0)],
            crs="epsg:26906",
        )
        nmsim_df_all, nmsim_df, new_pts = load_prediction_cache(
            roundtrip_pts, str(roundtrip_csv), altitude_m=600,
        )
        assert len(nmsim_df_all) == 2
        assert len(nmsim_df) == 2
        assert new_pts.empty
        assert nmsim_df_all.loc[0, "A"] == 45.0

        fractional_csv = tmp_path / "600m_O_10deg.csv"
        fractional_csv.write_text(
            "Xpos,Ypos,A,12.5,1000\n"
            "407202,7060771,450,213,400\n"
        )
        fractional_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0)],
            crs="epsg:26906",
        )
        _, nmsim_df, new_pts = load_prediction_cache(
            fractional_pts, str(fractional_csv), altitude_m=600,
        )
        assert new_pts.empty
        assert nmsim_df.loc[0, "12.5"] == 21.3
        assert nmsim_df.loc[0, "1000"] == 40.0

    def test_save_and_load_roundtrip_db(self, tmp_path: Path):
        csv_path = tmp_path / "600m_O_0deg.csv"
        predictions = pd.DataFrame({
            "Xpos": [407202.0],
            "Ypos": [7060771.0],
            "Zpos": [600.0],
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

        source_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0)],
            crs="epsg:26906",
        )
        loaded, hits, new_pts = load_prediction_cache(
            source_pts, str(csv_path), altitude_m=600,
        )
        assert new_pts.empty
        assert float(hits.iloc[0]["A"]) == 45.0
        assert float(hits.iloc[0]["12.5"]) == 21.3
        assert float(loaded.iloc[0]["12500"]) == -99.9

    def test_predict_with_cache_skips_predict_on_full_hit(self, tmp_path: Path):
        csv_path = tmp_path / "600m_O_0deg.csv"
        csv_path.write_text("Xpos,Ypos,A,1000\n407202,7060771,450,400\n")
        source_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0)],
            crs="epsg:26906",
        )
        calls: list[int] = []

        def predict_fn(pts):
            calls.append(len(pts))
            raise AssertionError("predict should not run on a full cache hit")

        pred_df, failed_pts = predict_with_cache(
            predict_fn, source_pts, str(csv_path), altitude_m=600, job_name="job",
        )
        assert calls == []
        assert failed_pts.empty
        assert float(pred_df.iloc[0]["A"]) == 45.0
