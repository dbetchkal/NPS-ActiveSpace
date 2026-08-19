"""Tests for NMSIM prediction CSV cache read/write."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator


class TestNmsimPredictionCache:
    def test_load_treats_empty_csv_as_missing(self, tmp_path: Path, caplog):
        import logging

        caplog.set_level(logging.WARNING)
        csv_path = tmp_path / "600m_O_10deg.csv"
        csv_path.write_text("")

        source_pts = gpd.GeoDataFrame(
            geometry=[Point(100.0, 200.0, 600.0)],
            crs="epsg:26906",
        )
        nmsim_df_all, nmsim_df, new_pts = ActiveSpaceGenerator.load_prev_nmsim_predictions(
            source_pts, str(csv_path), altitude_m=600,
        )

        assert nmsim_df_all.empty
        assert nmsim_df.empty
        assert len(new_pts) == 1
        assert not csv_path.exists()
        assert "0 bytes" in caplog.text
        assert "Removing file" in caplog.text

    def test_load_warns_on_missing_required_columns(self, tmp_path: Path, caplog):
        import logging

        caplog.set_level(logging.WARNING)
        csv_path = tmp_path / "600m_O_0deg.csv"
        csv_path.write_text("Xpos,Ypos\n407202,7060771\n")

        source_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0)],
            crs="epsg:26906",
        )
        ActiveSpaceGenerator.load_prev_nmsim_predictions(
            source_pts, str(csv_path), altitude_m=600,
        )

        assert "missing required columns" in caplog.text
        assert "'A'" in caplog.text
        assert not csv_path.exists()

    def test_save_skips_empty_dataframe(self, tmp_path: Path):
        csv_path = tmp_path / "cache.csv"
        ActiveSpaceGenerator.save_nmsim_predictions(pd.DataFrame(), str(csv_path))
        assert not csv_path.exists()

    def test_roundtrip_predictions(self, tmp_path: Path):
        csv_path = tmp_path / "600m_O_0deg.csv"
        csv_path.write_text(
            "Xpos,Ypos,A,1000\n"
            "407202,7060771,450,400\n"
            "408000,7061000,300,250\n"
        )
        source_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0), Point(408000.0, 7061000.0, 600.0)],
            crs="epsg:26906",
        )

        nmsim_df_all, nmsim_df, new_pts = ActiveSpaceGenerator.load_prev_nmsim_predictions(
            source_pts, str(csv_path), altitude_m=600,
        )
        assert len(nmsim_df_all) == 2
        assert len(nmsim_df) == 2
        assert new_pts.empty
        assert nmsim_df_all.loc[0, "A"] == 45.0

    def test_load_decodes_fractional_band_centibels(self, tmp_path: Path):
        csv_path = tmp_path / "600m_O_0deg.csv"
        csv_path.write_text(
            "Xpos,Ypos,A,12.5,1000\n"
            "407202,7060771,450,213,400\n"
        )
        source_pts = gpd.GeoDataFrame(
            geometry=[Point(407202.0, 7060771.0, 600.0)],
            crs="epsg:26906",
        )

        _, nmsim_df, new_pts = ActiveSpaceGenerator.load_prev_nmsim_predictions(
            source_pts, str(csv_path), altitude_m=600,
        )
        assert new_pts.empty
        assert nmsim_df.loc[0, "12.5"] == 21.3
        assert nmsim_df.loc[0, "1000"] == 40.0
