"""Tests for NMSIM prediction CSV cache read/write."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator


class TestNmsimPredictionCache:
    def test_load_treats_empty_csv_as_missing(self, tmp_path: Path):
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
