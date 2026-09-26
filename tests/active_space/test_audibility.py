"""Tests for audibility policy helpers."""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
import pandas.testing as pdt
from geopandas import testing as gdt

from nps_active_space.active_space.audibility import (
    audible_points_gdf,
    human_hearing_threshold,
    spectrum_is_audible,
)

BAND_COLUMNS = list(human_hearing_threshold.index)


def _spectral_pred_df(rows: list[dict[str, float]]) -> pd.DataFrame:
    base = {
        "Xpos": 0.0,
        "Ypos": 0.0,
        "Zpos": 0.0,
        "A": 0.0,
        **{band: 0.0 for band in BAND_COLUMNS},
    }
    return pd.DataFrame([{**base, **row} for row in rows])


def _spectral_ambience(level: float) -> pd.Series:
    return pd.Series({band: level for band in BAND_COLUMNS})


class TestSpectrumIsAudible:
    def test_broadband_above_ambience_is_audible(self) -> None:
        pred_df = pd.DataFrame({"A": [50.0, 40.0]})
        result = spectrum_is_audible(pred_df, 45.0)
        pdt.assert_series_equal(result, pd.Series([True, False], name=result.name))

    def test_spectral_band_above_threshold_is_audible(self) -> None:
        pred_df = _spectral_pred_df([{"1000": 45.0}])
        ambience = _spectral_ambience(40.0)
        assert spectrum_is_audible(pred_df, ambience).iloc[0]

    def test_spectral_all_bands_below_threshold_is_inaudible(self) -> None:
        pred_df = _spectral_pred_df([{"1000": 35.0, "4000": 30.0}])
        ambience = _spectral_ambience(40.0)
        assert not spectrum_is_audible(pred_df, ambience).iloc[0]

    def test_hearing_threshold_wins_over_quiet_ambience(self) -> None:
        pred_df = _spectral_pred_df([{"20": 70.0}])
        ambience = _spectral_ambience(10.0)
        assert not spectrum_is_audible(pred_df, ambience).iloc[0]

        pred_df = _spectral_pred_df([{"20": 80.0}])
        assert spectrum_is_audible(pred_df, ambience).iloc[0]
        assert human_hearing_threshold["20"] == 78.5


class TestAudiblePointsGdf:
    def test_empty_pred_df_returns_empty_gdf(self) -> None:
        crs = "EPSG:32606"
        result = audible_points_gdf(pd.DataFrame(), _spectral_ambience(40.0), crs)
        expected = gpd.GeoDataFrame(
            columns=["audible", "geometry"],
            geometry="geometry",
            crs=crs,
        )
        gdt.assert_geodataframe_equal(result, expected)

    def test_broadband_audible_points_encoded_as_ints(self) -> None:
        crs = "EPSG:32606"
        pred_df = pd.DataFrame({
            "Xpos": [1.0, 2.0],
            "Ypos": [3.0, 4.0],
            "Zpos": [5.0, 6.0],
            "A": [50.0, 40.0],
        })
        result = audible_points_gdf(pred_df, 45.0, crs)
        assert result["audible"].tolist() == [1, 0]
        assert result.crs.to_string().lower() == crs.lower()
        assert list(result.columns) == ["audible", "geometry"]

    def test_spectral_audible_point_geometry_from_positions(self) -> None:
        crs = "EPSG:32606"
        pred_df = _spectral_pred_df([
            {"Xpos": 100.0, "Ypos": 200.0, "Zpos": 300.0, "1000": 45.0},
            {"Xpos": 110.0, "Ypos": 210.0, "Zpos": 310.0, "1000": 35.0},
        ])
        result = audible_points_gdf(pred_df, _spectral_ambience(40.0), crs)
        assert result["audible"].tolist() == [1, 0]
        assert result.geometry.x.tolist() == [100.0, 110.0]
        assert result.geometry.y.tolist() == [200.0, 210.0]
        assert result.geometry.z.tolist() == [300.0, 310.0]
