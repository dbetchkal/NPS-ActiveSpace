"""Audibility policy: compare predicted levels to ambience and hearing threshold."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd

# dB threshold of human hearing at each 1/3 octave band, from ISO 226:2003
human_hearing_threshold = pd.Series({
    "20": 78.5, "25": 68.7, "31.5": 59.5, "40": 51.1, "50": 44.0, "63": 37.5, "80": 31.5,
    "100": 26.5, "125": 22.1, "160": 17.9, "200": 14.4, "250": 11.4, "315": 8.6, "400": 6.2,
    "500": 4.4, "630": 3.0, "800": 2.2, "1000": 2.4, "1250": 3.5, "1600": 1.7, "2000": -1.3,
    "2500": -4.2, "3150": -6.0, "4000": -5.4, "5000": -1.5, "6300": 6.0, "8000": 12.6,
    "10000": 13.9, "12500": 12.3
})


def spectrum_is_audible(
    pred_df: pd.DataFrame,
    ambience: float | int | pd.Series,
) -> pd.Series:
    """Boolean Series: True where predicted levels exceed ambience (and hearing threshold for spectral).

    NaN comparisons are False, so a missing AAM 12.5 kHz band never marks a point audible.
    """
    if isinstance(ambience, (float, int)):
        return pred_df.loc[:, "A"] > ambience

    threshold = np.maximum(ambience.loc["20":"12500"], human_hearing_threshold)
    return (pred_df.loc[:, "20":"12500"] > threshold.values).any(axis=1)


def audible_points_gdf(
    pred_df: pd.DataFrame,
    ambience: float | int | pd.Series,
    crs: str,
) -> gpd.GeoDataFrame:
    """If pred_df.empty, return empty GDF with columns ["audible", "geometry"] and the given crs.

    Otherwise classify with spectrum_is_audible and return GDF with audible as 0/1 ints and
    geometry from Xpos/Ypos/Zpos via gpd.points_from_xy.
    """
    if pred_df.empty:
        return gpd.GeoDataFrame(
            columns=["audible", "geometry"],
            geometry="geometry",
            crs=crs,
        )

    audible = spectrum_is_audible(pred_df, ambience)
    return gpd.GeoDataFrame(
        {"audible": audible.astype(int)},
        crs=crs,
        geometry=gpd.points_from_xy(pred_df["Xpos"], pred_df["Ypos"], pred_df["Zpos"]),
    )
