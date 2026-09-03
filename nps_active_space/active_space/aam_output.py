"""Map AAM POI output to NMSim-shaped prediction DataFrames."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd

from aam_translator.bands import band_number_for_frequency
from aam_translator.read_poi import PoiTimeHistory

from nps_active_space.active_space.aam_source import aam_source_id_from_omni
from nps_active_space.active_space.propagation_model import NMSIM_BAND_COLUMNS

__all__ = [
    "aam_source_id_from_omni",
    "poi_history_to_predictions_df",
]


def poi_history_to_predictions_df(
    history: PoiTimeHistory,
    source_pts: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """Map one POI zone (single receiver) to the NMSim TIS-shaped DataFrame."""
    n = history.n_samples
    if n != len(source_pts):
        raise ValueError(
            f"POI has {n} rows but source_pts has {len(source_pts)} points",
        )

    frame = pd.DataFrame({
        "Xpos": source_pts.geometry.x.values,
        "Ypos": source_pts.geometry.y.values,
        "Zpos": source_pts.geometry.z.values,
        "A": history.broadband("dBA"),
    })

    band_index = {bn: i for i, bn in enumerate(history.band_numbers)}
    for col in NMSIM_BAND_COLUMNS:
        band_num = band_number_for_frequency(float(col))
        if band_num in band_index:
            frame[col] = history.band_levels_db[:, band_index[band_num]]
        else:
            frame[col] = np.nan

    return frame
