"""Map AAM POI output and omni tuning files to NMSim-shaped prediction DataFrames."""

from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from aam_translator.bands import band_number_for_frequency
from aam_translator.read_poi import PoiTimeHistory

from nps_active_space.active_space.propagation_model import NMSIM_BAND_COLUMNS
from nps_active_space.utils.helpers import omni_to_gain

# NMSim omni tuning basename -> AAM NetCDF source id (FLATO*.nc in NCfiles/).
OMNI_TO_AAM_SOURCE_ID: dict[str, str] = {
    "O_+200": "FLATO200",
}


def aam_source_id_from_omni(omni_source: str) -> str:
    stem = Path(omni_source).stem
    if stem in OMNI_TO_AAM_SOURCE_ID:
        return OMNI_TO_AAM_SOURCE_ID[stem]
    if stem.startswith("FLATO"):
        return stem
    # NMSim omni tuning files (O_+005 = +0.5 dB, etc.) map to FLATO200; gain applied in predict().
    if re.fullmatch(r"O_[+-]\d{3}", stem):
        return "FLATO200"
    raise ValueError(
        f"no AAM source_id mapping for omni source {omni_source!r}; "
        f"known keys: {sorted(OMNI_TO_AAM_SOURCE_ID)}",
    )


def apply_omni_gain_offset(frame: pd.DataFrame, omni_source: str) -> pd.DataFrame:
    """Shift prediction levels by the NMSim omni tuning gain (dB)."""
    offset_db = omni_to_gain(omni_source)
    if offset_db == 0:
        return frame
    out = frame.copy()
    out["A"] = out["A"] + offset_db
    for col in NMSIM_BAND_COLUMNS:
        if col in out.columns:
            out[col] = out[col] + offset_db
    return out


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
