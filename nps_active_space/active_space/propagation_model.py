"""Propagation model interface for active space audibility prediction."""

from __future__ import annotations

import os
from typing import Protocol, TYPE_CHECKING

import geopandas as gpd
import pandas as pd

if TYPE_CHECKING:
    from nps_active_space.utils.models import Microphone

# 1/3-octave band columns expected by ``ActiveSpaceGenerator._find_audible_points``.
NMSIM_BAND_COLUMNS: tuple[str, ...] = (
    "10", "12.5", "15.8", "20", "25", "31.5", "40", "50", "63", "80", "100", "125",
    "160", "200", "250", "315", "400", "500", "630", "800", "1000", "1250", "1600",
    "2000", "2500", "3150", "4000", "5000", "6300", "8000", "10000", "12500",
)

# Site-relative dirs for incremental prediction caches (CSV keyed by alt/omni/heading).
NMSIM_OUTPUT_SUBDIR = "Output_Data/nmsim"
NMSIM_PREDICTIONS_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/predictions"
NMSIM_SCRATCH_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/scratch"

# AAM layout (see docs/aam_integration_notes.md):
#   Input_Data/AAM/terrain/{mic}/       — cached ELV/IMP (model input, long-lived)
#   Output_Data/aam/predictions/        — incremental spectral cache CSVs
#   Output_Data/aam/runs/{job}/         — per-batch .inp / .POI / run log (scratch)
AAM_INPUT_SUBDIR = "Input_Data/AAM"
AAM_TERRAIN_SUBDIR = f"{AAM_INPUT_SUBDIR}/terrain"
AAM_OUTPUT_SUBDIR = "Output_Data/aam"
AAM_PREDICTIONS_SUBDIR = f"{AAM_OUTPUT_SUBDIR}/predictions"
AAM_RUNS_SUBDIR = f"{AAM_OUTPUT_SUBDIR}/runs"
AAM_RUN_LOG_FILENAME = "active_space.log"


def prediction_cache_csv_path(
    root_dir: str,
    predictions_subdir: str,
    altitude_m: int,
    omni_stem: str,
    heading: int,
) -> str:
    """Path to the merged prediction cache CSV for one alt/omni/heading combination."""
    cache_dir = os.path.join(root_dir, predictions_subdir)
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, f"{altitude_m}m_{omni_stem}_{heading}deg.csv")


class PropagationModel(Protocol):
    """Acoustic propagation backend: fixed receiver, many source points in one batch."""

    max_points_per_run: int
    predictions_subdir: str

    def prepare_site(
        self,
        dem_src: str,
        study_area: gpd.GeoDataFrame,
        mic: "Microphone",
        *,
        project_dem: bool = True,
        suffix: str = "",
    ) -> object:
        """One-time per-mic terrain and receiver setup."""
        ...

    def predict(
        self,
        site: object,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None = None,
    ) -> pd.DataFrame:
        """Return spectral predictions with columns Xpos, Ypos, Zpos, A, and NMSIM band labels."""
        ...
