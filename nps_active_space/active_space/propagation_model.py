"""Propagation model interface for active space audibility prediction."""

from __future__ import annotations

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


class PropagationModel(Protocol):
    """Acoustic propagation backend: fixed receiver, many source points in one batch."""

    max_points_per_run: int

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
