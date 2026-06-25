"""Query MXAK AIS archives by date range and study area."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import geopandas as gpd
import pandas as pd

from nps_active_space.utils.ais.reader import MxakAis


def query_ais_mxak(
    ais_path: Path,
    start_date: str,
    end_date: str,
    mask: Optional[gpd.GeoDataFrame] = None,
    mask_buffer_distance: Optional[int] = None,
) -> MxakAis:
    """Query vessel tracks from MXAK AIS CSV archives for a date range and optional area.

    Parameters
    ----------
    ais_path : Path
        Absolute path to a directory with MXAK AIS CSV files (``MXAK-AIS-*-YYYYMMDD.csv``).
    start_date : str
        ISO date string (YYYY-mm-dd) indicating the beginning of the date range to query within.
    end_date : str
        ISO date string (YYYY-mm-dd) indicating the end of the date range to query within, inclusive.
    mask : gpd.GeoDataFrame, default None
        Geopandas.GeoDataFrame instance to spatially filter query results.
    mask_buffer_distance : int, default None
        If provided, buffer ``mask`` by this distance (meters) in Alaska Albers before clipping.

    Returns
    -------
    MxakAis
        Vessel track points from the MXAK archive.
    """
    if mask is not None:
        if mask_buffer_distance:
            ak_albers_mask = mask.to_crs(epsg=3338)
            mask.geometry = ak_albers_mask.buffer(mask_buffer_distance).to_crs(epsg=4326)
        mask = mask.to_crs("epsg:4326")

    all_files = sorted(ais_path.glob("*.csv"))
    assert len(all_files) > 0, f"No AIS files found in {ais_path}"

    start_dt = pd.Timestamp(start_date)
    end_dt = pd.Timestamp(end_date)
    dated_files: list[Path] = []
    for path in all_files:
        match = re.search(r"-(\d{8})(?:-\d+)?\.csv$", path.name)
        if match:
            file_date = pd.Timestamp(match.group(1))
            if start_dt <= file_date <= end_dt:
                dated_files.append(path)

    if dated_files:
        ais = MxakAis(dated_files)
    else:
        ais = MxakAis(ais_path)

    end_ts = pd.Timestamp(end_date) + pd.Timedelta(days=1)
    start_ts = pd.Timestamp(start_date)
    ais = ais.loc[(ais["TIME"] >= start_ts) & (ais["TIME"] < end_ts)]
    ais = ais.loc[~ais.geometry.is_empty]

    if mask is not None:
        ais = gpd.clip(ais, mask)

    assert not ais.empty, (
        f"No AIS data loaded for {start_date} to {end_date}, "
        "please check the time period and/or region mask"
    )
    return ais
