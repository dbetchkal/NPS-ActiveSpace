"""
Functions for reading from Alaska Marine Exchange (MXAK) AIS CSV archives.

Also included in this module is the :class:`MxakAis` class, which is a 
GeoDataFrame wrapper for the MXAK AIS CSV archives.
"""

from __future__ import annotations

import concurrent.futures
from pathlib import Path
from typing import Any, TypeAlias

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from nps_active_space.utils.ais.timestamp_parsing import parse_mxak_ais_timestamps

MODERN_TO_LEGACY_COLUMNS = {
    "bs_ts": "Base station time stamp",
    "mmsi": "MMSI",
    "imo": "IMO",
    "name": "Ship name",
    "nav_status": "Navigational status (text)",
    "lat": "Latitude",
    "lon": "Longitude",
    "cog": "Course over ground",
    "sog": "Speed over ground",
    "destination": "Destination",
    "shiptype": "Type of ship (text)",
    "draught": "Draught",
}

# Modern-only fields dropped when mapping to the legacy column layout.
MODERN_DROP_COLUMNS = ("callsign", "eta", "length", "width")

# Legacy fields dropped during normalization to the canonical track schema.
LEGACY_DROP_COLUMNS = (
    "IMO",
    "Size A",
    "Size B",
    "Size C",
    "Size D",
    "Draught",
    "Destination",
    "Heading",
    "Navigational status (text)",
    "Country (AIS)",
    "Target class (text)",
    "Data source type (text)",
    "Data source region",
)

LEGACY_TIMESTAMP_COLUMN = "Base station time stamp"

LEGACY_RENAME_TO_CANONICAL = {
    LEGACY_TIMESTAMP_COLUMN: "TIME",
    "Ship name": "ship_name",
    "Type of ship (text)": "shiptype",
    # Course over ground in MXAK exports; ``heading`` matches legacy ADSB/AIS schema names.
    "Course over ground": "heading",
    "Latitude": "lat",
    "Longitude": "lon",
    "Speed over ground": "velocity",
}

POINT_COLUMNS = ("MMSI", "TIME", "lat", "lon")
EVENT_GAP_SECONDS = 15 * 60
MIN_EVENT_POINTS = 3
WGS84 = "epsg:4326"

_SEGMENT_SCRATCH_COLUMNS = ("dur_secs", "diff_event", "cumsum")

MxakAisCsvPathInput: TypeAlias = str | Path | list[str | Path]


def parse_mxak_ais_csv_points(
    ais_file: str | Path, *, usecols: list[str] | None = None
) -> pd.DataFrame:
    """Parse one MXAK AIS CSV to normalized UTC point rows (no ``event_id`` segmentation)."""
    ais_file = _coerce_path(ais_file)
    raw_csv = pd.read_csv(
        ais_file,
        engine="c",
        usecols=usecols,
        low_memory=False,
    )

    if "mmsi" in raw_csv.columns:
        raw_csv = _normalize_modern_mxak_ais_headers(raw_csv)

    # 1090 MHz jet ADS-B points mixed into MXAK feeds carry MMSI < 1e8; drop them.
    mxak_points = raw_csv.loc[raw_csv["MMSI"] >= 100000000, :].copy()

    header_row_mask = mxak_points.iloc[:, 0].isin([LEGACY_TIMESTAMP_COLUMN])
    mxak_points = mxak_points[~header_row_mask]
    if LEGACY_TIMESTAMP_COLUMN not in mxak_points.columns:
        raise KeyError(
            f"MXAK AIS CSV {ais_file} lacks {LEGACY_TIMESTAMP_COLUMN!r} "
            f"(columns: {list(raw_csv.columns)})"
        )

    mxak_points = mxak_points.rename(columns=LEGACY_RENAME_TO_CANONICAL)
    mxak_points.drop(list(LEGACY_DROP_COLUMNS), axis=1, inplace=True, errors="ignore")

    mxak_points.drop_duplicates(inplace=True)
    mxak_points.dropna(subset=list(POINT_COLUMNS), how="any", axis=0, inplace=True)

    if mxak_points.empty:
        return mxak_points

    # Vessel z is fixed at sea level for now; tide-aware heights would use
    # NOAA CO-OPS (https://api.tidesandcurrents.noaa.gov/api/prod/).
    mxak_points["altitude"] = 0.0
    mxak_points["TIME"] = parse_mxak_ais_timestamps(mxak_points["TIME"])
    mxak_points["DATE"] = mxak_points["TIME"].dt.strftime("%Y%m%d")

    mxak_points.sort_values(["MMSI", "TIME"], inplace=True, ignore_index=True)
    return mxak_points


def segment_mxak_ais_events(sorted_mxak_points: pd.DataFrame) -> pd.DataFrame:
    """Assign ``event_id`` and drop short transits (requires sorted MMSI/TIME points)."""
    mxak_events = sorted_mxak_points.copy()
    mxak_events["dur_secs"] = mxak_events.groupby("MMSI")["TIME"].diff().dt.total_seconds()
    mxak_events["dur_secs"] = mxak_events["dur_secs"].fillna(0)

    mxak_events.drop_duplicates(subset=list(POINT_COLUMNS), keep="last", inplace=True)

    mxak_events["diff_event"] = mxak_events["dur_secs"] >= EVENT_GAP_SECONDS
    mxak_events["cumsum"] = mxak_events.groupby("MMSI")["diff_event"].cumsum()
    # The event_id for a transit includes the date, so if a transit crosses the 
    # midnight (00:00) boundary (for its timezone) it will be treated as two separate events. 
    mxak_events["event_id"] = (
        mxak_events["MMSI"].astype(str)
        + "_"
        + mxak_events["cumsum"].astype(str)
        + "_"
        + mxak_events["DATE"].astype(str)
    )

    # One position per vessel per second within an event (mirrors ADSB TIME+ICAO dedup).
    mxak_events.drop_duplicates(subset=["event_id", "TIME"], keep="last", inplace=True)
    mxak_events = mxak_events[
        mxak_events.groupby("event_id").event_id.transform(len) >= MIN_EVENT_POINTS
    ]
    return mxak_events.drop(columns=list(_SEGMENT_SCRATCH_COLUMNS))


def parse_mxak_ais_csv_events(
    ais_file: str | Path, *, usecols: list[str] | None = None
) -> pd.DataFrame:
    """Parse one MXAK AIS CSV to segmented transit rows with ``event_id``."""
    mxak_points = parse_mxak_ais_csv_points(ais_file, usecols=usecols)
    if mxak_points.empty:
        return mxak_points
    return segment_mxak_ais_events(mxak_points)


def load_mxak_ais_geodataframe(
    file_locations: MxakAisCsvPathInput,
) -> gpd.GeoDataFrame:
    """Load MXAK AIS from a directory or CSV path list into a WGS84 GeoDataFrame."""
    csv_paths = _coerce_mxak_ais_path_input(file_locations)
    if isinstance(csv_paths, Path):
        if not csv_paths.is_dir():
            raise FileNotFoundError(f"MXAK AIS directory not found: {csv_paths}")
        csv_paths = sorted(csv_paths.glob("*.csv"))
    elif not _is_mxak_ais_csv_file_list(csv_paths):
        raise ValueError("Expected a list of existing .csv MXAK AIS file paths.")

    if not csv_paths:
        raise ValueError(f"No MXAK AIS CSV files found at {file_locations}")

    with concurrent.futures.ThreadPoolExecutor() as pool:
        per_file_events = list(
            tqdm(
                pool.map(parse_mxak_ais_csv_events, csv_paths),
                total=len(csv_paths),
                unit="MXAK AIS files",
            )
        )

    mxak_events = pd.concat(per_file_events)
    mxak_events.drop_duplicates(subset=["TIME", "MMSI"], keep="last", inplace=True)
    return gpd.GeoDataFrame(
        mxak_events,
        geometry=gpd.points_from_xy(mxak_events["lon"], mxak_events["lat"]),
        crs=WGS84,
    )


def _coerce_path(path: str | Path) -> Path:
    return path if isinstance(path, Path) else Path(path)


def _normalize_modern_mxak_ais_headers(raw_csv: pd.DataFrame) -> pd.DataFrame:
    """Rename modern MXAK headers to the legacy layout used by canonicalization."""
    renamed = raw_csv.rename(columns=MODERN_TO_LEGACY_COLUMNS)
    return renamed.drop(columns=list(MODERN_DROP_COLUMNS), errors="ignore")


def _coerce_mxak_ais_path_input(file_locations: MxakAisCsvPathInput) -> Path | list[Path]:
    if isinstance(file_locations, list):
        return [_coerce_path(path) for path in file_locations]
    return _coerce_path(file_locations)


def _is_mxak_ais_csv_file_list(csv_paths: list[Path]) -> bool:
    return all(path.suffix == ".csv" and path.is_file() for path in csv_paths)


class MxakAis(gpd.GeoDataFrame):
    """GeoDataFrame wrapper for MXAK daily AIS CSV exports:
     https://www.mxak.org/services/vessel-tracking/.

    Parameters
    ----------
    data : str, Path, list[str | Path], or gpd.GeoDataFrame
        Directory of ``MXAK-AIS-*-YYYYMMDD.csv`` files, a list of such files,
        or an existing GeoDataFrame of MXAK AIS points.
    """

    def __init__(self, data: Any = None, *args: Any, **kwargs: Any):
        if isinstance(data, gpd.GeoDataFrame):
            data = data.to_crs(WGS84)
        elif isinstance(data, (str, Path)):
            data = load_mxak_ais_geodataframe(data)
        elif isinstance(data, list):
            if not data:
                raise ValueError("Expected at least one MXAK AIS CSV path.")
            if not all(isinstance(path, (str, Path)) for path in data):
                raise TypeError(
                    "MxakAis CSV path lists must contain only str or Path entries."
                )
            data = load_mxak_ais_geodataframe(data)
        super().__init__(data, *args, **kwargs)
