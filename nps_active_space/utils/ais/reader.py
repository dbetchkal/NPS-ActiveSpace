"""Alaska Marine Exchange (MXAK) AIS CSV archives."""

from __future__ import annotations

import concurrent.futures
import glob
import os
from os import PathLike
from typing import Any

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from nps_active_space.utils.ais.timestamp_parsing import parse_mxak_ais_timestamps

# Modern MXAK CSV headers (lowercase ``bs_ts``, ``mmsi``, …).
MODERN_CSV_COLUMNS = (
    "bs_ts",
    "mmsi",
    "callsign",
    "imo",
    "name",
    "nav_status",
    "lat",
    "lon",
    "cog",
    "sog",
    "destination",
    "eta",
    "shiptype",
    "draught",
    "length",
    "width",
)

# Legacy MXAK CSV headers; modern files are renamed to this set before normalization.
LEGACY_CSV_COLUMNS = (
    "Base station time stamp",
    "MMSI",
    "callsign",
    "IMO",
    "Ship name",
    "Navigational status (text)",
    "Latitude",
    "Longitude",
    "Course over ground",
    "Speed over ground",
    "Destination",
    "eta",
    "Type of ship (text)",
    "Draught",
    "length",
    "width",
)

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
    "Course over ground": "heading",
    "Latitude": "lat",
    "Longitude": "lon",
    "Speed over ground": "velocity",
}

POINT_COLUMNS = ("MMSI", "TIME", "lat", "lon")
EVENT_GAP_SECONDS = 900
MIN_EVENT_POINTS = 3
WGS84 = "epsg:4326"


def _is_csv_path_list(data: Any) -> bool:
    return (
        isinstance(data, list)
        and bool(data)
        and all(
            isinstance(path, (str, PathLike))
            and str(path).endswith(".csv")
            and os.path.isfile(path)
            for path in data
        )
    )


def _is_mxak_csv_input(data: Any) -> bool:
    return isinstance(data, str) or _is_csv_path_list(data)


class MxakAis(gpd.GeoDataFrame):
    """GeoDataFrame wrapper for MXAK daily AIS CSV exports (www.mxak.org).

    Parameters
    ----------
    data : list[str], str, or gpd.GeoDataFrame
        Directory of ``MXAK-AIS-*-YYYYMMDD.csv`` files, a list of such files,
        or an existing GeoDataFrame of MXAK AIS points.
    """

    def __init__(self, data: Any = None, *args: Any, **kwargs: Any):
        if isinstance(data, gpd.GeoDataFrame):
            data = data.to_crs(WGS84)
        elif _is_mxak_csv_input(data):
            data = self._read_paths(data)
        super().__init__(data, *args, **kwargs)

    def parse_mxak_file(self, ais_file_entry: str | PathLike, state: tuple = (None, None, 1)):
        """Parse one MXAK daily CSV into a normalized point table.

        ``state`` is accepted for ``ThreadPoolExecutor.map`` compatibility; only
        the filepath argument is used (``usecols=None`` reads all columns).
        """
        timestamps, columns, index_index = state

        df = pd.read_csv(
            str(ais_file_entry),
            engine="c",
            usecols=columns,
            low_memory=False,
        )

        if "mmsi" in df.columns:
            df.columns = list(LEGACY_CSV_COLUMNS)
            df.drop(list(MODERN_DROP_COLUMNS), axis=1, inplace=True, errors="ignore")

        # 1090 MHz jet ADS-B points mixed into MXAK feeds carry MMSI < 1e8; drop them.
        df = df.loc[df["MMSI"] >= 100000000, :].copy()

        mask = df.iloc[:, 0].isin([LEGACY_TIMESTAMP_COLUMN])
        df = df[~mask]
        if LEGACY_TIMESTAMP_COLUMN not in df.columns:
            raise KeyError

        df = df.rename(columns=LEGACY_RENAME_TO_CANONICAL)
        df.drop(list(LEGACY_DROP_COLUMNS), axis=1, inplace=True, errors="ignore")

        df.drop_duplicates(inplace=True)
        df.dropna(subset=list(POINT_COLUMNS), how="any", axis=0, inplace=True)

        if len(df) == 0:
            return df

        # Vessel z is fixed at sea level for now; tide-aware heights would use
        # NOAA CO-OPS (https://api.tidesandcurrents.noaa.gov/api/prod/).
        df["altitude"] = 0.0
        df["TIME"] = parse_mxak_ais_timestamps(df["TIME"])
        df["DATE"] = df["TIME"].dt.strftime("%Y%m%d")

        df.sort_values(["MMSI", "TIME"], inplace=True, ignore_index=True)

        df["dur_secs"] = df.groupby("MMSI")["TIME"].diff().dt.total_seconds()
        df["dur_secs"] = df["dur_secs"].fillna(0)

        df.drop_duplicates(subset=list(POINT_COLUMNS), keep="last")

        df["diff_event"] = df["dur_secs"] >= EVENT_GAP_SECONDS
        df["cumsum"] = df.groupby("MMSI")["diff_event"].cumsum()
        df["event_id"] = (
            df["MMSI"].astype(str)
            + "_"
            + df["cumsum"].astype(str)
            + "_"
            + df["DATE"].astype(str)
        )

        # One position per vessel per second within an event (mirrors ADSB TIME+ICAO dedup).
        df.drop_duplicates(subset=["event_id", "TIME"], keep="last", inplace=True)
        df = df[df.groupby("event_id").event_id.transform(len) >= MIN_EVENT_POINTS]

        return df

    def _read_paths(self, path_or_paths: str | list[str]) -> gpd.GeoDataFrame:
        """Read MXAK AIS points from a directory or explicit ``.csv`` file list.

        Parameters
        ----------
        path_or_paths : str or list[str]
            Directory of MXAK CSV files or explicit paths to ``.csv`` files.

        Returns
        -------
        gpd.GeoDataFrame
            WGS84 point geometry with normalized MXAK columns.

        Raises
        ------
        AssertionError
            If a path does not exist or a file is not ``.csv``.
        KeyError
            If a CSV lacks the expected MXAK timestamp column.
        """
        if isinstance(path_or_paths, str):
            assert os.path.isdir(path_or_paths), f"{path_or_paths} does not exist."
            csv_paths = glob.glob(os.path.join(path_or_paths, "*.csv"))
        else:
            assert _is_csv_path_list(path_or_paths), (
                "Expected a list of existing .csv MXAK AIS file paths."
            )
            csv_paths = path_or_paths

        with concurrent.futures.ThreadPoolExecutor() as pool:
            parts = list(
                tqdm(
                    pool.map(self.parse_mxak_file, csv_paths),
                    total=len(csv_paths),
                    unit="MXAK AIS files",
                )
            )

        data = pd.concat(parts)
        data.drop_duplicates(subset=["TIME", "MMSI"], keep="last", inplace=True)
        return gpd.GeoDataFrame(
            data,
            geometry=gpd.points_from_xy(data["lon"], data["lat"]),
            crs=WGS84,
        )
