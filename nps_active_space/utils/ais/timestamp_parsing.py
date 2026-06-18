"""MXAK AIS timestamp string parsing."""

from __future__ import annotations

import datetime as dt

import pandas as pd

_MXAK_TIMESTAMP_PATTERNS = {
    1: "%d %b %Y %H:%M:%S UTC",
    2: "%Y-%m-%d %H:%M:%S UTC",
    3: "%m-%d-%Y %H:%M:%S UTC",
    4: "%Y-%m-%d %H:%M:%S AKST",
    5: "%Y-%m-%d %H:%M:%S AKDT",
    6: "%Y-%m-%d %H:%M:%S",
}


def parse_mxak_ais_timestamps(time: pd.Series) -> pd.Series:
    """Parse MXAK AIS timestamp strings to UTC-naive datetimes.

    MXAK files use a consistent format within each file. Suffixes ``UTC``,
    ``AKST``, and ``AKDT`` are honored; bare ``YYYY-MM-DD HH:MM:SS`` values
    are treated as UTC (MXAK regional archive convention).
    """
    sample = str(time.iloc[0])
    time_pattern = sample.split(" ")

    if (len(time_pattern) == 5) and (time_pattern[-1] == "UTC"):
        parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[1])

    elif (len(time_pattern) == 3) and (time_pattern[-1] == "UTC"):
        try:
            parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[2])
        except ValueError:
            parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[3])

    elif (len(time_pattern) == 3) and (time_pattern[-1] == "AKST"):
        try:
            parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[4]) + dt.timedelta(hours=9)
        except ValueError:
            akst = time.loc[time.str.endswith("AKST")]
            akdt = time.loc[time.str.endswith("AKDT")]
            parsed = pd.Series(index=time.index, dtype="datetime64[ns]")
            parsed.loc[akst.index] = (
                pd.to_datetime(akst, format=_MXAK_TIMESTAMP_PATTERNS[4]) + dt.timedelta(hours=9)
            )
            parsed.loc[akdt.index] = (
                pd.to_datetime(akdt, format=_MXAK_TIMESTAMP_PATTERNS[5]) + dt.timedelta(hours=8)
            )

    elif (len(time_pattern) == 3) and (time_pattern[-1] == "AKDT"):
        try:
            parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[5]) + dt.timedelta(hours=8)
        except ValueError:
            akst = time.loc[time.str.endswith("AKST")]
            akdt = time.loc[time.str.endswith("AKDT")]
            parsed = pd.Series(index=time.index, dtype="datetime64[ns]")
            parsed.loc[akst.index] = (
                pd.to_datetime(akst, format=_MXAK_TIMESTAMP_PATTERNS[4]) + dt.timedelta(hours=9)
            )
            parsed.loc[akdt.index] = (
                pd.to_datetime(akdt, format=_MXAK_TIMESTAMP_PATTERNS[5]) + dt.timedelta(hours=8)
            )

    elif len(time_pattern) == 2:
        parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[6])

    else:
        raise ValueError("The file's timestamp format is not recognized!")

    return parsed.dt.tz_localize("UTC").dt.tz_localize(None)
