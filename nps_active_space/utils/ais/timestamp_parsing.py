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

# Alaska standard/daylight offsets from UTC (fixed hours, not IANA zones).
_AKST_UTC_OFFSET = dt.timedelta(hours=9)  # AKST is UTC-9
_AKDT_UTC_OFFSET = dt.timedelta(hours=8)  # AKDT is UTC-8


def _parse_ak_local_with_dst_fallback(
    time: pd.Series,
    primary_format: str,
    primary_offset: dt.timedelta,
) -> pd.Series:
    """Parse Alaska-local timestamps, with a fallback for mixed DST-transition rows."""
    try:
        return pd.to_datetime(time, format=primary_format) + primary_offset
    except ValueError:
        # MXAK files are usually one suffix per file, but transition days can mix
        # AKST and AKDT rows (November AKST→AKDT, March AKDT→AKST).
        akst = time.loc[time.str.endswith("AKST")]
        akdt = time.loc[time.str.endswith("AKDT")]
        parsed = pd.Series(pd.NaT, index=time.index)
        result_dtype = None
        if not akst.empty:
            akst_parsed = (
                pd.to_datetime(akst, format=_MXAK_TIMESTAMP_PATTERNS[4])
                + _AKST_UTC_OFFSET
            )
            parsed.loc[akst.index] = akst_parsed
            result_dtype = akst_parsed.dtype
        if not akdt.empty:
            akdt_parsed = (
                pd.to_datetime(akdt, format=_MXAK_TIMESTAMP_PATTERNS[5])
                + _AKDT_UTC_OFFSET
            )
            parsed.loc[akdt.index] = akdt_parsed
            result_dtype = akdt_parsed.dtype
        # .loc assignment widens to datetime64[ns]; restore pd.to_datetime resolution.
        return parsed.astype(result_dtype)


def parse_mxak_ais_timestamps(time: pd.Series) -> pd.Series:
    """Parse MXAK AIS timestamp strings to UTC-naive datetimes.

    MXAK has shipped several string formats over the years. Within a single
    file every row uses the same format, so we classify from the first row.
    Suffixes ``UTC``, ``AKST``, and ``AKDT`` are honored; bare
    ``YYYY-MM-DD HH:MM:SS`` values are treated as UTC (MXAK regional archive
    convention).
    """
    # Classify format from the first row; all rows in a file share one pattern.
    time_pattern = str(time.iloc[0]).split(" ")
    suffix = time_pattern[-1] if len(time_pattern) >= 2 else ""

    if len(time_pattern) == 5 and suffix == "UTC":
        parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[1])

    elif len(time_pattern) == 3 and suffix == "UTC":
        for fmt in (_MXAK_TIMESTAMP_PATTERNS[2], _MXAK_TIMESTAMP_PATTERNS[3]):
            try:
                parsed = pd.to_datetime(time, format=fmt)
                break
            except ValueError:
                continue
        else:
            raise ValueError("The file's timestamp format is not recognized!")

    elif len(time_pattern) == 3 and suffix in ("AKST", "AKDT"):
        if suffix == "AKST":
            parsed = _parse_ak_local_with_dst_fallback(
                time, _MXAK_TIMESTAMP_PATTERNS[4], _AKST_UTC_OFFSET
            )
        else:
            parsed = _parse_ak_local_with_dst_fallback(
                time, _MXAK_TIMESTAMP_PATTERNS[5], _AKDT_UTC_OFFSET
            )

    elif len(time_pattern) == 2:
        # Bare datetime with no suffix — MXAK treats these as UTC.
        parsed = pd.to_datetime(time, format=_MXAK_TIMESTAMP_PATTERNS[6])

    else:
        raise ValueError("The file's timestamp format is not recognized!")

    if parsed.isna().any():
        bad_count = int(parsed.isna().sum())
        example = time.loc[parsed.isna()].iloc[0]
        raise ValueError(
            f"Failed to parse {bad_count} MXAK AIS timestamp(s); "
            f"first unparseable value: {example!r}"
        )

    return parsed.dt.tz_localize("UTC").dt.tz_localize(None)
