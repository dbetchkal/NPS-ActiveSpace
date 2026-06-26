"""Map monitoring-site coordinates to NVSPL wall-clock time."""

from __future__ import annotations

from zoneinfo import ZoneInfo

import pandas as pd
from timezonefinder import TimezoneFinder

_TZ_FINDER = TimezoneFinder()


def site_timezone_name(lat: float, lon: float) -> str:
    """Return the IANA timezone name for a WGS84 monitoring-site location."""
    tz_name = _TZ_FINDER.timezone_at(lat=lat, lng=lon)
    if not tz_name:
        raise ValueError(f"Could not determine timezone for site at ({lat}, {lon}).")
    return tz_name


def utc_naive_to_site_naive(timestamps: pd.Series, site_tz: str) -> pd.Series:
    """Convert UTC-naive timestamps to site-local naive (NVSPL clock)."""
    utc = pd.to_datetime(timestamps).dt.tz_localize("UTC")
    return utc.dt.tz_convert(ZoneInfo(site_tz)).dt.tz_localize(None)
