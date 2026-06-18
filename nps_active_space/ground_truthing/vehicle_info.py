import geopandas as gpd
import pandas as pd

from nps_active_space.utils.models import FAAReleasable


def load_faa(
    faa_path: str | None,
    faa_corrections_path: str | None,
    database_type: str,
    tracks: gpd.GeoDataFrame,
) -> pd.DataFrame | None:
    """Load FAA releasable data for the tracks in the session, if applicable."""
    if faa_path is None or faa_corrections_path is None:
        return None
    if database_type == "ADSB":
        icaos = tracks["ICAO_address"].unique()
        faa = FAAReleasable(faa_path, faa_corrections_path, icao_addresses=icaos).data
        if faa.empty:
            return None
        return faa
    elif database_type == "GPS":
        track_ids = tracks["track_id"].unique()
        n_numbers = list(map(lambda s: s.split("_")[0][1:], track_ids))
        faa = FAAReleasable(faa_path, faa_corrections_path, n_numbers=n_numbers).data
        if faa.empty:
            return None
        return faa
    return None


def lookup_aircraft(
    faa: pd.DataFrame | None,
    database_type: str,
    track_id: str,
    points: gpd.GeoDataFrame,
) -> tuple[pd.Series | None, str | None, str | None]:
    """Look up FAA row and help text for a track."""
    if faa is None:
        return None, None, None
    faa_row = None
    aircraft_help_text = None
    if database_type == "ADSB":
        icao = points.iloc[0]["ICAO_address"]
        matching = faa[faa["MODE S CODE HEX"] == icao]
        if not matching.empty:
            faa_row = matching.iloc[0]
            aircraft_help_text = f"N-Number: {faa_row['N-NUMBER']}"
    elif database_type == "GPS":
        n_number = track_id.split("_")[0][1:]
        matching = faa[faa["N-NUMBER"] == n_number]
        if not matching.empty:
            faa_row = matching.iloc[0]
            aircraft_help_text = f"ICAO Address: {faa_row["MODE S CODE HEX"]}"
    aircraft_type = faa_row["TYPE AIRCRAFT"] if faa_row is not None else None
    return faa_row, aircraft_help_text, aircraft_type


def lookup_vessel(
    database_type: str,
    track_id: str,
    points: gpd.GeoDataFrame,
) -> tuple[str | None, str | None, str | None]:
    """Look up help text, vessel type, and vessel name for a track."""
    if database_type != "AIS":
        return None, None, None

    mmsi = track_id.split("_")[0]
    help_text = f"MMSI: {mmsi}"

    vessel_name = None
    if "ship_name" in points.columns and points["ship_name"].notna().any():
        vessel_name = points["ship_name"].dropna().iloc[0]

    vessel_type = None
    if "shiptype" in points.columns and points["shiptype"].notna().any():
        vessel_type = points["shiptype"].dropna().iloc[0]

    return help_text, vessel_type, vessel_name
