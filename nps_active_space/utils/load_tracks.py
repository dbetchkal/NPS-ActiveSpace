from pathlib import Path
from typing import NamedTuple

import geopandas as gpd
import nps_active_space.utils.config as cfg
from nps_active_space.utils.ais import query_ais_mxak
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.helpers import create_overflights_engine, query_adsb, query_tracks
from nps_active_space.utils.models import Adsb, EarlyAdsb, Microphone, Tracks
from nps_active_space.utils.time_utils import site_timezone_name, utc_naive_to_site_naive


class LoadedTracks(NamedTuple):
    """Tracks for a deployment window plus optional FAA lookup paths (aircraft only)."""

    tracks: Tracks
    faa_path: str | None
    faa_corrections_path: str | None


def _faa_paths() -> tuple[str, str]:
    return (
        cfg.read("project", "FAA_Releasable_db"),
        cfg.read("project", "FAA_type_corrections"),
    )


def _apply_site_local_timestamps(tracks: Tracks, microphone: Microphone | None) -> None:
    """Convert UTC-naive ``point_dt`` to site-local naive using deployment lat/lon."""
    if microphone is None:
        return
    site_tz = site_timezone_name(microphone.lat, microphone.lon)
    tracks["point_dt"] = utc_naive_to_site_naive(tracks["point_dt"], site_tz)


def _load_adsb_tracks(
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    *,
    include_faa_paths: bool,
) -> LoadedTracks:
    raw_tracks = query_adsb(
        adsb_path=cfg.read("data", "adsb"),
        start_date=start_date,
        end_date=end_date,
        mask=study_area,
    )
    tracks = Tracks(raw_tracks, id_col="flight_id", datetime_col="TIME", z_col="altitude")
    # TSV Unix → naive datetime is used on the NVSPL clock as-is (ground-truthing convention).
    # Legacy .txt EarlyAdsb is also naive local. Do not apply site timezone shift for ADSB.
    if not isinstance(raw_tracks, (Adsb, EarlyAdsb)):
        raise TypeError(f"Unexpected ADSB query result type: {type(raw_tracks)!r}")
    if include_faa_paths:
        faa_path, faa_corrections_path = _faa_paths()
        return LoadedTracks(tracks, faa_path, faa_corrections_path)
    return LoadedTracks(tracks, None, None)


def _load_gps_tracks(
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    *,
    include_faa_paths: bool,
) -> LoadedTracks:
    engine = create_overflights_engine(cfg.read("database:overflights"))
    raw_tracks = query_tracks(
        engine=engine,
        start_date=start_date,
        end_date=end_date,
        mask=study_area,
    )
    tracks = Tracks(raw_tracks, id_col="flight_id", datetime_col="ak_datetime", z_col="altitude_m")
    if include_faa_paths:
        faa_path, faa_corrections_path = _faa_paths()
        return LoadedTracks(tracks, faa_path, faa_corrections_path)
    return LoadedTracks(tracks, None, None)


def _load_ais_tracks(
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    microphone: Microphone | None,
) -> LoadedTracks:
    raw_tracks = query_ais_mxak(
        ais_path=Path(cfg.read("data", "ais")),
        start_date=start_date,
        end_date=end_date,
        mask=study_area,
    )
    tracks = Tracks(raw_tracks, id_col="event_id", datetime_col="TIME", z_col="altitude")
    _apply_site_local_timestamps(tracks, microphone)
    return LoadedTracks(tracks, None, None)


def load_tracks(
    source: TrackSource,
    *,
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    microphone: Microphone | None = None,
    include_faa_paths: bool = True,
) -> LoadedTracks:
    """Load tracks for a deployment window.

    ADSB TSV and legacy ``.txt`` keep parser naive ``point_dt`` on the NVSPL clock (no site
    timezone shift). AIS converts UTC-naive timestamps to site-local naive when
    ``microphone`` is provided. GPS ``ak_datetime`` is already site-local from the
    overflights database.

    ADSB and GPS return FAA lookup paths when ``include_faa_paths`` is true (the
    default for ground truthing). Spatial filtering uses ``study_area`` only;
    ``microphone`` is not used for extent.
    """
    match source:
        case TrackSource.ADSB:
            return _load_adsb_tracks(
                start_date,
                end_date,
                study_area,
                include_faa_paths=include_faa_paths,
            )
        case TrackSource.GPS:
            return _load_gps_tracks(
                start_date, end_date, study_area, include_faa_paths=include_faa_paths
            )
        case TrackSource.AIS:
            return _load_ais_tracks(start_date, end_date, study_area, microphone)
        case _:
            raise ValueError(f"Unknown track source: {source}")
