from pathlib import Path
from typing import NamedTuple

import geopandas as gpd
import nps_active_space.utils.config as cfg
from nps_active_space.utils.ais import query_ais_mxak
from nps_active_space.utils.helpers import create_overflights_engine, query_adsb, query_tracks
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.models import Microphone, Tracks
from nps_active_space.utils.time_utils import site_timezone_name, utc_naive_to_site_naive


class GroundTruthingTracks(NamedTuple):
    tracks: Tracks
    faa_path: str | None
    faa_corrections_path: str | None


LoadedTracks = GroundTruthingTracks


def _faa_paths() -> tuple[str, str]:
    return (
        cfg.read('project', 'FAA_Releasable_db'),
        cfg.read('project', 'FAA_type_corrections'),
    )


def _optional_faa_paths(include_faa_paths: bool) -> tuple[str | None, str | None]:
    if include_faa_paths:
        return _faa_paths()
    return None, None


def _load_adsb_tracks(
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    *,
    include_faa_paths: bool,
) -> GroundTruthingTracks:
    raw_tracks = query_adsb(
        adsb_path=cfg.read('data', 'adsb'),
        start_date=start_date,
        end_date=end_date,
        mask=study_area,
    )
    tracks = Tracks(raw_tracks, id_col='flight_id', datetime_col='TIME', z_col='altitude')
    faa_path, faa_corrections_path = _optional_faa_paths(include_faa_paths)
    return GroundTruthingTracks(tracks, faa_path, faa_corrections_path)


def _load_gps_tracks(
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    *,
    include_faa_paths: bool,
) -> GroundTruthingTracks:
    engine = create_overflights_engine(cfg.read('database:overflights'))
    raw_tracks = query_tracks(
        engine=engine,
        start_date=start_date,
        end_date=end_date,
        mask=study_area,
    )
    tracks = Tracks(raw_tracks, 'flight_id', 'ak_datetime', 'altitude_m')
    faa_path, faa_corrections_path = _optional_faa_paths(include_faa_paths)
    return GroundTruthingTracks(tracks, faa_path, faa_corrections_path)


def _load_ais_tracks(
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    microphone: Microphone | None,
) -> GroundTruthingTracks:
    raw_tracks = query_ais_mxak(
        ais_path=Path(cfg.read("data", "ais")),
        start_date=start_date,
        end_date=end_date,
        mask=study_area,
    )
    tracks = Tracks(raw_tracks, id_col='event_id', datetime_col='TIME', z_col='altitude')
    if microphone is not None:
        site_tz = site_timezone_name(microphone.lat, microphone.lon)
        tracks["point_dt"] = utc_naive_to_site_naive(tracks["point_dt"], site_tz)
    return GroundTruthingTracks(tracks, None, None)


def load_tracks(
    source: TrackSource,
    *,
    start_date: str,
    end_date: str,
    study_area: gpd.GeoDataFrame,
    microphone: Microphone | None = None,
    include_faa_paths: bool = True,
) -> GroundTruthingTracks:
    """Load tracks for a deployment window.

    ADSB and GPS return FAA lookup paths when ``include_faa_paths`` is true (the
    default for ground truthing). Viz passes ``False`` so viewing does not
    require FAA config keys.
    """
    match source:
        case TrackSource.ADSB:
            return _load_adsb_tracks(
                start_date, end_date, study_area, include_faa_paths=include_faa_paths
            )
        case TrackSource.GPS:
            return _load_gps_tracks(
                start_date, end_date, study_area, include_faa_paths=include_faa_paths
            )
        case TrackSource.AIS:
            return _load_ais_tracks(start_date, end_date, study_area, microphone)
        case _:
            raise ValueError(f"Unknown track source: {source}")
