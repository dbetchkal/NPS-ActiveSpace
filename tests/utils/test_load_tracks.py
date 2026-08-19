from __future__ import annotations

from unittest.mock import MagicMock

import geopandas as gpd
import pandas as pd
import pytest
from pandas.testing import assert_series_equal
from shapely.geometry import box

import nps_active_space.utils.config as cfg
import nps_active_space.utils.load_tracks as load_tracks_mod
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.load_tracks import LoadedTracks, load_tracks
from nps_active_space.utils.models import Adsb, EarlyAdsb, Microphone, Tracks
from nps_active_space.utils.time_utils import site_timezone_name, utc_naive_to_site_naive

UTM8 = "epsg:32608"
START = "2024-05-24"
END = "2024-05-24"


@pytest.fixture
def study_area() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(geometry=[box(500_000, 6_000_000, 510_000, 6_010_000)], crs=UTM8)


@pytest.fixture
def microphone() -> Microphone:
    return Microphone("GLBALSTL2024", lat=58.4, lon=-136.0, z=12.0)


def _points_gdf(
    *,
    id_col: str,
    id_value: str,
    time_col: str,
    z_col: str,
    z_value: float,
) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            id_col: [id_value, id_value],
            time_col: [
                pd.Timestamp("2024-05-24 10:00:00"),
                pd.Timestamp("2024-05-24 10:05:00"),
            ],
            z_col: [z_value, z_value + 100.0],
        },
        geometry=gpd.points_from_xy([500_100.0, 500_200.0], [6_000_100.0, 6_000_200.0]),
        crs=UTM8,
    )


def _faa_config_read(section: str, option: str | None = None) -> str:
    if section == "project" and option == "FAA_Releasable_db":
        return "/faa/MASTER.txt"
    if section == "project" and option == "FAA_type_corrections":
        return "/faa/corrections.json"
    raise KeyError(f"{section}.{option}")


class TestLoadTracks:
    def test_adsb_passes_query_kwargs_and_normalizes_columns(
        self, monkeypatch: pytest.MonkeyPatch, study_area: gpd.GeoDataFrame, microphone: Microphone
    ) -> None:
        raw = _points_gdf(
            id_col="flight_id",
            id_value="abc123",
            time_col="TIME",
            z_col="altitude",
            z_value=5000.0,
        )
        calls: list[dict] = []

        def fake_read(section: str, option: str | None = None) -> str:
            if section == "data" and option == "adsb":
                return "/data/adsb"
            return _faa_config_read(section, option)

        def fake_query_adsb(**kwargs) -> Adsb:
            calls.append(kwargs)
            adsb = Adsb.__new__(Adsb)
            gpd.GeoDataFrame.__init__(adsb, raw)
            return adsb

        monkeypatch.setattr(cfg, "read", fake_read)
        monkeypatch.setattr(load_tracks_mod, "query_adsb", fake_query_adsb)

        loaded = load_tracks(
            TrackSource.ADSB,
            start_date=START,
            end_date=END,
            study_area=study_area,
            microphone=microphone,
        )

        assert calls == [
            {
                "adsb_path": "/data/adsb",
                "start_date": START,
                "end_date": END,
                "mask": study_area,
            }
        ]
        assert isinstance(loaded.tracks, Tracks)
        assert loaded.faa_path == "/faa/MASTER.txt"
        assert loaded.tracks["track_id"].tolist() == ["abc123", "abc123"]
        assert list(loaded.tracks["z"]) == pytest.approx([5000.0, 5100.0])
        assert_series_equal(loaded.tracks["point_dt"], raw["TIME"], check_names=False)

    def test_adsb_early_format_keeps_parser_timestamps(
        self, monkeypatch: pytest.MonkeyPatch, study_area: gpd.GeoDataFrame, microphone: Microphone
    ) -> None:
        raw = _points_gdf(
            id_col="flight_id",
            id_value="legacy1",
            time_col="TIME",
            z_col="altitude",
            z_value=5000.0,
        )
        local_times = raw["TIME"].copy()

        def fake_read(section: str, option: str | None = None) -> str:
            if section == "data" and option == "adsb":
                return "/data/adsb"
            return _faa_config_read(section, option)

        def fake_query_adsb(**kwargs) -> EarlyAdsb:
            early = EarlyAdsb.__new__(EarlyAdsb)
            gpd.GeoDataFrame.__init__(early, raw)
            return early

        monkeypatch.setattr(cfg, "read", fake_read)
        monkeypatch.setattr(load_tracks_mod, "query_adsb", fake_query_adsb)

        loaded = load_tracks(
            TrackSource.ADSB,
            start_date=START,
            end_date=END,
            study_area=study_area,
            microphone=microphone,
        )

        assert_series_equal(loaded.tracks["point_dt"], local_times, check_names=False)

    def test_gps_passes_query_kwargs(
        self, monkeypatch: pytest.MonkeyPatch, study_area: gpd.GeoDataFrame, microphone: Microphone
    ) -> None:
        raw = _points_gdf(
            id_col="flight_id",
            id_value="gps1",
            time_col="ak_datetime",
            z_col="altitude_m",
            z_value=8000.0,
        )
        engine = MagicMock()
        track_calls: list[dict] = []

        def fake_read(section: str, option: str | None = None):
            if section == "database:overflights" and option is None:
                return {
                    "name": "db",
                    "host": "localhost",
                    "username": "u",
                    "password": "p",
                    "port": 5432,
                }
            return _faa_config_read(section, option)

        def fake_query_tracks(**kwargs) -> gpd.GeoDataFrame:
            track_calls.append(kwargs)
            return raw

        monkeypatch.setattr(cfg, "read", fake_read)
        monkeypatch.setattr(load_tracks_mod, "create_overflights_engine", lambda _cfg: engine)
        monkeypatch.setattr(load_tracks_mod, "query_tracks", fake_query_tracks)

        loaded = load_tracks(
            TrackSource.GPS,
            start_date=START,
            end_date=END,
            study_area=study_area,
            microphone=microphone,
        )

        assert track_calls[0]["engine"] is engine
        assert track_calls[0]["start_date"] == START
        assert track_calls[0]["mask"] is study_area
        assert loaded.tracks["track_id"].iloc[0] == "gps1"
        assert loaded.faa_path == "/faa/MASTER.txt"
        assert loaded.faa_corrections_path == "/faa/corrections.json"

    def test_adsb_skips_faa_paths_when_disabled(
        self, monkeypatch: pytest.MonkeyPatch, study_area: gpd.GeoDataFrame
    ) -> None:
        raw = _points_gdf(
            id_col="flight_id",
            id_value="abc123",
            time_col="TIME",
            z_col="altitude",
            z_value=5000.0,
        )

        def fake_read(section: str, option: str | None = None) -> str:
            if section == "data" and option == "adsb":
                return "/data/adsb"
            raise KeyError(f"{section}.{option}")

        def fake_query_adsb(**kwargs) -> Adsb:
            adsb = Adsb.__new__(Adsb)
            gpd.GeoDataFrame.__init__(adsb, raw)
            return adsb

        monkeypatch.setattr(cfg, "read", fake_read)
        monkeypatch.setattr(load_tracks_mod, "query_adsb", fake_query_adsb)

        loaded = load_tracks(
            TrackSource.ADSB,
            start_date=START,
            end_date=END,
            study_area=study_area,
            include_faa_paths=False,
        )

        assert loaded.faa_path is None
        assert loaded.faa_corrections_path is None

    def test_ais_without_microphone_keeps_utc_timestamps(
        self, monkeypatch: pytest.MonkeyPatch, study_area: gpd.GeoDataFrame
    ) -> None:
        raw = _points_gdf(
            id_col="event_id",
            id_value="evt1",
            time_col="TIME",
            z_col="altitude",
            z_value=0.0,
        )
        utc_times = raw["TIME"].copy()

        monkeypatch.setattr(cfg, "read", lambda section, option=None: "/data/ais")
        monkeypatch.setattr(load_tracks_mod, "query_ais_mxak", lambda **kwargs: raw)

        loaded = load_tracks(
            TrackSource.AIS,
            start_date=START,
            end_date=END,
            study_area=study_area,
        )

        assert_series_equal(loaded.tracks["point_dt"], utc_times, check_names=False)

    def test_ais_converts_timestamps_to_site_local(
        self, monkeypatch: pytest.MonkeyPatch, study_area: gpd.GeoDataFrame, microphone: Microphone
    ) -> None:
        raw = _points_gdf(
            id_col="event_id",
            id_value="evt1",
            time_col="TIME",
            z_col="altitude",
            z_value=0.0,
        )
        utc_times = raw["TIME"].copy()
        ais_calls: list[dict] = []

        monkeypatch.setattr(cfg, "read", lambda section, option=None: "/data/ais")
        monkeypatch.setattr(
            load_tracks_mod,
            "query_ais_mxak",
            lambda **kwargs: ais_calls.append(kwargs) or raw,
        )

        loaded = load_tracks(
            TrackSource.AIS,
            start_date=START,
            end_date=END,
            study_area=study_area,
            microphone=microphone,
        )

        assert ais_calls[0]["ais_path"].as_posix() == "/data/ais"
        assert ais_calls[0]["mask"] is study_area
        site_tz = site_timezone_name(microphone.lat, microphone.lon)
        expected = utc_naive_to_site_naive(utc_times, site_tz)
        assert_series_equal(loaded.tracks["point_dt"], expected, check_names=False)
        assert loaded.faa_path is None
