from types import SimpleNamespace
from unittest.mock import patch

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point, box

from nps_active_space.scripts.run_audible_transits import (
    AudibleTransitsAIS,
    init_audible_transits,
)
from nps_active_space.utils.models import Tracks


def _ais_metadata() -> dict:
    return {
        "unit": "GLBA",
        "site": "TEST",
        "activespace year": 2025,
        "gain": -20.0,
        "study start": "2025-01-07",
        "study end": "2025-01-07",
        "database type": "AIS",
    }


def _ais_paths(tmp_path) -> dict:
    return {"project": str(tmp_path), "AIS": str(tmp_path / "ais")}


def _ais_track_points() -> Tracks:
    raw = gpd.GeoDataFrame(
        {
            "event_id": ["368018710_0_20250107", "368018710_0_20250107"],
            "TIME": [
                pd.Timestamp("2025-01-07 04:58:04"),
                pd.Timestamp("2025-01-07 05:08:04"),
            ],
            "altitude": [0.0, 0.0],
            "MMSI": [368018710, 368018710],
            "shiptype": ["Tug", "Tug"],
        },
        geometry=[Point(-134.97, 58.40), Point(-134.98, 58.41)],
        crs="EPSG:4326",
    )
    return Tracks(raw, id_col="event_id", datetime_col="TIME", z_col="altitude")


class TestInitAudibleTransitsAis:
    def test_returns_ais_subclass(self, tmp_path):
        listener = init_audible_transits(_ais_metadata(), _ais_paths(tmp_path))
        assert isinstance(listener, AudibleTransitsAIS)
        assert listener.database_type == "AIS"
        assert listener.three_dimensional_run is True


class TestExtractAircraftInfoAis:
    def test_maps_mmsi_and_shiptype(self, tmp_path):
        listener = AudibleTransitsAIS(_ais_metadata(), _ais_paths(tmp_path))
        listener.tracks = _ais_track_points()
        listener.extract_aircraft_info()

        assert listener.tracks["n_number"].tolist() == ["368018710", "368018710"]
        assert listener.tracks["aircraft_type"].tolist() == ["Tug", "Tug"]

    def test_falls_back_to_track_id_without_mmsi_column(self, tmp_path):
        listener = AudibleTransitsAIS(_ais_metadata(), _ais_paths(tmp_path))
        tracks = _ais_track_points()
        tracks = tracks.drop(columns=["MMSI"])
        listener.tracks = tracks
        listener.extract_aircraft_info()

        assert listener.tracks["n_number"].tolist() == ["368018710", "368018710"]
        assert listener.tracks["aircraft_type"].tolist() == ["Tug", "Tug"]


class TestAisPipelineHooks:
    def test_detect_takeoffs_is_noop(self, tmp_path):
        listener = AudibleTransitsAIS(_ais_metadata(), _ais_paths(tmp_path))
        listener.tracks = _ais_track_points()
        listener.detect_takeoffs_and_landings()  # should not raise


class TestLoadTracksFromDatabaseAis:
    def test_builds_tracks_schema_from_query(self, tmp_path):
        listener = AudibleTransitsAIS(_ais_metadata(), _ais_paths(tmp_path))
        listener.study_start = "2025-01-07"
        listener.study_end = "2025-01-07"
        listener.three_dimensional_run = False
        listener.active_layer = gpd.GeoDataFrame(
            geometry=[box(-135.0, 58.3, -134.9, 58.5)], crs="EPSG:4326"
        )
        listener.active_3d = None
        listener.mic = SimpleNamespace(lat=58.40, lon=-134.97)

        mock_raw = gpd.GeoDataFrame(
            {
                "event_id": ["368018710_0_20250107", "368018710_0_20250107"],
                "TIME": [
                    pd.Timestamp("2025-01-07 12:58:04"),
                    pd.Timestamp("2025-01-07 13:08:04"),
                ],
                "altitude": [0.0, 0.0],
                "shiptype": ["Tug", "Tug"],
            },
            geometry=[Point(-134.97, 58.40), Point(-134.98, 58.41)],
            crs="EPSG:4326",
        )

        with patch(
            "nps_active_space.scripts.run_audible_transits.query_ais_mxak",
            return_value=mock_raw,
        ) as query_mock:
            loaded = listener.load_tracks_from_database(buffer=1000)

        query_mock.assert_called_once_with(
            listener.paths["AIS"],
            "2025-01-07",
            "2025-01-07",
            mask=listener.active_layer,
            mask_buffer_distance=1000,
        )
        assert list(loaded.columns) == [
            "track_id",
            "point_dt",
            "z",
            "shiptype",
            "geometry",
        ]
        assert loaded["track_id"].tolist() == [
            "368018710_0_20250107",
            "368018710_0_20250107",
        ]
        assert loaded["z"].tolist() == [0.0, 0.0]
        assert all(geom.has_z for geom in loaded.geometry)
        assert len(listener.tracks) == len(loaded)
