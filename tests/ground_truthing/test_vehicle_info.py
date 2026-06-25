import geopandas as gpd
import pandas as pd
from pandas.testing import assert_series_equal
from shapely.geometry import Point

from nps_active_space.ground_truthing import vehicle_info
from nps_active_space.ground_truthing.vehicle_info import TrackVehicleInfo
from nps_active_space.utils.enums import TrackSource


def _faa_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "MODE S CODE HEX": ["A98046"],
            "N-NUMBER": ["12345"],
            "TYPE AIRCRAFT": ["Fixed-wing"],
        }
    )


def _ais_points(
    *,
    ship_name: str | None = "ANNA T",
    shiptype: str | None = "Tug",
) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "ship_name": [ship_name],
            "shiptype": [shiptype],
        },
        geometry=[Point(-134.97, 58.40)],
        crs="EPSG:4326",
    )


def _adsb_points(icao: str = "A98046") -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"ICAO_address": [icao]},
        geometry=[Point(-149.0, 63.7)],
        crs="EPSG:4326",
    )


class TestLookupVessel:
    def test_ais_returns_name_and_type(self):
        help_text, vessel_type, vessel_name = vehicle_info.lookup_vessel(
            TrackSource.AIS, "368018710_0_20250107", _ais_points()
        )
        assert help_text == "MMSI: 368018710"
        assert vessel_type == "Tug"
        assert vessel_name == "ANNA T"

    def test_ais_omits_missing_name_and_type(self):
        help_text, vessel_type, vessel_name = vehicle_info.lookup_vessel(
            TrackSource.AIS, "368018710_0_20250107", _ais_points(ship_name=None, shiptype=None)
        )
        assert help_text == "MMSI: 368018710"
        assert vessel_type is None
        assert vessel_name is None

    def test_non_ais_track_source_returns_none(self):
        assert vehicle_info.lookup_vessel(TrackSource.ADSB, "x", _ais_points()) == (
            None,
            None,
            None,
        )


class TestLookupAircraft:
    def test_adsb_returns_faa_row_and_n_number_help_text(self):
        faa = _faa_table()
        faa_row, help_text, aircraft_type = vehicle_info.lookup_aircraft(
            faa, TrackSource.ADSB, "A98046_0_20250623", _adsb_points()
        )
        assert_series_equal(faa_row, faa.iloc[0])
        assert help_text == "N-Number: 12345"
        assert aircraft_type == "Fixed-wing"

    def test_gps_returns_faa_row_and_icao_help_text(self):
        faa = _faa_table()
        points = gpd.GeoDataFrame(geometry=[Point(0, 0)], crs="EPSG:4326")
        faa_row, help_text, aircraft_type = vehicle_info.lookup_aircraft(
            faa, TrackSource.GPS, "N12345_0_20250623", points
        )
        assert_series_equal(faa_row, faa.iloc[0])
        assert help_text == "ICAO Address: A98046"
        assert aircraft_type == "Fixed-wing"

    def test_none_faa_returns_nones(self):
        assert vehicle_info.lookup_aircraft(
            None, TrackSource.ADSB, "A98046_0_20250623", _adsb_points()
        ) == (None, None, None)

    def test_adsb_no_match_returns_nones(self):
        faa = _faa_table()
        missing_icao_points = _adsb_points(icao="DEADBE")
        assert vehicle_info.lookup_aircraft(
            faa, TrackSource.ADSB, "A98046_0_20250623", missing_icao_points
        ) == (None, None, None)

    def test_gps_no_match_returns_nones(self):
        faa = _faa_table()
        missing_n_number_track_id = "N99999_0_20250623"
        points = gpd.GeoDataFrame(geometry=[Point(0, 0)], crs="EPSG:4326")
        assert vehicle_info.lookup_aircraft(
            faa, TrackSource.GPS, missing_n_number_track_id, points
        ) == (None, None, None)

    def test_ais_returns_nones(self):
        faa = _faa_table()
        assert vehicle_info.lookup_aircraft(
            faa, TrackSource.AIS, "368018710_0_20250107", _ais_points()
        ) == (None, None, None)


class TestLookupTrackVehicle:
    def test_ais_returns_vessel_metadata(self):
        info = vehicle_info.lookup_track_vehicle(
            None, TrackSource.AIS, "368018710_0_20250107", _ais_points()
        )
        assert info == TrackVehicleInfo(None, "MMSI: 368018710", "Tug", "ANNA T")

    def test_adsb_returns_aircraft_metadata(self):
        faa = _faa_table()
        faa_row = faa.iloc[0]
        info = vehicle_info.lookup_track_vehicle(
            faa, TrackSource.ADSB, "A98046_0_20250623", _adsb_points()
        )
        expected = TrackVehicleInfo(faa_row, "N-Number: 12345", "Fixed-wing", None)
        assert_series_equal(info.faa_row, faa_row)
        assert (info.help_text, info.vehicle_type, info.vessel_name) == (
            expected.help_text,
            expected.vehicle_type,
            expected.vessel_name,
        )

    def test_gps_returns_aircraft_metadata(self):
        faa = _faa_table()
        faa_row = faa.iloc[0]
        points = gpd.GeoDataFrame(geometry=[Point(0, 0)], crs="EPSG:4326")
        info = vehicle_info.lookup_track_vehicle(
            faa, TrackSource.GPS, "N12345_0_20250623", points
        )
        expected = TrackVehicleInfo(faa_row, "ICAO Address: A98046", "Fixed-wing", None)
        assert_series_equal(info.faa_row, faa_row)
        assert (info.help_text, info.vehicle_type, info.vessel_name) == (
            expected.help_text,
            expected.vehicle_type,
            expected.vessel_name,
        )

    def test_missing_faa_returns_empty_track_vehicle_info(self):
        info = vehicle_info.lookup_track_vehicle(
            None, TrackSource.ADSB, "A98046_0_20250623", _adsb_points()
        )
        assert info == TrackVehicleInfo(None, None, None, None)
