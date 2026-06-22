import geopandas as gpd
from shapely.geometry import Point

from nps_active_space.ground_truthing import vehicle_info
from nps_active_space.utils.enums import TrackSource


def _ais_points() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "ship_name": ["ANNA T"],
            "shiptype": ["Tug"],
        },
        geometry=[Point(-134.97, 58.40)],
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

    def test_non_ais_returns_none(self):
        assert vehicle_info.lookup_vessel(TrackSource.ADSB, "x", _ais_points()) == (
            None,
            None,
            None,
        )
