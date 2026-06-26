from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from nps_active_space.utils.ais.reader import MxakAis

REPO = Path(__file__).resolve().parents[3]
AIS_FIXTURE = REPO / "example_data" / "MXAK-AIS-GLBA" / "MXAK-AIS-GLBA-20250107.csv"
ANNA_MMSI = 368018710

METADATA_COLS = ["TIME", "MMSI", "lat", "lon", "ship_name", "shiptype", "altitude", "event_id"]


def _tabular(gdf: gpd.GeoDataFrame, columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(gdf)[columns].reset_index(drop=True)


@pytest.fixture(scope="module")
def mxak_fixture_day() -> MxakAis:
    return MxakAis([AIS_FIXTURE])


class TestMxakAisParser:
    def test_parse_anna_t_first_waypoint(self, mxak_fixture_day: MxakAis):
        anna = mxak_fixture_day.loc[mxak_fixture_day["MMSI"] == ANNA_MMSI].sort_values("TIME")
        actual = _tabular(anna.head(1), METADATA_COLS)
        expected = pd.DataFrame(
            {
                "TIME": [pd.Timestamp("2025-01-07 04:58:04")],
                "MMSI": [ANNA_MMSI],
                "lat": [58.40340333],
                "lon": [-134.97689333],
                "ship_name": ["ANNA T"],
                "shiptype": ["Tug"],
                "altitude": [0.0],
                "event_id": ["368018710_0_20250107"],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_segmentation_scratch_columns_are_dropped(self, mxak_fixture_day: MxakAis):
        assert not set(mxak_fixture_day.columns).intersection({"dur_secs", "diff_event", "cumsum"})

    def test_missing_directory_raises_file_not_found(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="MXAK AIS directory not found"):
            MxakAis(tmp_path / "missing_ais_dir")

    def test_missing_timestamp_column_raises_key_error(self, tmp_path: Path):
        csv_path = tmp_path / "MXAK-AIS-TEST-20250107.csv"
        csv_path.write_text("MMSI,Latitude,Longitude\n123456789,58.0,-135.0\n")
        with pytest.raises(KeyError, match="Base station time stamp"):
            MxakAis([csv_path])

    def test_unparseable_timestamp_raises_value_error(self, tmp_path: Path):
        csv_path = tmp_path / "MXAK-AIS-TEST-20250107.csv"
        csv_path.write_text(
            "bs_ts,mmsi,callsign,imo,name,nav_status,lat,lon,cog,sog,"
            "destination,eta,shiptype,draught,length,width\n"
            "2024-11-03 01:30:00 AKST,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.40,-134.97,190.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "not-a-timestamp,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.41,-134.98,191.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "2024-11-03 01:35:00 AKST,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.42,-134.99,192.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
        )
        with pytest.raises(ValueError, match="Failed to parse 1 MXAK AIS timestamp"):
            MxakAis([csv_path])

    def test_event_id_format(self, mxak_fixture_day: MxakAis):
        anna = mxak_fixture_day.loc[mxak_fixture_day["MMSI"] == ANNA_MMSI]
        assert anna["event_id"].unique().tolist() == ["368018710_0_20250107"]

    def test_time_is_naive_utc(self, mxak_fixture_day: MxakAis):
        assert mxak_fixture_day["TIME"].dt.tz is None

    def test_column_slice_does_not_reparse_files(self, mxak_fixture_day: MxakAis):
        subset = mxak_fixture_day[["TIME", "lat", "lon"]]
        assert len(subset) == len(mxak_fixture_day)
        assert list(subset.columns) == ["TIME", "lat", "lon"]

    def test_accepts_str_csv_path(self, tmp_path: Path):
        csv_path = tmp_path / "MXAK-AIS-TEST-20250107.csv"
        csv_path.write_text(
            "bs_ts,mmsi,callsign,imo,name,nav_status,lat,lon,cog,sog,"
            "destination,eta,shiptype,draught,length,width\n"
            "2025-01-07 05:00:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.40,-134.97,190.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "2025-01-07 05:05:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.42,-134.99,192.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "2025-01-07 05:10:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.43,-135.00,193.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
        )
        str_csv_path = str(csv_path)
        ais = MxakAis([str_csv_path])
        assert len(ais) == 3

    def test_csv_path_list_rejects_non_path_entries(self, tmp_path: Path):
        csv_path = tmp_path / "MXAK-AIS-TEST-20250107.csv"
        csv_path.write_text("bs_ts,mmsi,lat,lon\n")
        with pytest.raises(TypeError, match="str or Path"):
            MxakAis([csv_path, 42])

    def test_dedups_same_second_different_position(self, tmp_path):
        csv_path = tmp_path / "MXAK-AIS-TEST-20250107.csv"
        csv_path.write_text(
            "bs_ts,mmsi,callsign,imo,name,nav_status,lat,lon,cog,sog,"
            "destination,eta,shiptype,draught,length,width\n"
            "2025-01-07 05:00:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.40,-134.97,190.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "2025-01-07 05:00:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.41,-134.98,191.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "2025-01-07 05:05:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.42,-134.99,192.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
            "2025-01-07 05:10:00,368018710,WDJ8677,9173501,ANNA T,"
            "Under way using engine,58.43,-135.00,193.0,9.9,SKAGWAY,"
            "2025-12-20 02:00:00,Tug,6,32,11\n"
        )
        ais = MxakAis([csv_path])
        actual = _tabular(ais, ["TIME", "lat", "lon"])
        expected = pd.DataFrame(
            {
                "TIME": [
                    pd.Timestamp("2025-01-07 05:00:00"),
                    pd.Timestamp("2025-01-07 05:05:00"),
                    pd.Timestamp("2025-01-07 05:10:00"),
                ],
                "lat": [58.41, 58.42, 58.43],
                "lon": [-134.98, -134.99, -135.00],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)
