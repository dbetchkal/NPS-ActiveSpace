from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from nps_active_space.utils.helpers import query_adsb
from nps_active_space.utils.models import Adsb

ADSB_FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "example_data"
    / "ADS-B"
    / "healy_repeater"
    / "20250623.TSV"
)
ADSB_DIR = ADSB_FIXTURE.parent
FIXTURE_DAY_START = pd.Timestamp("2025-06-23")
FIXTURE_DAY_END = pd.Timestamp("2025-06-24")

PARSE_CHECK_COLS = ["lat", "lon", "altitude", "flight_id", "TIME"]


def _tabular(gdf: gpd.GeoDataFrame, columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(gdf)[columns].reset_index(drop=True)


@pytest.fixture(scope="module")
def adsb_fixture_day() -> Adsb:
    return Adsb(
        [str(ADSB_FIXTURE)],
        start_dt=FIXTURE_DAY_START,
        end_dt=FIXTURE_DAY_END,
    )


class TestAdsbParser:
    def test_parse_scales_coordinates_and_flight_id(self, adsb_fixture_day: Adsb):
        adsb = adsb_fixture_day
        assert not adsb.empty
        assert "ICAO_address" in adsb.columns
        assert "flight_id" in adsb.columns

        cfs = adsb[adsb["ICAO_address"] == "A98046"]
        assert not cfs.empty
        actual = _tabular(cfs.head(1), PARSE_CHECK_COLS)
        expected = pd.DataFrame(
            {
                "lat": [63.7391488],
                "lon": [-149.02944],
                "altitude": [5791.2],
                "flight_id": ["A98046_0_20250623"],
                "TIME": [pd.Timestamp("2025-06-23 12:53:23")],
            }
        )
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_flight_id_format(self, adsb_fixture_day: Adsb):
        for flight_id in adsb_fixture_day["flight_id"].unique():
            icao, cumsum, date = flight_id.split("_")
            assert len(icao) == 6
            assert cumsum.isdigit()
            assert date == "20250623"

    def test_drops_pre_2019_timestamps(self, adsb_fixture_day: Adsb):
        assert adsb_fixture_day["TIME"].min() >= pd.Timestamp("2019-01-01")
        assert adsb_fixture_day["TIME"].min() >= FIXTURE_DAY_START


class TestQueryAdsb:
    def test_date_filter_inclusive_single_day(self):
        adsb = query_adsb(str(ADSB_DIR), "2025-06-23", "2025-06-23")
        assert not adsb.empty
        assert adsb["TIME"].min() >= pd.Timestamp("2025-06-23")
        assert adsb["TIME"].max() < pd.Timestamp("2025-06-24")

    def test_date_filter_excludes_other_days(self):
        with pytest.raises(AssertionError, match="No ADSB data loaded"):
            query_adsb(str(ADSB_DIR), "2025-01-07", "2025-01-08")
