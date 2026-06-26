from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from nps_active_space.utils.ais.query import query_ais_mxak
from nps_active_space.utils.ais.reader import MxakAis

REPO = Path(__file__).resolve().parents[3]
AIS_FIXTURE = REPO / "example_data" / "MXAK-AIS-GLBA" / "MXAK-AIS-GLBA-20250107.csv"
AIS_DIR = AIS_FIXTURE.parent

METADATA_COLS = ["TIME", "MMSI", "lat", "lon", "ship_name", "shiptype", "altitude", "event_id"]


def _tabular(gdf: gpd.GeoDataFrame, columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(gdf)[columns].reset_index(drop=True)


@pytest.fixture(scope="module")
def mxak_fixture_day() -> MxakAis:
    return MxakAis([AIS_FIXTURE])


class TestQueryAisMxak:
    def test_matches_direct_parse_for_fixture_day(self, mxak_fixture_day: MxakAis):
        queried = query_ais_mxak(AIS_DIR, "2025-01-07", "2025-01-07")
        assert_frame_equal(
            _tabular(queried, METADATA_COLS),
            _tabular(mxak_fixture_day, METADATA_COLS),
            check_dtype=False,
        )

    def test_date_filter_excludes_other_days(self):
        with pytest.raises(AssertionError, match="No AIS data loaded"):
            query_ais_mxak(AIS_DIR, "2025-01-08", "2025-01-08")
