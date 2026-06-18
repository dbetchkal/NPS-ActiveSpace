from pathlib import Path

import pandas as pd
import pytest

from nps_active_space.utils.ais import MxakAis
from nps_active_space.utils.time_utils import site_timezone_name, utc_naive_to_site_naive
from nps_active_space.utils.models import Nvspl

REPO = Path(__file__).resolve().parents[3]
AIS_MAY24 = REPO / "example_data" / "AIS" / "MXAK-AIS-GLBA-20240524.csv"
NVSPL_DIR = REPO / "example_data" / "NVSPL" / "GLBALSTL_2024"
SHIP_VISITS = REPO / "example_data" / "GLBALSTL_ship_visits.csv"

GLBALSTL_LAT, GLBALSTL_LON = 58.78, -136.32


class TestMay24Alignment:
    """Cross-check May 24 GLBALSTL ship visits against AIS (UTC) and NVSPL (local)."""

    @pytest.fixture
    def visits_may24(self) -> pd.DataFrame:
        visits = pd.read_csv(SHIP_VISITS, parse_dates=["CPATime"])
        return visits[visits["CPATime"].dt.date == pd.Timestamp("2024-05-24").date()]

    def test_ais_cpa_matches_ship_visits_after_site_conversion(self, visits_may24):
        ais = MxakAis([str(AIS_MAY24)])
        site_tz = site_timezone_name(GLBALSTL_LAT, GLBALSTL_LON)

        for _, row in visits_may24.iterrows():
            ship = ais[ais["MMSI"] == row.MMSI].copy()
            assert not ship.empty, f"No AIS data for MMSI {row.MMSI}"
            ship["local_time"] = utc_naive_to_site_naive(ship["TIME"], site_tz)

            target = pd.Timestamp(row.CPATime)
            idx = (ship["local_time"] - target).abs().idxmin()
            delta_sec = abs((ship.loc[idx, "local_time"] - target).total_seconds())
            assert delta_sec < 20 * 60, (
                f"{row.ShipName}: AIS local CPA {ship.loc[idx, 'local_time']} "
                f"vs visits {target} (delta {delta_sec:.0f}s)"
            )

    def test_nvspl_noise_peak_near_visit_cpa(self, visits_may24):
        for _, row in visits_may24.iterrows():
            cpa = pd.Timestamp(row.CPATime)
            hour_file = NVSPL_DIR / f"NVSPL_GLBALSTL_2024_05_24_{cpa.hour:02d}.txt"
            nvspl = Nvspl([str(hour_file)])
            window = nvspl.loc[
                cpa - pd.Timedelta(minutes=15) : cpa + pd.Timedelta(minutes=15),
                "dbA",
            ]
            hour_peak = nvspl["dbA"].max()
            window_peak = window.max()
            assert window_peak >= hour_peak - 5, (
                f"{row.ShipName}: NVSPL peak in ±15 min of CPA ({cpa}) "
                f"should be near hour max (window={window_peak:.1f}, hour={hour_peak:.1f})"
            )
