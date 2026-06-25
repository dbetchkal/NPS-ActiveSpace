"""Clock alignment between MXAK AIS (UTC), ship-visit CPA (site-local), and NVSPL (site-local).

These tests do not prove acoustic detection — they guard timezone handling and loose
co-timing between bundled example files. AIS transit segmentation (``event_id``) is
tested indirectly via :class:`~nps_active_space.utils.ais.reader.MxakAis`; clock
checks use :func:`~nps_active_space.utils.ais.reader.parse_mxak_ais_csv_points` so they
stay stable when ``EVENT_GAP_SECONDS`` changes.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from nps_active_space.utils.ais.reader import EVENT_GAP_SECONDS, parse_mxak_ais_csv_points
from nps_active_space.utils.models import Nvspl
from nps_active_space.utils.time_utils import site_timezone_name, utc_naive_to_site_naive

REPO = Path(__file__).resolve().parents[3]
AIS_MAY24 = REPO / "example_data" / "AIS" / "MXAK-AIS-GLBA-20240524.csv"
NVSPL_DIR = REPO / "example_data" / "NVSPL" / "GLBALSTL_2024"
SHIP_VISITS = REPO / "example_data" / "GLBALSTL_ship_visits.csv"

GLBALSTL_LAT, GLBALSTL_LON = 58.78, -136.32
SITE_TZ = "America/Juneau"

# Same transit scale as MXAK event splitting; used for CPA windows and loose AIS timing.
CPA_WINDOW_MINUTES = EVENT_GAP_SECONDS // 60
AIS_CPA_TOLERANCE_SEC = 2 * EVENT_GAP_SECONDS
UTC_MISALIGNMENT_FLOOR_SEC = 60 * 60
NVSPL_PEAK_TOLERANCE_DB = 5

MAY24_MMSIS = (367365630, 246648000, 367578110)


def nearest_ais_local_to_cpa(
    ais_times_utc: pd.Series,
    cpa_local: pd.Timestamp,
    site_tz: str,
) -> tuple[pd.Timestamp, float]:
    """Nearest AIS site-local timestamp to a ship-visit CPA and |delta| in seconds."""
    local = utc_naive_to_site_naive(ais_times_utc, site_tz)
    idx = (local - cpa_local).abs().idxmin()
    nearest = local.loc[idx]
    delta_sec = abs((nearest - cpa_local).total_seconds())
    return nearest, delta_sec


def nvspl_dba_peaks_near_cpa(
    nvspl: pd.DataFrame,
    cpa_local: pd.Timestamp,
    window_minutes: int = CPA_WINDOW_MINUTES,
) -> tuple[float, float]:
    """Peak dBA in a ±window around CPA and peak dBA over the loaded NVSPL series."""
    window = nvspl.loc[
        cpa_local - pd.Timedelta(minutes=window_minutes) : cpa_local + pd.Timedelta(minutes=window_minutes),
        "dbA",
    ]
    return float(window.max()), float(nvspl["dbA"].max())


def _minimal_nvspl_df(index: pd.DatetimeIndex, dba: pd.Series) -> pd.DataFrame:
    """Minimal NVSPL-shaped frame (only dbA varies) for unit tests."""
    data = pd.DataFrame(index=index)
    data["SiteID"] = "TEST"
    data["dbA"] = dba.reindex(index).astype(float)
    data["dbC"] = 0.0
    data["Voltage"] = 0.0
    data["WindSpeed"] = 0.0
    data["WindDir"] = 0.0
    data["TempIns"] = 0.0
    data["TempOut"] = 0.0
    data["Humidity"] = 0.0
    data["INVID"] = ""
    data["INSID"] = ""
    data["GChar1"] = ""
    data["AdjustmentsApplied"] = ""
    data["CalibrationAdjustment"] = 0.0
    data["GPSTimeAdjustment"] = 0.0
    data["GainAdjustment"] = 0.0
    data["Status"] = ""
    return data


class TestUtcToSiteLocalConversion:
    """Synthetic checks that AIS UTC timestamps convert to site-local CPA time."""

    def test_nearest_ais_point_matches_local_cpa(self):
        cpa_local = pd.Timestamp("2024-05-24 07:40:00")
        cpa_utc = pd.Timestamp("2024-05-24 15:40:00")  # 07:40 AKDT
        ais_times_utc = pd.Series(
            [
                cpa_utc - pd.Timedelta(minutes=30),
                cpa_utc - pd.Timedelta(minutes=10),
                cpa_utc,
                cpa_utc + pd.Timedelta(minutes=10),
                cpa_utc + pd.Timedelta(minutes=30),
            ]
        )

        nearest, delta_sec = nearest_ais_local_to_cpa(ais_times_utc, cpa_local, SITE_TZ)

        assert nearest == cpa_local
        assert delta_sec == 0.0

    def test_utc_ais_without_conversion_is_hours_off_cpa(self):
        """Guards against comparing raw UTC AIS to site-local ship-visit CPA."""
        cpa_local = pd.Timestamp("2024-05-24 07:40:00")
        cpa_utc = pd.Timestamp("2024-05-24 15:40:00")
        ais_times_utc = pd.Series([cpa_utc])

        wrong_delta_sec = abs((ais_times_utc.iloc[0] - cpa_local).total_seconds())
        _, correct_delta_sec = nearest_ais_local_to_cpa(ais_times_utc, cpa_local, SITE_TZ)

        assert wrong_delta_sec == 8 * 3600
        assert correct_delta_sec == 0.0


class TestNvsplCpaWindow:
    """Synthetic NVSPL CPA-window selection (site-local timestamps, no files)."""

    def test_cpa_window_captures_injected_noise_peak(self):
        cpa = pd.Timestamp("2024-05-24 07:40:00")
        index = pd.date_range(cpa - pd.Timedelta(hours=1), cpa + pd.Timedelta(hours=1), freq="1min")
        dba = pd.Series(50.0, index=index)
        dba.loc[cpa] = 90.0

        nvspl = _minimal_nvspl_df(index, dba)
        window_peak, series_peak = nvspl_dba_peaks_near_cpa(nvspl, cpa)

        assert window_peak == 90.0
        assert series_peak == 90.0

    def test_noise_peak_outside_window_is_not_required_in_cpa_window(self):
        cpa = pd.Timestamp("2024-05-24 07:40:00")
        index = pd.date_range(cpa - pd.Timedelta(hours=1), cpa + pd.Timedelta(hours=1), freq="1min")
        dba = pd.Series(55.0, index=index)
        dba.loc[cpa - pd.Timedelta(minutes=30)] = 90.0  # peak outside ±EVENT_GAP window

        nvspl = _minimal_nvspl_df(index, dba)
        window_peak, series_peak = nvspl_dba_peaks_near_cpa(nvspl, cpa)

        assert window_peak == 55.0
        assert series_peak == 90.0
        assert window_peak < series_peak - NVSPL_PEAK_TOLERANCE_DB


class TestMay24ExampleDataClockAlignment:
    """Example-data regression: MXAK UTC → site-local must co-time with ship-visit CPA.

    Uses parsed AIS points only (no ``event_id`` filter) so results do not depend on
    transit-segmentation thresholds.
    """

    @pytest.fixture
    def visits_may24(self) -> pd.DataFrame:
        visits = pd.read_csv(SHIP_VISITS, parse_dates=["CPATime"])
        return visits[visits["CPATime"].dt.date == pd.Timestamp("2024-05-24").date()]

    @pytest.fixture(scope="class")
    def ais_points_may24(self) -> pd.DataFrame:
        return parse_mxak_ais_csv_points(AIS_MAY24)

    @pytest.mark.parametrize("mmsi", MAY24_MMSIS)
    def test_ais_cpa_matches_ship_visits_after_site_conversion(
        self,
        visits_may24: pd.DataFrame,
        ais_points_may24: pd.DataFrame,
        mmsi: int,
    ):
        row = visits_may24.loc[visits_may24["MMSI"] == mmsi].iloc[0]
        site_tz = site_timezone_name(GLBALSTL_LAT, GLBALSTL_LON)
        cpa_local = pd.Timestamp(row.CPATime)
        ship = ais_points_may24.loc[ais_points_may24["MMSI"] == mmsi]
        assert not ship.empty, f"No AIS data for MMSI {mmsi}"

        nearest, converted_delta_sec = nearest_ais_local_to_cpa(ship["TIME"], cpa_local, site_tz)
        local = utc_naive_to_site_naive(ship["TIME"], site_tz)
        idx = (local - nearest).abs().idxmin()
        raw_delta_sec = abs((ship["TIME"].loc[idx] - cpa_local).total_seconds())

        assert raw_delta_sec > UTC_MISALIGNMENT_FLOOR_SEC, (
            f"{row.ShipName}: raw UTC AIS is not hours off CPA — conversion may be skipped"
        )
        assert converted_delta_sec < AIS_CPA_TOLERANCE_SEC, (
            f"{row.ShipName}: site-local AIS {nearest} vs visit CPA {cpa_local} "
            f"(delta {converted_delta_sec:.0f}s, limit {AIS_CPA_TOLERANCE_SEC}s)"
        )


class TestMay24ExampleDataSplSanity:
    """Loose sanity check: NVSPL dBA peaks near ship-visit CPA on bundled May 24 data."""

    @pytest.fixture
    def visits_may24(self) -> pd.DataFrame:
        visits = pd.read_csv(SHIP_VISITS, parse_dates=["CPATime"])
        return visits[visits["CPATime"].dt.date == pd.Timestamp("2024-05-24").date()]

    @pytest.mark.parametrize("mmsi", MAY24_MMSIS)
    def test_nvspl_noise_peak_near_visit_cpa(
        self,
        visits_may24: pd.DataFrame,
        mmsi: int,
    ):
        row = visits_may24.loc[visits_may24["MMSI"] == mmsi].iloc[0]
        cpa = pd.Timestamp(row.CPATime)
        hour_file = NVSPL_DIR / f"NVSPL_GLBALSTL_2024_05_24_{cpa.hour:02d}.txt"
        nvspl = Nvspl([str(hour_file)])
        window_peak, hour_peak = nvspl_dba_peaks_near_cpa(nvspl, cpa)
        assert window_peak >= hour_peak - NVSPL_PEAK_TOLERANCE_DB, (
            f"{row.ShipName}: NVSPL peak in ±{CPA_WINDOW_MINUTES} min of CPA ({cpa}) "
            f"should be near hour max (window={window_peak:.1f}, hour={hour_peak:.1f})"
        )
