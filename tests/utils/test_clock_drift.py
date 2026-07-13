import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

from nps_active_space.utils.clock_drift import correct_clock_drift


def _make_tracks(times):
    return gpd.GeoDataFrame(
        {"point_dt": pd.to_datetime(times)},
        geometry=[Point(0, 0) for _ in times],
    )


class TestCorrectClockDrift:

    def test_inplace_false_returns_corrected_copy(self, tmp_path):
        # drift is 10s at 00:00 and 20s at 01:00, so a point at 00:30 should
        # pick up the interpolated 15s.
        drift_file = tmp_path / "drift.csv"
        drift_file.write_text(
            "Time,Seconds\n"
            "2020-01-01 00:00:00,10\n"
            "2020-01-01 01:00:00,20\n"
        )

        tracks = _make_tracks(["2020-01-01 00:00:00", "2020-01-01 00:30:00"])
        original = tracks["point_dt"].tolist()

        result = correct_clock_drift(tracks, str(drift_file), inplace=False)

        expected = pd.to_datetime(["2020-01-01 00:00:10", "2020-01-01 00:30:15"])
        assert result["point_dt"].tolist() == list(expected)
        # inplace=False shouldn't touch the caller's frame
        assert tracks["point_dt"].tolist() == original
