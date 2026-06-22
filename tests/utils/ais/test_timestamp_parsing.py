import pandas as pd

from nps_active_space.utils.ais.timestamp_parsing import parse_mxak_ais_timestamps


class TestParseMxakAisTimestamps:
    def test_bare_format_is_utc(self):
        time = pd.Series(["2024-05-24 15:40:00", "2024-05-24 15:41:00"])
        out = parse_mxak_ais_timestamps(time)
        assert out.iloc[0] == pd.Timestamp("2024-05-24 15:40:00")
        assert out.dt.tz is None

    def test_akdt_suffix_converts_to_utc(self):
        time = pd.Series(["2024-05-24 07:40:00 AKDT"])
        out = parse_mxak_ais_timestamps(time)
        assert out.iloc[0] == pd.Timestamp("2024-05-24 15:40:00")

    def test_akst_suffix_converts_to_utc(self):
        time = pd.Series(["2024-01-15 07:40:00 AKST"])
        out = parse_mxak_ais_timestamps(time)
        assert out.iloc[0] == pd.Timestamp("2024-01-15 16:40:00")
