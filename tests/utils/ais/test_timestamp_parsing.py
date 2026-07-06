import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from nps_active_space.utils.ais.timestamp_parsing import parse_mxak_ais_timestamps


class TestParseMxakAisTimestamps:
    def test_bare_format_is_utc(self):
        time = pd.Series(["2024-05-24 15:40:00", "2024-05-24 15:41:00"])
        expected = pd.Series(
            [
                pd.Timestamp("2024-05-24 15:40:00"),
                pd.Timestamp("2024-05-24 15:41:00"),
            ]
        )
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected)
        assert out.dt.tz is None

    def test_utc_pattern1_day_month_year(self):
        """Pattern 1: ``%d %b %Y %H:%M:%S UTC``."""
        time = pd.Series(
            [
                "24 May 2024 15:40:00 UTC",
                "24 May 2024 15:41:00 UTC",
            ]
        )
        expected = pd.Series(
            [
                pd.Timestamp("2024-05-24 15:40:00"),
                pd.Timestamp("2024-05-24 15:41:00"),
            ]
        )
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected)

    def test_utc_pattern2_iso_date(self):
        """Pattern 2: ``%Y-%m-%d %H:%M:%S UTC`` (single format per file)."""
        time = pd.Series(
            [
                "2024-05-24 15:40:00 UTC",
                "2024-05-24 15:41:00 UTC",
            ]
        )
        expected = pd.Series(
            [
                pd.Timestamp("2024-05-24 15:40:00"),
                pd.Timestamp("2024-05-24 15:41:00"),
            ]
        )
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected)

    def test_utc_pattern3_us_date(self):
        """Pattern 3: ``%m-%d-%Y %H:%M:%S UTC``."""
        time = pd.Series(
            [
                "05-24-2024 15:40:00 UTC",
                "05-24-2024 15:41:00 UTC",
            ]
        )
        expected = pd.Series(
            [
                pd.Timestamp("2024-05-24 15:40:00"),
                pd.Timestamp("2024-05-24 15:41:00"),
            ]
        )
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected)

    def test_akdt_suffix_converts_to_utc(self):
        time = pd.Series(["2024-05-24 07:40:00 AKDT"])
        expected = pd.Series([pd.Timestamp("2024-05-24 15:40:00")])
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected)

    def test_akst_suffix_converts_to_utc(self):
        time = pd.Series(["2024-01-15 07:40:00 AKST"])
        expected = pd.Series([pd.Timestamp("2024-01-15 16:40:00")])
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected)

    def test_mixed_akst_akdt_november_dst_fallback(self):
        """November fall-back day: mixed AKST/AKDT rows trigger DST fallback."""
        # First row is AKST (primary format); AKDT row fails primary parse.
        time = pd.Series(
            [
                "2024-11-03 01:30:00 AKST",
                "2024-11-03 01:30:00 AKDT",
            ]
        )
        expected = pd.Series(
            [
                pd.Timestamp("2024-11-03 10:30:00"),  # 01:30 AKST + 9 h
                pd.Timestamp("2024-11-03 09:30:00"),  # 01:30 AKDT + 8 h
            ]
        )
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected, check_dtype=False)

    def test_mixed_akst_akdt_march_dst_fallback(self):
        """March spring-forward day: mixed AKDT/AKST rows trigger DST fallback."""
        # First row is AKDT (primary format); AKST row fails primary parse.
        time = pd.Series(
            [
                "2024-03-10 03:30:00 AKDT",
                "2024-03-10 01:30:00 AKST",
            ]
        )
        expected = pd.Series(
            [
                pd.Timestamp("2024-03-10 11:30:00"),  # 03:30 AKDT + 8 h
                pd.Timestamp("2024-03-10 10:30:00"),  # 01:30 AKST + 9 h
            ]
        )
        out = parse_mxak_ais_timestamps(time)
        assert_series_equal(out, expected, check_dtype=False)

    def test_unrecognized_format_raises_value_error(self):
        time = pd.Series(["not-a-timestamp"])
        with pytest.raises(ValueError, match="not recognized"):
            parse_mxak_ais_timestamps(time)

    def test_mixed_akst_with_unparseable_row_raises(self):
        """DST fallback only re-parses AKST/AKDT rows; other suffixes must not become NaT."""
        time = pd.Series(
            [
                "2024-11-03 01:30:00 AKST",
                "2024-11-03 01:30:00 UTC",
            ]
        )
        with pytest.raises(ValueError, match="Failed to parse 1 MXAK AIS timestamp"):
            parse_mxak_ais_timestamps(time)
