import pandas as pd

from nps_active_space.utils.time_utils import site_timezone_name, utc_naive_to_site_naive

GLBALSTL_LAT, GLBALSTL_LON = 58.79696, -136.35875


class TestSiteConversion:
    def test_glba_site_is_juneau(self):
        assert site_timezone_name(GLBALSTL_LAT, GLBALSTL_LON) == "America/Juneau"

    def test_utc_to_site_naive_may24(self):
        utc = pd.Series([pd.Timestamp("2024-05-24 15:40:00")])
        local = utc_naive_to_site_naive(utc, "America/Juneau")
        assert local.iloc[0] == pd.Timestamp("2024-05-24 07:40:00")
