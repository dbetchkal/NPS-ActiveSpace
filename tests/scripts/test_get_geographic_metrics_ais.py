import pandas as pd

from nps_active_space.scripts.get_geographic_metrics import filter_geo_metric_tracks


class TestFilterGeoMetricTracks:
    def _sample_tracks(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "track_id": ["a", "b", "c"],
                "aircraft_type": ["Fixed-wing", "Jet", "Tug"],
            }
        )

    def test_aircraft_keeps_fixed_wing_only(self):
        tracks = self._sample_tracks()
        filtered = filter_geo_metric_tracks(tracks, "ADSB")
        assert filtered["track_id"].tolist() == ["a"]

    def test_ais_keeps_all_vessel_types(self):
        tracks = self._sample_tracks()
        filtered = filter_geo_metric_tracks(tracks, "AIS")
        assert filtered["track_id"].tolist() == ["a", "b", "c"]
