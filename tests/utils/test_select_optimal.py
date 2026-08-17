import geopandas as gpd
from shapely.geometry import Point, Polygon

from nps_active_space.utils.computation import select_optimal


class TestSelectOptimal:
    def test_returns_best_omni_when_all_fbeta_scores_are_zero(self):
        valid_points = gpd.GeoDataFrame(
            {"audible": [True, False]},
            geometry=[Point(0, 0), Point(1, 1)],
            crs="EPSG:4326",
        )
        distant_active_space = gpd.GeoDataFrame(
            geometry=[Polygon([(10, 10), (11, 10), (11, 11), (10, 11)])],
            crs="EPSG:4326",
        )
        active_space_polygons = [
            ("O_+100", distant_active_space),
            ("O_+200", distant_active_space),
        ]

        best_omni, max_fbeta, _, _, detection_results = select_optimal(
            unit="DENA",
            site="TRLA",
            year=2025,
            valid_points=valid_points,
            active_space_polygons=active_space_polygons,
            verbose=False,
            plot=False,
        )

        assert best_omni is not None
        assert max_fbeta == 0.0
        assert len(detection_results) == 2
