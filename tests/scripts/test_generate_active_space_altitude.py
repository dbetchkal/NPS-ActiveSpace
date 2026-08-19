import geopandas as gpd
import pytest
from shapely.geometry import Point

from nps_active_space.scripts.generate_active_space import infer_annotation_altitude


def _valid_points(rows: list[dict]) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")


class TestInferAnnotationAltitude:
    def test_sea_level_audible_points_yield_zero_altitude(self):
        valid_points = _valid_points(
            [
                {"annotation_idx": 0, "audible": True, "geometry": Point(0, 0, 0)},
                {"annotation_idx": 0, "audible": True, "geometry": Point(0, 0, 0)},
            ]
        )

        assert infer_annotation_altitude(valid_points) == 0

    def test_positive_audible_points_still_inferred(self):
        valid_points = _valid_points(
            [
                {"annotation_idx": 0, "audible": True, "geometry": Point(0, 0, 100)},
                {"annotation_idx": 1, "audible": True, "geometry": Point(0, 0, 200)},
            ]
        )

        assert infer_annotation_altitude(valid_points) == 150

    def test_raises_when_no_audible_altitude_points(self):
        valid_points = _valid_points(
            [
                {"annotation_idx": 0, "audible": False, "geometry": Point(0, 0, 0)},
                {"annotation_idx": 1, "audible": True, "geometry": Point(0, 0, -5)},
            ]
        )

        with pytest.raises(ValueError, match="pass -l/--altitude explicitly"):
            infer_annotation_altitude(valid_points)
