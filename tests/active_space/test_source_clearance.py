import numpy as np
import pytest
from shapely.geometry import Point

from nps_active_space.active_space.source_clearance import (
    SOURCE_SURFACE_CLEARANCE_M,
    apply_surface_clearance_m,
    with_point_z,
)


class TestApplySurfaceClearance:
    def test_sea_level_on_flat_water_is_lifted(self):
        z_m, underground = apply_surface_clearance_m(
            np.array([0.0]),
            np.array([0.0]),
        )
        assert not underground[0]
        assert z_m[0] == SOURCE_SURFACE_CLEARANCE_M

    def test_bathymetry_below_requested_msl_keeps_msl(self):
        z_m, underground = apply_surface_clearance_m(
            np.array([0.0]),
            np.array([-4.0]),
        )
        assert not underground[0]
        assert z_m[0] == 0.0

    def test_land_well_above_sea_level_stays_underground(self):
        z_m, underground = apply_surface_clearance_m(
            np.array([0.0]),
            np.array([50.0]),
        )
        assert underground[0]
        assert z_m[0] == 0.0

    def test_dem_noise_just_below_surface_is_lifted(self):
        z_m, underground = apply_surface_clearance_m(
            np.array([0.0]),
            np.array([0.4]),
        )
        assert not underground[0]
        assert z_m[0] == pytest.approx(0.4 + SOURCE_SURFACE_CLEARANCE_M)

    def test_airborne_msl_unchanged(self):
        z_m, underground = apply_surface_clearance_m(
            np.array([1000.0]),
            np.array([50.0]),
        )
        assert not underground[0]
        assert z_m[0] == 1000.0

    def test_missing_surface_is_underground(self):
        _z_m, underground = apply_surface_clearance_m(
            np.array([0.0]),
            np.array([np.nan]),
        )
        assert underground[0]


class TestWithPointZ:
    def test_replaces_z_only(self):
        import geopandas as gpd

        pts = gpd.GeoDataFrame(
            geometry=[Point(1.0, 2.0, 0.0)],
            crs="EPSG:32606",
        )
        out = with_point_z(pts, np.array([2.0]))
        assert out.geometry.iloc[0].x == 1.0
        assert out.geometry.iloc[0].y == 2.0
        assert out.geometry.iloc[0].z == 2.0
