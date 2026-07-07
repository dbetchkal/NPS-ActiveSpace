import numpy as np
import numpy.testing as npt

from nps_active_space.utils.computation import (
    barometric_pressure,
    contiguous_regions,
    coords_to_utm,
)


class TestContiguousRegions:
    def test_single_block(self):
        cond = np.array([False, True, True, True, False])
        result = contiguous_regions(cond)
        npt.assert_array_equal(result, [[1, 4]])

    def test_two_blocks(self):
        cond = np.array([True, True, False, False, True, True, True])
        result = contiguous_regions(cond)
        npt.assert_array_equal(result, [[0, 2], [4, 7]])

    def test_all_true(self):
        cond = np.array([True, True, True])
        result = contiguous_regions(cond)
        npt.assert_array_equal(result, [[0, 3]])

    def test_all_false(self):
        cond = np.array([False, False, False])
        result = contiguous_regions(cond)
        assert result.shape == (0, 2)

    def test_single_element_true(self):
        cond = np.array([True])
        result = contiguous_regions(cond)
        npt.assert_array_equal(result, [[0, 1]])

    def test_alternating(self):
        cond = np.array([True, False, True, False, True])
        result = contiguous_regions(cond)
        npt.assert_array_equal(result, [[0, 1], [2, 3], [4, 5]])

    def test_numeric_condition(self):
        arr = np.array([0, 5, 6, 2, 0, 8])
        result = contiguous_regions(arr > 3)
        npt.assert_array_equal(result, [[1, 3], [5, 6]])


class TestCoordsToUtm:
    def test_denali_alaska(self):
        result = coords_to_utm(63.07, -151.01)
        assert result == "epsg:26905"

    def test_new_york(self):
        result = coords_to_utm(40.71, -74.01)
        assert result == "epsg:26918"

    def test_southern_hemisphere(self):
        result = coords_to_utm(-33.87, 151.21)
        assert result == "epsg:32756"

    def test_prime_meridian(self):
        result = coords_to_utm(51.48, -0.01)
        assert result == "epsg:26930"

    def test_zone_boundary(self):
        result = coords_to_utm(45.0, -126.0)
        assert result == "epsg:26910"


class TestBarometricPressure:
    def test_sea_level(self):
        p = barometric_pressure(0)
        npt.assert_almost_equal(p, 101.325, decimal=2)

    def test_decreases_with_altitude(self):
        p0 = barometric_pressure(0)
        p1000 = barometric_pressure(1000)
        p5000 = barometric_pressure(5000)
        assert p0 > p1000 > p5000

    def test_known_altitude(self):
        p = barometric_pressure(1500)
        assert 80 < p < 90
