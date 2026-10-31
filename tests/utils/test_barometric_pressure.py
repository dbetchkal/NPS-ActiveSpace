import numpy as np
import pytest

from nps_active_space.utils.computation import barometric_pressure


class TestBarometricPressure:
    """Reference values are the International Standard Atmosphere below the tropopause."""

    @pytest.mark.parametrize(
        "altitude_m, expected_kpa",
        [
            (0, 101.325),
            (1000, 89.875),
            (2000, 79.495),
            (5000, 54.020),
            (11000, 22.632),
        ],
    )
    def test_matches_standard_atmosphere(self, altitude_m: float, expected_kpa: float):
        assert barometric_pressure(altitude_m) == pytest.approx(expected_kpa, abs=1e-2)

    def test_decreases_with_altitude(self):
        altitudes_m = np.array([0, 500, 1000, 5000, 11000])
        pressures_kpa = np.array([barometric_pressure(h) for h in altitudes_m])
        assert np.all(np.diff(pressures_kpa) < 0)
