import numpy as np
import pytest

from nps_active_space.utils.computation import climb_angle


class TestClimbAngle:

    def test_non_unit_vectors(self):
        # 3-4-5 triangle, so the climb angle is asin(4/5)
        assert climb_angle(np.array([3., 0., 4.])) == pytest.approx(53.130102, abs=1e-5)
        assert climb_angle(np.array([10., 0., 10.])) == pytest.approx(45., abs=1e-5)
        assert climb_angle(np.array([1., 1., 1.])) == pytest.approx(35.264390, abs=1e-5)

    def test_unit_vectors(self):
        assert climb_angle(np.array([0., 0., 1.])) == pytest.approx(90.)
        assert climb_angle(np.array([1., 0., 0.])) == pytest.approx(0.)
        assert climb_angle(np.array([0., 0., -1.])) == pytest.approx(-90.)
