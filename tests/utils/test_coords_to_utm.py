import pytest

from nps_active_space.utils.computation import coords_to_utm


class TestCoordsToUtm:
    """Unit tests for the coords_to_utm() helper."""

    def test_alaska_denali(self):
        """Denali area (~63.1N, 151.0W) should land in UTM zone 5N."""
        result = coords_to_utm(63.1, -151.0)
        assert result == "epsg:26905"

    def test_conus_colorado(self):
        """Colorado (~39.7N, 104.9W) should land in UTM zone 13N."""
        result = coords_to_utm(39.7, -104.9)
        assert result == "epsg:26913"

    def test_southern_hemisphere_australia(self):
        """Sydney area (~33.9S, 151.2E) should use the 327xx (south) series."""
        result = coords_to_utm(-33.9, 151.2)
        # UTM zone = (151.2 + 180) // 6 + 1 = 56
        assert result == "epsg:32756"

    def test_zone_boundary_east_edge(self):
        """A point right at a zone boundary (lon = -90.0) sits at
        the start of UTM zone 16.  Verify the boundary is handled."""
        result = coords_to_utm(45.0, -90.0)
        assert result == "epsg:26916"

    def test_zone_boundary_west_edge(self):
        """Longitude just west of -90 should still be in UTM zone 15."""
        result = coords_to_utm(45.0, -90.1)
        assert result == "epsg:26915"

    def test_prime_meridian(self):
        """Longitude 0 should fall in UTM zone 31."""
        result = coords_to_utm(51.5, 0.0)
        assert result == "epsg:26931"

    def test_returned_string_format_northern(self):
        """Northern hemisphere results should match the epsg:269XX pattern."""
        result = coords_to_utm(10.0, -70.0)
        assert result.startswith("epsg:269")
        zone_num = int(result.split("epsg:269")[1])
        assert 1 <= zone_num <= 60

    def test_returned_string_format_southern(self):
        """Southern hemisphere results should match the epsg:327XX pattern."""
        result = coords_to_utm(-10.0, -70.0)
        assert result.startswith("epsg:327")
        zone_num = int(result.split("epsg:327")[1])
        assert 1 <= zone_num <= 60

    def test_far_east_zone_60(self):
        """Longitude 179.0 should be in UTM zone 60."""
        result = coords_to_utm(30.0, 179.0)
        assert result == "epsg:26960"

    def test_far_west_zone_1(self):
        """Longitude -179.0 should be in UTM zone 1."""
        result = coords_to_utm(30.0, -179.0)
        assert result == "epsg:26901"
