"""Tests for 3D layered active-space fitting."""

from pathlib import Path

import geopandas as gpd
from shapely.geometry import Polygon

from nps_active_space.active_space.layered_active_space import LayeredActiveSpace

_MINI_GEOJSON = """{
  "type": "FeatureCollection",
  "features": [{
    "type": "Feature",
    "properties": {},
    "geometry": {
      "type": "Polygon",
      "coordinates": [[
        [-149.5, 63.6], [-149.4, 63.6], [-149.4, 63.7], [-149.5, 63.7], [-149.5, 63.6]
      ]]
    }
  }]
}"""


def _write_layer_gain_geojson(layer_dir: Path, usy: str, stem: str) -> None:
    layer_dir.mkdir(parents=True, exist_ok=True)
    (layer_dir / f"{usy}_{stem}.geojson").write_text(_MINI_GEOJSON)


class TestLayeredActiveSpaceGainDiscovery:
    def test_discovers_five_db_ladder_not_half_db(self, tmp_path: Path) -> None:
        usy = "DENASUSH2026"
        study_area = gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
            crs="EPSG:4326",
        )
        layer_dirs = {
            1200: str(tmp_path / f"{usy}_1200m"),
            1500: str(tmp_path / f"{usy}_1500m"),
        }
        stems = ["O_-100", "O_-050", "O_+000", "O_+050"]
        for layer in layer_dirs.values():
            for stem in stems:
                _write_layer_gain_geojson(Path(layer), usy, stem)

        layered = LayeredActiveSpace(usy, layer_dirs, study_area)
        assert layered.gain_values == [-10.0, -5.0, 0.0, 5.0]
        assert layered._resolve_fit_gain_values(-10.0, 40.0) == [-10.0, -5.0, 0.0, 5.0]

    def test_fit_gain_values_skip_missing_half_db_steps(self, tmp_path: Path) -> None:
        usy = "DENASUSH2026"
        study_area = gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
            crs="EPSG:4326",
        )
        layer_dir = str(tmp_path / f"{usy}_1200m")
        _write_layer_gain_geojson(Path(layer_dir), usy, "O_-100")
        _write_layer_gain_geojson(Path(layer_dir), usy, "O_+000")

        layered = LayeredActiveSpace(usy, {1200: layer_dir}, study_area)
        assert -9.5 not in layered.gain_values
        assert layered._resolve_fit_gain_values(-10.0, 0.0) == [-10.0, 0.0]
