import pytest
from shapely.geometry import Polygon

import geopandas as gpd

from nps_active_space.utils.helpers import load_layered_activespace


class TestLoadLayeredActivespace:
    def test_finds_layer_dirs_with_posix_paths(self, tmp_path):
        project_dir = tmp_path / "project"
        site_dir = project_dir / "GLBALSTL" / "Output_Data" / "ACTIVESPACES"
        for altitude in (0, 300):
            layer_dir = site_dir / f"GLBALSTL2024_{altitude}m"
            layer_dir.mkdir(parents=True)
            gpd.GeoDataFrame(
                geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
                crs="EPSG:4326",
            ).to_file(layer_dir / "GLBALSTL2024_O_+100.geojson", driver="GeoJSON")

        study_shp = project_dir / "GLBALSTL" / "GLBALSTL_study_area.shp"
        study_shp.parent.mkdir(parents=True, exist_ok=True)
        gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])],
            crs="EPSG:4326",
        ).to_file(study_shp)

        model = load_layered_activespace(str(project_dir), "GLBA", "LSTL", 2024, gain=10.0)
        assert sorted(model.layer_dirs) == [0, 300]

    def test_raises_when_no_layers(self, tmp_path):
        project_dir = tmp_path / "project"
        (project_dir / "GLBALSTL" / "Output_Data" / "ACTIVESPACES").mkdir(parents=True)
        gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
            crs="EPSG:4326",
        ).to_file(project_dir / "GLBALSTL" / "GLBALSTL_study_area.shp")

        with pytest.raises(FileNotFoundError, match="No active space layers with GeoJSON"):
            load_layered_activespace(str(project_dir), "GLBA", "LSTL", 2024)

    def test_ignores_empty_layer_directories(self, tmp_path):
        project_dir = tmp_path / "project"
        site_dir = project_dir / "GLBALSTL" / "Output_Data" / "ACTIVESPACES"
        (site_dir / "GLBALSTL2024_0m").mkdir(parents=True)
        (site_dir / "GLBALSTL2024_300m").mkdir(parents=True)
        gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
            crs="EPSG:4326",
        ).to_file(project_dir / "GLBALSTL" / "GLBALSTL_study_area.shp")

        with pytest.raises(FileNotFoundError, match="No active space layers with GeoJSON"):
            load_layered_activespace(str(project_dir), "GLBA", "LSTL", 2024)
