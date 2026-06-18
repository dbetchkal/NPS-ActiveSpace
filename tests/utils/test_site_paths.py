import pytest

from nps_active_space.utils.helpers import get_deployment, load_DEM, load_studyarea


class TestSitePathErrors:
    def test_load_dem_missing_file(self, tmp_path):
        (tmp_path / "GLBABART" / "Input_Data" / "01_ELEVATION").mkdir(parents=True)
        with pytest.raises(FileNotFoundError, match="elevation DEM for GLBABART"):
            load_DEM(str(tmp_path), "GLBA", "BART")

    def test_load_studyarea_missing_file(self, tmp_path):
        (tmp_path / "GLBABART").mkdir(parents=True)
        with pytest.raises(FileNotFoundError, match="study area shapefile for GLBABART2025"):
            load_studyarea(str(tmp_path), "GLBA", "BART", 2025)

    def test_get_deployment_missing_sit(self, tmp_path):
        site_dir = tmp_path / "GLBABART" / "Input_Data" / "05_SITES"
        site_dir.mkdir(parents=True)
        with pytest.raises(FileNotFoundError, match="Microphone site file not found"):
            get_deployment(str(tmp_path), "GLBA", "BART", 2025, elevation=False)

    def test_get_deployment_invalid_sit(self, tmp_path):
        sit_path = (
            tmp_path / "GLBABART" / "Input_Data" / "05_SITES" / "GLBABART2025.sit"
        )
        sit_path.parent.mkdir(parents=True)
        sit_path.write_text("too\nshort")
        with pytest.raises(ValueError, match="invalid format"):
            get_deployment(str(tmp_path), "GLBA", "BART", 2025, elevation=False)
