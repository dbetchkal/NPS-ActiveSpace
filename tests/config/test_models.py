from pathlib import Path, PureWindowsPath

import pytest

from nps_active_space.config.models import (
    ActiveSpaceConfig,
    DatabaseConfig,
    FileDataConfig,
    ProjectConfig,
)
from tests.config.helpers import make_config, toml_str


class TestConfigFieldParsing:
    @pytest.fixture
    def paths(self, tmp_path: Path) -> dict[str, Path]:
        data_dir = tmp_path / "data"
        tools_dir = tmp_path / "tools"
        data_dir.mkdir()
        tools_dir.mkdir()
        return {
            "site_metadata": data_dir / "metadata.csv",
            "nvspl_archive": data_dir / "nvspl",
            "adsb": data_dir / "adsb",
            "dem": data_dir / "dem.tif",
            "mennitt": data_dir / "mennitt.tif",
            "project_dir": tmp_path / "projects",
            "nmsim_binary": tools_dir / "NMSIM.exe",
            "faa_db": data_dir / "MASTER.txt",
        }

    @pytest.fixture
    def full_config(
        self, tmp_config_dir: Path, paths: dict[str, Path]
    ) -> tuple[ActiveSpaceConfig, dict[str, Path]]:
        toml = f'''
[database]
name     = 'my_db'
username = 'user'
password = 'pass'
port     = '5432'
host     = 'localhost'

[data]
site_metadata = {toml_str(paths["site_metadata"])}
nvspl_archive = {toml_str(paths["nvspl_archive"])}
adsb          = {toml_str(paths["adsb"])}
dem           = {toml_str(paths["dem"])}
mennitt       = {toml_str(paths["mennitt"])}

[project]
dir               = {toml_str(paths["project_dir"])}
nmsim_binary      = {toml_str(paths["nmsim_binary"])}
FAA_Releasable_db = {toml_str(paths["faa_db"])}
'''
        return make_config(tmp_config_dir, "full", toml), paths

    def test_returns_active_space_config(
        self, full_config: tuple[ActiveSpaceConfig, dict[str, Path]]
    ):
        cfg, _ = full_config
        assert isinstance(cfg, ActiveSpaceConfig)

    def test_database_fields(self, full_config: tuple[ActiveSpaceConfig, dict[str, Path]]):
        cfg, _ = full_config
        assert cfg.database == DatabaseConfig(
            name="my_db", username="user", password="pass", port="5432", host="localhost"
        )

    def test_data_paths(self, full_config: tuple[ActiveSpaceConfig, dict[str, Path]]):
        cfg, paths = full_config
        assert cfg.data.site_metadata == paths["site_metadata"]
        assert cfg.data.nvspl_archive == paths["nvspl_archive"]
        assert cfg.data.adsb == paths["adsb"]
        assert cfg.data.dem == paths["dem"]
        assert cfg.data.mennitt == paths["mennitt"]

    def test_project_paths(self, full_config: tuple[ActiveSpaceConfig, dict[str, Path]]):
        cfg, paths = full_config
        assert cfg.project.dir == paths["project_dir"]
        assert cfg.project.nmsim_binary == paths["nmsim_binary"]
        assert cfg.project.FAA_Releasable_db == paths["faa_db"]

    def test_optional_project_field_defaults_to_none(
        self, full_config: tuple[ActiveSpaceConfig, dict[str, Path]]
    ):
        cfg, _ = full_config
        assert cfg.project.FAA_type_corrections is None


class TestEmptyAndWhitespacePathsCoerceToNone:
    TOML = '''
[data]
site_metadata = ''
nvspl_archive = '   '
dem           = '\t'
'''

    def test_empty_string_is_none(self, tmp_config_dir: Path):
        cfg = make_config(tmp_config_dir, "empty", self.TOML)
        assert cfg.data.site_metadata is None

    def test_whitespace_only_is_none(self, tmp_config_dir: Path):
        cfg = make_config(tmp_config_dir, "empty", self.TOML)
        assert cfg.data.nvspl_archive is None

    def test_tab_only_is_none(self, tmp_config_dir: Path):
        cfg = make_config(tmp_config_dir, "empty", self.TOML)
        assert cfg.data.dem is None


class TestWindowsPaths:
    TOML = r'''
[data]
site_metadata = 'C:\Windows\Path\To\metadata.csv'
dem           = 'D:\Data\dem_file.tif'

[project]
dir = 'C:\Project\Dir'
'''

    def test_loads_without_error(self, tmp_config_dir: Path):
        cfg = make_config(tmp_config_dir, "win", self.TOML)
        assert isinstance(cfg, ActiveSpaceConfig)

    def test_windows_path_string_round_trips(self, tmp_config_dir: Path):
        cfg = make_config(tmp_config_dir, "win", self.TOML)
        assert str(cfg.data.site_metadata) == r"C:\Windows\Path\To\metadata.csv"
        assert str(cfg.data.dem) == r"D:\Data\dem_file.tif"
        assert str(cfg.project.dir) == r"C:\Project\Dir"

    def test_windows_path_native_semantics(self, tmp_config_dir: Path):
        cfg = make_config(tmp_config_dir, "win", self.TOML)
        win_path = PureWindowsPath(cfg.data.site_metadata)
        assert win_path.drive == "C:"
        assert win_path.name == "metadata.csv"


class TestPartialConfig:
    @pytest.fixture
    def project_only_cfg(
        self, tmp_config_dir: Path, tmp_path: Path
    ) -> tuple[ActiveSpaceConfig, Path]:
        project_dir = tmp_path / "projects"
        toml = f'''
[project]
dir = {toml_str(project_dir)}
'''
        return make_config(tmp_config_dir, "partial_proj", toml), project_dir

    @pytest.fixture
    def data_only_cfg(self, tmp_config_dir: Path, tmp_path: Path) -> ActiveSpaceConfig:
        data_dir = tmp_path / "data"
        toml = f'''
[data]
site_metadata = {toml_str(data_dir / "metadata.csv")}
nvspl_archive = {toml_str(data_dir / "nvspl")}
'''
        return make_config(tmp_config_dir, "partial_data", toml)

    def test_loads_with_only_project_section(self, project_only_cfg):
        cfg, _ = project_only_cfg
        assert isinstance(cfg, ActiveSpaceConfig)

    def test_missing_database_section_gives_empty_defaults(self, project_only_cfg):
        cfg, _ = project_only_cfg
        assert cfg.database == DatabaseConfig()

    def test_missing_data_section_gives_all_none(self, project_only_cfg):
        cfg, _ = project_only_cfg
        assert cfg.data == FileDataConfig()

    def test_missing_project_section_gives_all_none(self, data_only_cfg):
        assert data_only_cfg.project == ProjectConfig()

    def test_present_field_still_parsed(self, project_only_cfg):
        cfg, project_dir = project_only_cfg
        assert cfg.project.dir == project_dir
