import pytest
from pathlib import Path
from unittest.mock import patch

from nps_active_space.utils import config
from nps_active_space.utils.config import ActiveSpaceConfig

@pytest.fixture
def temp_config_dir(tmp_path):
    """Create a temporary config directory resembling the actual project structure."""
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    return config_dir

def test_config_parsing_unix_paths(temp_config_dir):
    config_content = """
[database]
name = "my_db"
username = "user"
password = "pass"
port = "5432"
host = "localhost"

[data]
site_metadata = "/unix/path/to/metadata.csv"
nvspl_archive = "/unix/path/to/nvspl"

[project]
dir = "/unix/path/to/project"
"""
    (temp_config_dir / "test_unix.toml").write_text(config_content)
    
    with patch("nps_active_space.utils.config.ACTIVE_SPACE_DIR", str(temp_config_dir.parent)):
        cfg = config.load_config("test_unix")
    
    assert cfg is not None
    assert isinstance(cfg, ActiveSpaceConfig)
    assert cfg.database.name == "my_db"
    assert cfg.database.port == "5432"
    assert cfg.data.site_metadata == Path("/unix/path/to/metadata.csv")
    assert cfg.data.nvspl_archive == Path("/unix/path/to/nvspl")
    assert cfg.project.dir == Path("/unix/path/to/project")

def test_config_parsing_windows_paths(temp_config_dir):
    config_content = """
[database]
name = "windows_db"

[data]
site_metadata = "C:\\\\Windows\\\\Path\\\\To\\\\metadata.csv"
dem = "D:\\\\Data\\\\dem_file.tif"

[project]
dir = "C:\\\\Project\\\\Dir"
"""
    (temp_config_dir / "test_win.toml").write_text(config_content)
    
    with patch("nps_active_space.utils.config.ACTIVE_SPACE_DIR", str(temp_config_dir.parent)):
        cfg = config.load_config("test_win")
    
    assert cfg is not None
    assert isinstance(cfg, ActiveSpaceConfig)
    assert cfg.database.name == "windows_db"
    assert cfg.data.site_metadata == Path("C:\\Windows\\Path\\To\\metadata.csv")
    assert cfg.data.dem == Path("D:\\Data\\dem_file.tif")
    assert cfg.project.dir == Path("C:\\Project\\Dir")

def test_empty_paths_handled_correctly(temp_config_dir):
    config_content = """
[data]
site_metadata = ""
nvspl_archive = "  "
"""
    (temp_config_dir / "test_empty.toml").write_text(config_content)
    
    with patch("nps_active_space.utils.config.ACTIVE_SPACE_DIR", str(temp_config_dir.parent)):
        cfg = config.load_config("test_empty")
        
    assert isinstance(cfg, ActiveSpaceConfig)
    assert cfg.data.site_metadata is None
    assert cfg.data.nvspl_archive is None
