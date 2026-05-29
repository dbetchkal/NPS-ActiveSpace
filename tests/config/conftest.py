from pathlib import Path

import pytest


@pytest.fixture
def tmp_config_dir(tmp_path: Path) -> Path:
    """Temporary directory that mimics a config folder (e.g. NPS_ACTIVE_SPACE_CONFIG_DIR)."""
    config_dir = tmp_path / "configs"
    config_dir.mkdir()
    return config_dir
