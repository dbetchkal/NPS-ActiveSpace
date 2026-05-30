import os
from pathlib import Path
from unittest.mock import patch

from nps_active_space.config import CONFIG_DIR_ENV_VAR, ActiveSpaceConfig, load_config


def make_config(tmp_config_dir: Path, name: str, toml: str) -> ActiveSpaceConfig:
    """Write *toml* to *name*.toml under *tmp_config_dir*, then load and return the config."""
    config_file = tmp_config_dir / f"{name}.toml"
    config_file.write_text(toml)
    with patch.dict(os.environ, {CONFIG_DIR_ENV_VAR: str(tmp_config_dir)}, clear=False):
        return load_config(name)


def toml_str(path: Path) -> str:
    """Escape a Path as a TOML single-quoted literal string."""
    return "'" + str(path).replace("'", "''") + "'"
