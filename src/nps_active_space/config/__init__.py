"""Environment TOML schema, discovery, and loading."""

from nps_active_space.config.io import (
    config_setup_hint,
    init_config,
    load_config,
    show_config,
)
from nps_active_space.config.models import (
    ActiveSpaceConfig,
    DatabaseConfig,
    FileDataConfig,
    ProjectConfig,
)
from nps_active_space.config.paths import (
    BUNDLED_TEMPLATE_NAME,
    CONFIG_DIR_ENV_VAR,
    USER_CONFIG_SUBDIR,
    bundled_config_template_path,
    config_filename,
    resolve_config_dir,
    resolve_config_path,
    user_config_dir,
)

__all__ = [
    "CONFIG_DIR_ENV_VAR",
    "USER_CONFIG_SUBDIR",
    "BUNDLED_TEMPLATE_NAME",
    "user_config_dir",
    "resolve_config_dir",
    "resolve_config_path",
    "bundled_config_template_path",
    "config_filename",
    "init_config",
    "config_setup_hint",
    "load_config",
    "ActiveSpaceConfig",
    "DatabaseConfig",
    "FileDataConfig",
    "ProjectConfig",
    "show_config",
]
