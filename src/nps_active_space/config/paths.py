import os
from pathlib import Path

CONFIG_DIR_ENV_VAR = "NPS_ACTIVE_SPACE_CONFIG_DIR"
USER_CONFIG_SUBDIR = ".nps_active_space"
BUNDLED_TEMPLATE_NAME = "template.toml"

_CONFIG_PACKAGE_DIR = Path(__file__).resolve().parent


def user_config_dir() -> Path:
    return Path.home() / USER_CONFIG_SUBDIR


def resolve_config_dir(*, create_user_dir: bool = True) -> Path:
    if env_config_dir := os.environ.get(CONFIG_DIR_ENV_VAR):
        return Path(env_config_dir).expanduser().resolve()
    user_dir = user_config_dir()
    if create_user_dir:
        user_dir.mkdir(parents=True, exist_ok=True)
    return user_dir


def bundled_config_template_path() -> Path:
    return _CONFIG_PACKAGE_DIR / BUNDLED_TEMPLATE_NAME


def config_filename(environment: str) -> str:
    return f"{environment}.toml"


def resolve_config_path(environment: str) -> Path:
    """Return the TOML path for *environment*.

    Search order:

    1. If :envvar:`NPS_ACTIVE_SPACE_CONFIG_DIR` is set, only
       ``{that_dir}/{environment}.toml``.
    2. Otherwise ``~/.nps_active_space/{environment}.toml``.

    CLI ``-e`` / ``--environment`` selects the config name (e.g. ``-e DENA`` ->
    ``DENA.toml``).
    """
    config_dir = resolve_config_dir(create_user_dir=False)
    path = config_dir / config_filename(environment)
    if path.is_file():
        return path.resolve()

    if os.environ.get(CONFIG_DIR_ENV_VAR):
        label = f"{CONFIG_DIR_ENV_VAR} ({config_dir})"
    else:
        label = str(config_dir)
    raise FileNotFoundError(
        f"No config file for environment {environment!r} at {path}\n"
        f"  ({label})\n"
        f"Run: nps-init-config -e {environment}"
    )
