import os
import shutil
import tomllib
from pathlib import Path

from nps_active_space.config.models import ActiveSpaceConfig
from nps_active_space.config.paths import (
    CONFIG_DIR_ENV_VAR,
    bundled_config_template_path,
    config_filename,
    resolve_config_dir,
    resolve_config_path,
    user_config_dir,
)


def init_config(
    environment: str,
    *,
    config_dir: Path | None = None,
    force: bool = False,
) -> Path:
    """Copy the bundled template to ``{config_dir}/{environment}.toml``.

    Returns the path to the new file. Raises :class:`FileExistsError` if the
    destination exists unless *force* is true.
    """
    template = bundled_config_template_path()
    if not template.is_file():
        raise FileNotFoundError(f"Bundled config template not found: {template}")

    dest_dir = config_dir if config_dir is not None else resolve_config_dir()
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / config_filename(environment)

    if dest.exists() and not force:
        raise FileExistsError(
            f"Config already exists: {dest}\n"
            "Use --force to overwrite, or choose a different -e name."
        )

    shutil.copy2(template, dest)
    return dest.resolve()


def config_setup_hint(config_dir: Path) -> str | None:
    """Optional one-line guidance after ``init_config``; ``None`` if redundant."""
    resolved_dir = config_dir.expanduser().resolve()
    if env_dir := os.environ.get(CONFIG_DIR_ENV_VAR):
        env_path = Path(env_dir).expanduser().resolve()
        if env_path == resolved_dir:
            return None
        return (
            f"Set {CONFIG_DIR_ENV_VAR} to {resolved_dir} so pipeline scripts "
            "load configs from that folder."
        )

    if resolved_dir == user_config_dir().resolve():
        return (
            f"Using default config folder {resolved_dir}. "
            f"Optional: set {CONFIG_DIR_ENV_VAR} to a shared folder "
            "(e.g. a network drive) if your team keeps configs in one place."
        )

    return (
        f"Set {CONFIG_DIR_ENV_VAR} to {resolved_dir} so pipeline scripts "
        "load configs from that folder."
    )


def load_config(environment: str) -> ActiveSpaceConfig:
    config_file = resolve_config_path(environment)
    with open(config_file, "rb") as f:
        config_fields = tomllib.load(f)

    return ActiveSpaceConfig(**config_fields)


def show_config(config: ActiveSpaceConfig) -> None:
    print("\n" + "=" * 50)
    print("CURRENT CONFIGURATION")
    print("=" * 50)

    for sec_name, sec_obj in [
        ("DATABASE", config.database),
        ("DATA", config.data),
        ("PROJECT", config.project),
    ]:
        print(f"\n[{sec_name}]")
        for name, val in sec_obj.model_dump().items():
            print(f"  {name:<20} : {val}")
    print("\n" + "=" * 50 + "\n")
