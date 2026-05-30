import os

import pytest
from pathlib import Path
from unittest.mock import patch

from nps_active_space.config import (
    CONFIG_DIR_ENV_VAR,
    load_config,
    resolve_config_dir,
    resolve_config_path,
)


class TestConfigPathResolution:
    def test_env_dir_selects_named_config(self, tmp_path: Path):
        config_dir = tmp_path / "configs"
        config_dir.mkdir()
        config_file = config_dir / "DENA.toml"
        config_file.write_text("[database]\nname = 'from_env_dir'\n")
        with patch.dict(os.environ, {CONFIG_DIR_ENV_VAR: str(config_dir)}, clear=False):
            assert resolve_config_path("DENA") == config_file.resolve()
            cfg = load_config("DENA")
        assert cfg.database.name == "from_env_dir"

    def test_env_dir_does_not_fall_through_when_file_missing(self, tmp_path: Path):
        config_dir = tmp_path / "empty_configs"
        config_dir.mkdir()
        with patch.dict(os.environ, {CONFIG_DIR_ENV_VAR: str(config_dir)}, clear=False):
            with pytest.raises(FileNotFoundError, match="NPS_ACTIVE_SPACE_CONFIG_DIR"):
                resolve_config_path("missing_env")

    def test_user_config_dir_when_no_env_var(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.delenv(CONFIG_DIR_ENV_VAR, raising=False)
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("USERPROFILE", str(tmp_path))
        user_dir = tmp_path / ".nps_active_space"
        user_dir.mkdir()
        config_file = user_dir / "DENA.toml"
        config_file.write_text("[database]\nname = 'from_home'\n")
        assert resolve_config_path("DENA") == config_file.resolve()

    def test_missing_config_suggests_init_command(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.delenv(CONFIG_DIR_ENV_VAR, raising=False)
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("USERPROFILE", str(tmp_path))
        with pytest.raises(FileNotFoundError, match="nps-init-config -e DENA"):
            resolve_config_path("DENA")

    def test_resolve_config_dir_creates_user_dir(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.delenv(CONFIG_DIR_ENV_VAR, raising=False)
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("USERPROFILE", str(tmp_path))
        user_dir = tmp_path / ".nps_active_space"
        assert not user_dir.exists()
        assert resolve_config_dir() == user_dir
        assert user_dir.is_dir()
