import tomllib
from pathlib import Path

import pytest

from nps_active_space.config import bundled_config_template_path, init_config
from nps_active_space.config.models import ActiveSpaceConfig


class TestBundledTemplate:
    def test_template_exists(self):
        path = bundled_config_template_path()
        assert path.is_file()

    def test_template_round_trips_to_defaults(self):
        with open(bundled_config_template_path(), "rb") as f:
            data = tomllib.load(f)
        assert ActiveSpaceConfig(**data) == ActiveSpaceConfig()


class TestInitConfig:
    def test_writes_environment_toml(self, tmp_path: Path):
        dest = init_config("DENA", config_dir=tmp_path)
        assert dest == (tmp_path / "DENA.toml").resolve()
        assert dest.read_text() == bundled_config_template_path().read_text()

    def test_refuses_overwrite_without_force(self, tmp_path: Path):
        init_config("DENA", config_dir=tmp_path)
        with pytest.raises(FileExistsError, match="already exists"):
            init_config("DENA", config_dir=tmp_path)

    def test_force_overwrites(self, tmp_path: Path):
        init_config("DENA", config_dir=tmp_path)
        (tmp_path / "DENA.toml").write_text("[database]\nname = 'stale'\n")
        init_config("DENA", config_dir=tmp_path, force=True)
        with open(tmp_path / "DENA.toml", "rb") as f:
            assert tomllib.load(f)["database"]["name"] == ""
