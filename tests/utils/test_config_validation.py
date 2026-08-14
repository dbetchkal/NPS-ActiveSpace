"""Tests for nps_active_space.utils.config.validate()."""

import textwrap
from pathlib import Path

import pytest

import nps_active_space.utils.config as cfg


def _write_config(tmp_path, text):
    """Write *text* into a .config file and initialize cfg from it."""
    config_file = tmp_path / "test.config"
    config_file.write_text(textwrap.dedent(text))
    # Point the module-level _config at the file we just wrote.
    from configparser import ConfigParser
    parser = ConfigParser()
    parser.read(str(config_file))
    cfg._config = parser


# -- helpers for building valid configs -----------------------------------

VALID_CONFIG = """\
[database:overflights]
name = mydb
username = user
password = secret
port = 5432
host = localhost

[data]
site_metadata =
nvspl_archive =
adsb =
ais =
dem =
dem_elevation_units = feet
mennitt =

[project]
dir = {project_dir}
nmsim =
faa_releasable_db =
faa_type_corrections =
"""


class TestValidateHappyPath:
    """A well-formed config should produce zero errors."""

    def test_valid_config_returns_no_errors(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        _write_config(tmp_path, VALID_CONFIG.format(project_dir=project_dir))
        errors = cfg.validate(verbose=False)
        assert errors == []

    def test_blank_optional_paths_are_fine(self, tmp_path):
        """Keys like nmsim or adsb can be blank without triggering errors."""
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        _write_config(tmp_path, VALID_CONFIG.format(project_dir=project_dir))
        errors = cfg.validate(verbose=False)
        assert errors == []


class TestValidateUninitialized:
    """validate() called before initialize() should report that."""

    def test_none_config(self):
        cfg._config = None
        errors = cfg.validate(verbose=False)
        assert len(errors) == 1
        assert "No configuration loaded" in errors[0]

    def test_empty_config(self, tmp_path):
        """An empty file produces a ConfigParser with no sections."""
        _write_config(tmp_path, "")
        errors = cfg.validate(verbose=False)
        assert len(errors) == 1
        assert "No configuration loaded" in errors[0]


class TestValidateMissingSections:
    def test_missing_project_section(self, tmp_path):
        _write_config(tmp_path, """\
        [database:overflights]
        name =
        username =
        password =
        port =
        host =

        [data]
        site_metadata =
        nvspl_archive =
        adsb =
        ais =
        dem =
        dem_elevation_units = feet
        mennitt =
        """)
        errors = cfg.validate(verbose=False)
        assert any("Missing required section [project]" in e for e in errors)
        assert len(errors) == 1

    def test_missing_data_section(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        _write_config(tmp_path, f"""\
        [database:overflights]
        name =
        username =
        password =
        port =
        host =

        [project]
        dir = {project_dir}
        nmsim =
        FAA_Releasable_db =
        FAA_type_corrections =
        """)
        errors = cfg.validate(verbose=False)
        assert any("Missing required section [data]" in e for e in errors)
        assert len(errors) == 1


class TestValidateMissingKeys:
    def test_missing_key_in_project(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        _write_config(tmp_path, f"""\
        [database:overflights]
        name =
        username =
        password =
        port =
        host =

        [data]
        site_metadata =
        nvspl_archive =
        adsb =
        ais =
        dem =
        dem_elevation_units = feet
        mennitt =

        [project]
        dir = {project_dir}
        nmsim =
        """)
        errors = cfg.validate(verbose=False)
        missing = [e for e in errors if "Missing key" in e and "[project]" in e]
        key_names = " ".join(missing)
        assert "faa_releasable_db" in key_names.lower()
        assert "faa_type_corrections" in key_names.lower()
        assert len(errors) == 2


class TestValidateUnknownSectionsAndKeys:
    def test_unknown_section(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        _write_config(tmp_path, VALID_CONFIG.format(project_dir=project_dir) + """\
[bogus]
foo = bar
""")
        errors = cfg.validate(verbose=False)
        assert any("Unknown section [bogus]" in e for e in errors)
        assert len(errors) == 1

    def test_unknown_key(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        text = VALID_CONFIG.format(project_dir=project_dir).replace(
            "[project]",
            "[project]\ntypo_key = something"
        )
        _write_config(tmp_path, text)
        errors = cfg.validate(verbose=False)
        assert any("Unknown key 'typo_key' in [project]" in e for e in errors)
        assert len(errors) == 1


class TestValidateRequiredValues:
    def test_empty_project_dir(self, tmp_path):
        _write_config(tmp_path, VALID_CONFIG.format(project_dir=""))
        errors = cfg.validate(verbose=False)
        assert any(
            "[project] dir must not be empty" in e
            for e in errors
        )
        assert len(errors) == 1


class TestValidatePathExistence:
    def test_nonexistent_path(self, tmp_path):
        _write_config(tmp_path, VALID_CONFIG.format(
            project_dir="/does/not/exist/at/all"
        ))
        errors = cfg.validate(verbose=False)
        assert any(
            "path does not exist" in e and "[project] dir" in e
            for e in errors
        )
        assert len(errors) == 1

    def test_nonexistent_data_path(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        text = VALID_CONFIG.format(project_dir=project_dir).replace(
            "nvspl_archive =",
            "nvspl_archive = /no/such/archive"
        )
        _write_config(tmp_path, text)
        errors = cfg.validate(verbose=False)
        assert any(
            "path does not exist" in e and "nvspl_archive" in e
            for e in errors
        )
        assert len(errors) == 1

    def test_existing_path_no_error(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        archive = tmp_path / "archive"
        archive.mkdir()
        text = VALID_CONFIG.format(project_dir=project_dir).replace(
            "nvspl_archive =",
            f"nvspl_archive = {archive}"
        )
        _write_config(tmp_path, text)
        errors = cfg.validate(verbose=False)
        assert errors == []


class TestValidateVerboseOutput:
    def test_verbose_prints(self, tmp_path, capsys):
        cfg._config = None
        cfg.validate(verbose=True)
        captured = capsys.readouterr()
        assert "No configuration loaded" in captured.out

    def test_quiet_does_not_print(self, tmp_path, capsys):
        cfg._config = None
        cfg.validate(verbose=False)
        captured = capsys.readouterr()
        assert captured.out == ""


class TestValidateDemElevationUnits:
    def test_invalid_dem_elevation_units(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        text = VALID_CONFIG.format(project_dir=project_dir).replace(
            "dem_elevation_units = feet",
            "dem_elevation_units = yards",
        )
        _write_config(tmp_path, text)
        errors = cfg.validate(verbose=False)
        assert any(
            "[data] dem_elevation_units must be 'feet' or 'meters', got 'yards'"
            in e
            for e in errors
        )
        assert len(errors) == 1

    def test_blank_dem_elevation_units_is_fine(self, tmp_path):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        text = VALID_CONFIG.format(project_dir=project_dir).replace(
            "dem_elevation_units = feet",
            "dem_elevation_units =",
        )
        _write_config(tmp_path, text)
        errors = cfg.validate(verbose=False)
        assert errors == []


class TestInitializeValidation:
    def _write_env_config(self, tmp_path, name: str, body: str) -> None:
        config_dir = tmp_path / "config"
        config_dir.mkdir(parents=True, exist_ok=True)
        (config_dir / f"{name}.config").write_text(textwrap.dedent(body))

    def test_initialize_raises_on_invalid_config(self, tmp_path, monkeypatch):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        self._write_env_config(
            tmp_path,
            "invalid_units",
            VALID_CONFIG.format(project_dir=project_dir).replace(
                "dem_elevation_units = feet",
                "dem_elevation_units = bad",
            ),
        )
        monkeypatch.setattr(cfg, "ACTIVE_SPACE_DIR", str(tmp_path))
        with pytest.raises(ValueError, match="Configuration validation failed"):
            cfg.initialize("invalid_units", validate_config=True)

    def test_initialize_skips_validation_when_disabled(self, tmp_path, monkeypatch):
        project_dir = tmp_path / "projects"
        project_dir.mkdir()
        self._write_env_config(
            tmp_path,
            "invalid_units",
            VALID_CONFIG.format(project_dir=project_dir).replace(
                "dem_elevation_units = feet",
                "dem_elevation_units = bad",
            ),
        )
        monkeypatch.setattr(cfg, "ACTIVE_SPACE_DIR", str(tmp_path))
        cfg.initialize("invalid_units", validate_config=False)
        errors = cfg.validate(verbose=False)
        assert len(errors) >= 1


@pytest.fixture
def repo_root():
    return Path(__file__).resolve().parents[2]


class TestExampleConfigs:
    @pytest.mark.parametrize("environment", ["DENA_example", "GLBA_example"])
    def test_shipped_example_configs_validate(
        self, environment, repo_root, monkeypatch
    ):
        monkeypatch.chdir(repo_root)
        cfg.initialize(environment, validate_config=False)
        errors = cfg.validate(verbose=False)
        assert errors == []

    def test_full_example_data_fixture_validates(self, repo_root, monkeypatch):
        from configparser import ConfigParser

        monkeypatch.chdir(repo_root)
        fixture_path = repo_root / "tests" / "fixtures" / "full_example_data.config"
        parser = ConfigParser()
        parser.read(str(fixture_path))
        cfg._config = parser
        errors = cfg.validate(verbose=False)
        assert errors == []
