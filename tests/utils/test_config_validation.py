"""Tests for nps_active_space.utils.config.validate()."""

import textwrap

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
        mennitt =
        """)
        errors = cfg.validate(verbose=False)
        assert any("Missing required section [project]" in e for e in errors)

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


class TestValidateRequiredValues:
    def test_empty_project_dir(self, tmp_path):
        _write_config(tmp_path, VALID_CONFIG.format(project_dir=""))
        errors = cfg.validate(verbose=False)
        assert any(
            "[project] dir must not be empty" in e
            for e in errors
        )


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
