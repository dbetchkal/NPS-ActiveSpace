from __future__ import annotations

from configparser import ConfigParser
from typing import Any, Dict, List, Optional, Union

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.utils.enums import SourceElevationUnits
import os

__all__ = [
    'initialize',
    'read',
    'read_dem_elevation_units',
    'validate',
]

_config = None

# Schema definition: maps each expected section to its known keys.
# Keys whose values represent filesystem paths are listed in _PATH_KEYS
# so that validate() can check whether those paths exist on disk.
EXPECTED_SCHEMA = {
    'database:overflights': ['name', 'username', 'password', 'port', 'host'],
    'data': ['site_metadata', 'nvspl_archive', 'adsb', 'ais', 'dem', 'mennitt'],
    'project': ['dir', 'nmsim', 'faa_releasable_db', 'faa_type_corrections'],
}

# Keys whose non-empty values should point to an existing file or directory.
_PATH_KEYS = {
    'data': {'site_metadata', 'nvspl_archive', 'adsb', 'ais', 'dem', 'mennitt'},
    'project': {'dir', 'nmsim', 'faa_releasable_db', 'faa_type_corrections'},
}

# The minimum set of keys that must have a non-empty value for the config
# to be considered usable. Other keys are allowed to be blank.
_REQUIRED_VALUES = {
    'project': {'dir'},
}


def initialize(environment: str):
    """
    Initialize a connection to a configuration file.

    Parameters
    ----------
    environment : str
        The name of the environment of the configuration file to read.

    Example
    -------
    import nps_active_space.utils.config as cfg

    cfg.initialize('production')
    """
    global _config

    config_file = os.path.join(ACTIVE_SPACE_DIR, "config", f"{environment}.config")
    _config = ConfigParser()
    _config.read(config_file)


def read(section: str, option: Optional[str] = None) -> Union[Dict, Any]:
    """
    Read in a specific section or option from a specific section from the loaded configuration file.
    The global configuration file variable must be initialized prior to reading from it.

    Parameters
    ----------
    section : str
        Section of the configuration file to read from.
    option : str, default None
        Option in the requested section of the configuration file to return.

    Returns
    ------
    If a section is requested, a dictionary of all options and values.
    If an option is requested, the option value.

    Raises
    ------
    AssertionError is config file has not been initialized

    Example
    -------
    import nps_active_space.utils.config as cfg

    cfg.initialize('../', 'production')
    user = cfg.read('database', 'username')
    """
    assert _config, "Config file initialization required before reading."

    if option:
        return _config.get(section, option)
    else:
        return dict(_config.items(section))


def read_dem_elevation_units() -> SourceElevationUnits:
    """
    Read ``dem_elevation_units`` from the loaded config's ``[data]`` section.

    Defaults to ``"feet"`` when the option is missing or blank.
    """
    assert _config, "Config file initialization required before reading."

    if _config.has_option("data", "dem_elevation_units"):
        raw = _config.get("data", "dem_elevation_units").strip()
        if raw:
            return SourceElevationUnits.parse_config_value(raw)
    return SourceElevationUnits.FEET


def validate(verbose: bool = True) -> List[str]:
    """
    Validate the currently loaded configuration.

    Checks performed:
    1. Config has been initialized (not None / empty).
    2. All expected sections are present.
    3. All expected keys exist within each section.
    4. Required values (like ``[project] dir``) are non-empty.
    5. Non-empty path values point to an existing file or directory.
    6. Unknown sections or keys are flagged (catches typos).

    Parameters
    ----------
    verbose : bool, default True
        When True, each error is printed to stdout as it is found.

    Returns
    -------
    list of str
        Human-readable error messages. An empty list means the config
        is valid.

    Example
    -------
    import nps_active_space.utils.config as cfg

    cfg.initialize('DENA_example')
    errors = cfg.validate()
    if errors:
        print(f"Found {len(errors)} config problem(s)")
    """
    errors: List[str] = []

    if _config is None or not _config.sections():
        msg = (
            "No configuration loaded. Call initialize() first, and make "
            "sure the .config file exists and is not empty."
        )
        errors.append(msg)
        if verbose:
            print(msg)
        return errors

    loaded_sections = set(_config.sections())
    expected_sections = set(EXPECTED_SCHEMA.keys())

    # --- unknown sections (typo detection) ---
    for section in sorted(loaded_sections - expected_sections):
        msg = f"Unknown section [{section}] - possible typo?"
        errors.append(msg)
        if verbose:
            print(msg)

    # --- missing sections ---
    for section in sorted(expected_sections - loaded_sections):
        msg = f"Missing required section [{section}]"
        errors.append(msg)
        if verbose:
            print(msg)

    # --- per-section key checks ---
    for section, expected_keys in EXPECTED_SCHEMA.items():
        if section not in loaded_sections:
            continue  # already reported above

        actual_keys = set(dict(_config.items(section)).keys())
        expected_set = set(expected_keys)

        # unknown keys
        for key in sorted(actual_keys - expected_set):
            msg = f"Unknown key '{key}' in [{section}] - possible typo?"
            errors.append(msg)
            if verbose:
                print(msg)

        # missing keys
        for key in sorted(expected_set - actual_keys):
            msg = f"Missing key '{key}' in [{section}]"
            errors.append(msg)
            if verbose:
                print(msg)

        # required-value checks
        required = _REQUIRED_VALUES.get(section, set())
        for key in sorted(required):
            if key not in actual_keys:
                continue  # already flagged as missing key
            value = _config.get(section, key).strip()
            if not value:
                msg = (
                    f"[{section}] {key} must not be empty - "
                    f"set it to a valid path in your .config file"
                )
                errors.append(msg)
                if verbose:
                    print(msg)

        # path existence checks (only for non-empty values)
        path_keys = _PATH_KEYS.get(section, set())
        for key in sorted(path_keys):
            if key not in actual_keys:
                continue
            value = _config.get(section, key).strip()
            if value and not os.path.exists(value):
                msg = (
                    f"[{section}] {key} = '{value}' - "
                    f"path does not exist"
                )
                errors.append(msg)
                if verbose:
                    print(msg)

    return errors
