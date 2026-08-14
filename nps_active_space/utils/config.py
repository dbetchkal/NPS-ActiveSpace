from __future__ import annotations

from configparser import ConfigParser
from typing import Any, Dict, List, Optional, Set, Union

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
    'data': ['site_metadata', 'nvspl_archive', 'adsb', 'ais', 'dem', 'dem_elevation_units', 'mennitt'],
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


def initialize(environment: str, *, validate_config: bool = True):
    """
    Initialize a connection to a configuration file.

    Parameters
    ----------
    environment : str
        The name of the environment of the configuration file to read.
    validate_config : bool, default True
        When True, run :func:`validate` after loading and raise if errors are found.

    Raises
    ------
    ValueError
        If ``validate_config`` is True and validation reports one or more errors.

    Example
    -------
    import nps_active_space.utils.config as cfg

    cfg.initialize('production')
    """
    global _config

    config_file = os.path.join(ACTIVE_SPACE_DIR, "config", f"{environment}.config")
    _config = ConfigParser()
    _config.read(config_file)

    if validate_config:
        errors = validate(verbose=True)
        if errors:
            raise ValueError(
                f"Configuration validation failed with {len(errors)} error(s): "
                f"{errors[0]}"
            )


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


def _check_loaded(config: Optional[ConfigParser]) -> List[str]:
    """Return a single error if the config has not been loaded or is empty."""
    if config is None or not config.sections():
        return [
            "No configuration loaded. Call initialize() first, and make "
            "sure the .config file exists and is not empty."
        ]
    return []


def _check_unknown_sections(
    loaded: Set[str], expected: Set[str]
) -> List[str]:
    """Flag sections present in the file but not in the schema."""
    errors: List[str] = []
    for section in sorted(loaded - expected):
        errors.append(f"Unknown section [{section}] - possible typo?")
    return errors


def _check_missing_sections(
    loaded: Set[str], expected: Set[str]
) -> List[str]:
    """Flag required sections that are absent from the file."""
    errors: List[str] = []
    for section in sorted(expected - loaded):
        errors.append(f"Missing required section [{section}]")
    return errors


def _check_section_keys(
    config: ConfigParser,
    section: str,
    expected_keys: List[str],
    path_keys: Set[str],
    required_values: Set[str],
) -> List[str]:
    """Validate keys, required values, and path existence for one section."""
    errors: List[str] = []
    actual_keys = set(dict(config.items(section)).keys())
    expected_set = set(expected_keys)

    for key in sorted(actual_keys - expected_set):
        errors.append(f"Unknown key '{key}' in [{section}] - possible typo?")

    for key in sorted(expected_set - actual_keys):
        errors.append(f"Missing key '{key}' in [{section}]")

    for key in sorted(required_values):
        if key not in actual_keys:
            continue
        value = config.get(section, key).strip()
        if not value:
            errors.append(
                f"[{section}] {key} must not be empty - "
                f"set it to a valid path in your .config file"
            )

    for key in sorted(path_keys):
        if key not in actual_keys:
            continue
        value = config.get(section, key).strip()
        if value and not os.path.exists(value):
            errors.append(
                f"[{section}] {key} = '{value}' - "
                f"path does not exist"
            )

    return errors


def _check_dem_elevation_units(config: ConfigParser) -> List[str]:
    """Validate non-empty ``dem_elevation_units`` in ``[data]``."""
    errors: List[str] = []
    if not config.has_section("data"):
        return errors
    if not config.has_option("data", "dem_elevation_units"):
        return errors
    raw = config.get("data", "dem_elevation_units").strip()
    if not raw:
        return errors
    try:
        SourceElevationUnits.parse_config_value(raw)
    except ValueError:
        errors.append(
            f"[data] dem_elevation_units must be 'feet' or 'meters', got '{raw}'"
        )
    return errors


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
    7. Non-empty ``dem_elevation_units`` values are ``feet`` or ``meters``.

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
    load_errors = _check_loaded(_config)
    if load_errors:
        if verbose:
            for msg in load_errors:
                print(msg)
        return load_errors

    loaded_sections = set(_config.sections())
    expected_sections = set(EXPECTED_SCHEMA.keys())

    errors: List[str] = []
    errors.extend(_check_unknown_sections(loaded_sections, expected_sections))
    errors.extend(_check_missing_sections(loaded_sections, expected_sections))

    for section, expected_keys in EXPECTED_SCHEMA.items():
        if section not in loaded_sections:
            continue
        errors.extend(_check_section_keys(
            _config,
            section,
            expected_keys,
            _PATH_KEYS.get(section, set()),
            _REQUIRED_VALUES.get(section, set()),
        ))

    errors.extend(_check_dem_elevation_units(_config))

    if verbose:
        for msg in errors:
            print(msg)

    return errors
