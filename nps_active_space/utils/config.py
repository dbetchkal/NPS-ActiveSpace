from __future__ import annotations

from configparser import ConfigParser
from typing import Any, Dict, Optional, Union

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.utils.enums import SourceElevationUnits
import os

__all__ = [
    'initialize',
    'read',
    'read_dem_elevation_units',
]

_config = None


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
