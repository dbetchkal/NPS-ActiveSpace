from configparser import ConfigParser
from typing import Any, Dict, Optional, Union

__all__ = [
    'initialize',
    'read'
]

_config = None


def initialize(config_dir: str, environment: str):
    """Initialize a connection to a configuration file.

    Parameters
    ----------
        config_dir : str
            The name of the directory where the config file is stored.
        environment: str
            The name of the environment of the configuration file to read.

    Example
    -------
        import _DENA.resource.config as cfg

        cfg.initialize('../', 'production')
    """
    global _config

    config_file = f"{config_dir}/{environment}.config"
    _config = ConfigParser()
    _config.read(config_file)


def read(section: str, option: Optional[str] = None) -> Union[Dict, Any]:
    """Read in a specific section or option from a specific section from the loaded configuration file.
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

    Example
    -------
        import _DENA.resource.config as cfg

        cfg.initialize('../', 'production')
        user = cfg.read('database', 'username')
    """
    assert _config, "Config file initialization required before reading."

    if option:
        return _config.get(section, option)
    else:
        return dict(_config.items(section))
