from configparser import ConfigParser
import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Union
from nps_active_space import ACTIVE_SPACE_DIR
import os

__all__ = [
    'initialize',
    'read',
    'show',
    'NPSConfig',
    'active_config'
]

@dataclass
class DatabaseConfig:
    name: str = ""
    username: str = ""
    password: str = ""
    port: str = ""
    host: str = ""

@dataclass
class FileDataConfig:
    site_metadata: Optional[Path] = None
    nvspl_archive: Optional[Path] = None
    adsb: Optional[Path] = None
    dem: Optional[Path] = None
    mennitt: Optional[Path] = None

@dataclass
class ProjectConfig:
    dir: Optional[Path] = None
    nmsim: Optional[Path] = None
    FAA_Releasable_db: Optional[Path] = None
    FAA_type_corrections: Optional[Path] = None

@dataclass
class ActiveSpaceConfig:
    database: DatabaseConfig = field(default_factory=DatabaseConfig)
    data: FileDataConfig = field(default_factory=FileDataConfig)
    project: ProjectConfig = field(default_factory=ProjectConfig)

_config: Optional[ConfigParser] = None
active_config: Optional[ActiveSpaceConfig] = None


def _parse_path(val: str) -> Optional[Path]:
    val = val.strip() if val else ""
    return Path(val) if val else None


def initialize(environment: str, verbose: bool = False):
    """
    Initialize a connection to a configuration file.

    Parameters
    ----------
    environment : str
        The name of the environment of the configuration file to read.
    verbose : bool, default False
        If True, print the configuration values after loading.

    Example
    -------
    import nps_active_space.utils.config as cfg

    cfg.initialize('production', verbose=True)
    """
    global _config, active_config

    config_file = os.path.join(ACTIVE_SPACE_DIR, "config", f"{environment}.config")
    _config = ConfigParser()
    _config.read(config_file)

    db_sec = "database:overflights"
    active_config = ActiveSpaceConfig(
        database=DatabaseConfig(
            name=_config.get(db_sec, "name", fallback=""),
            username=_config.get(db_sec, "username", fallback=""),
            password=_config.get(db_sec, "password", fallback=""),
            port=_config.get(db_sec, "port", fallback=""),
            host=_config.get(db_sec, "host", fallback="")
        ),
        data=FileDataConfig(
            site_metadata=_parse_path(_config.get("data", "site_metadata", fallback="")),
            nvspl_archive=_parse_path(_config.get("data", "nvspl_archive", fallback="")),
            adsb=_parse_path(_config.get("data", "adsb", fallback="")),
            dem=_parse_path(_config.get("data", "dem", fallback="")),
            mennitt=_parse_path(_config.get("data", "mennitt", fallback=""))
        ),
        project=ProjectConfig(
            dir=_parse_path(_config.get("project", "dir", fallback="")),
            nmsim=_parse_path(_config.get("project", "nmsim", fallback="")),
            FAA_Releasable_db=_parse_path(_config.get("project", "FAA_Releasable_db", fallback="")),
            FAA_type_corrections=_parse_path(_config.get("project", "FAA_type_corrections", fallback=""))
        )
    )

    if verbose:
        show()


def read(section: str, option: Optional[str] = None, as_path: bool = False) -> Union[Dict[str, Any], Any]:
    """
    Legacy method to read config values. 
    It is recommended to use `active_config` directly.

    Read in a specific section or option from a specific section from the loaded configuration file.
    The global configuration file variable must be initialized prior to reading from it.

    Parameters
    ----------
    section : str
        Section of the configuration file to read from.
    option : str, default None
        Option in the requested section of the configuration file to return.
    as_path : bool, default False
        If True, return the value(s) as Path object(s).

    Returns
    ------
    If a section is requested, a dictionary of all options and values.
    If an option is requested, the option value.
    Values will be converted to Path objects if as_path is True and value is not empty.

    Raises
    ------
    AssertionError if config file has not been initialized

    Example
    -------
    import nps_active_space.utils.config as cfg

    cfg.initialize('production')
    user = cfg.read('database', 'username')
    data_paths = cfg.read('data', as_path=True)
    """
    assert _config, "Config file initialization required before reading."

    if option:
        val = _config.get(section, option)
        return Path(val) if as_path and val else val
    else:
        items = dict(_config.items(section))
        if as_path:
            return {k: Path(v) if v else v for k, v in items.items()}
        return items


def show():
    """
    Print a formatted view of the current configuration.
    
    Raises
    ------
    AssertionError if config file has not been initialized
    """
    assert active_config, "Config file initialization required before reading."
    
    print("\n" + "="*50)
    print("CURRENT CONFIGURATION")
    print("="*50)
    
    for sec_name, sec_obj in [("DATABASE:OVERFLIGHTS", active_config.database), 
                              ("DATA", active_config.data), 
                              ("PROJECT", active_config.project)]:
        print(f"\n[{sec_name}]")
        for f in dataclasses.fields(sec_obj):
            val = getattr(sec_obj, f.name)
            print(f"  {f.name:<20} : {val}")
    print("\n" + "="*50 + "\n")
