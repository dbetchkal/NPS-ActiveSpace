import tomllib
from pathlib import Path
from pydantic import BaseModel, Field, field_validator
from nps_active_space import ACTIVE_SPACE_DIR

__all__ = [
    'load_config',
    'ActiveSpaceConfig',
    'DatabaseConfig',
    'FileDataConfig',
    'ProjectConfig',
    'show_config'
]

class DatabaseConfig(BaseModel):
    name: str = ""
    username: str = ""
    password: str = ""
    port: str = ""
    host: str = ""

class FileDataConfig(BaseModel):
    site_metadata: Path | None = None
    nvspl_archive: Path | None = None
    adsb: Path | None = None
    dem: Path | None = None
    mennitt: Path | None = None

    @field_validator("*", mode="before")
    @classmethod
    def empty_str_to_none(cls, v):
        if isinstance(v, str) and not v.strip():
            return None
        return v

class ProjectConfig(BaseModel):
    dir: Path | None = None
    nmsim_binary: Path | None = None
    FAA_Releasable_db: Path | None = None
    FAA_type_corrections: Path | None = None

    @field_validator("*", mode="before")
    @classmethod
    def empty_str_to_none(cls, v):
        if isinstance(v, str) and not v.strip():
            return None
        return v

class ActiveSpaceConfig(BaseModel):
    database: DatabaseConfig = Field(default_factory=DatabaseConfig)
    data: FileDataConfig = Field(default_factory=FileDataConfig)
    project: ProjectConfig = Field(default_factory=ProjectConfig)

def load_config(environment: str) -> ActiveSpaceConfig:
    """
    Load a configuration file for a given environment and return the parsed ActiveSpaceConfig.

    Parameters
    ----------
    environment : str
        The name of the environment (e.g., 'production').

    Returns
    -------
    ActiveSpaceConfig
        The parsed and typed configuration model.
    """
    config_file = Path(ACTIVE_SPACE_DIR) / "config" / f"{environment}.toml"
    with open(config_file, "rb") as f:
        data = tomllib.load(f)
    
    # Map the legacy/alternative "database:overflights" section to "database" if present
    if "database:overflights" in data:
        data["database"] = data.pop("database:overflights")
        
    return ActiveSpaceConfig(**data)

def show_config(config: ActiveSpaceConfig) -> None:
    """Print a formatted view of the configuration."""
    print("\n" + "="*50)
    print("CURRENT CONFIGURATION")
    print("="*50)
    
    for sec_name, sec_obj in [("DATABASE", config.database), 
                              ("DATA", config.data), 
                              ("PROJECT", config.project)]:
        print(f"\n[{sec_name}]")
        for name, val in sec_obj.model_dump().items():
            print(f"  {name:<20} : {val}")
    print("\n" + "="*50 + "\n")
