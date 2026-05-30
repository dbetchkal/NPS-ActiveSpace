from pathlib import Path

from pydantic import BaseModel, Field, field_validator


def _empty_str_to_none(v):
    if isinstance(v, str) and not v.strip():
        return None
    return v


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
        return _empty_str_to_none(v)


class ProjectConfig(BaseModel):
    dir: Path | None = None
    nmsim_binary: Path | None = None
    FAA_Releasable_db: Path | None = None
    FAA_type_corrections: Path | None = None

    @field_validator("*", mode="before")
    @classmethod
    def empty_str_to_none(cls, v):
        return _empty_str_to_none(v)


class ActiveSpaceConfig(BaseModel):
    database: DatabaseConfig = Field(default_factory=DatabaseConfig)
    data: FileDataConfig = Field(default_factory=FileDataConfig)
    project: ProjectConfig = Field(default_factory=ProjectConfig)
