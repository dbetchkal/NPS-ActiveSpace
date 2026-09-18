from __future__ import annotations

from enum import StrEnum


FEET_TO_METERS = 0.3048


class TrackSource(StrEnum):
    """External flight-track feed used for ground truthing and metrics."""

    GPS = 'GPS'
    ADSB = 'ADSB'
    AIS = 'AIS'


class AcousticModel(StrEnum):
    """Propagation backend for active-space generation, fit, and viz."""

    NMSIM = "nmsim"
    AAM = "aam"

    @classmethod
    def parse(cls, value: str | AcousticModel) -> AcousticModel:
        """Parse a CLI, env, or CSV model name."""
        if isinstance(value, cls):
            return value
        try:
            return cls(str(value).strip().lower())
        except ValueError as exc:
            names = ", ".join(repr(member.value) for member in cls)
            raise ValueError(
                f"Unknown acoustic model {value!r}; expected {names}"
            ) from exc


class SourceElevationUnits(StrEnum):
    """Vertical units of a source DEM raster before conversion to meters."""

    FEET = "feet"
    METERS = "meters"

    def to_meters_scale(self) -> float:
        return FEET_TO_METERS if self is SourceElevationUnits.FEET else 1.0

    @classmethod
    def parse_config_value(cls, units: str) -> SourceElevationUnits:
        """Parse a ``dem_elevation_units`` config string."""
        try:
            return cls(units.strip().lower())
        except ValueError as exc:
            raise ValueError(
                f"dem_elevation_units must be 'feet' or 'meters', got {units!r}"
            ) from exc
