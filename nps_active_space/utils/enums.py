from enum import StrEnum


class TrackSource(StrEnum):
    """External flight-track feed used for ground truthing and metrics."""

    GPS = 'GPS'
    ADSB = 'ADSB'
    AIS = 'AIS'
