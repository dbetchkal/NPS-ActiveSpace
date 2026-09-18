"""AAM runtime tuning (chunk size, subprocess timeout)."""

from __future__ import annotations

import os

# Generator ``predict()`` may pass up to ``DEFAULT_MAX_POINTS_PER_PREDICT`` (4000) points;
# AAM ``ONE TRACK`` is capped separately — see ``resolve_aam_chunk_size()``.
DEFAULT_AAM_CHUNK_SIZE = 400
AAM_RUN_TIMEOUT_S = 600

__all__ = [
    "AAM_RUN_TIMEOUT_S",
    "DEFAULT_AAM_CHUNK_SIZE",
    "resolve_aam_chunk_size",
]


def resolve_aam_chunk_size() -> int:
    """Points per AAM process (Windows native or Docker+Wine).

    Default is ``DEFAULT_AAM_CHUNK_SIZE`` (400), AAM's ``ONE TRACK`` cap /
    ``aam_translator.MAX_TRACK_POINTS``. Override with env ``AAM_CHUNK_SIZE``.
    """
    return max(1, int(os.environ.get("AAM_CHUNK_SIZE", str(DEFAULT_AAM_CHUNK_SIZE))))
