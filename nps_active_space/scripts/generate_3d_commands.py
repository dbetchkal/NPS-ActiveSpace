"""Helpers for building per-layer generate_active_space.py command lines."""

from __future__ import annotations

import shlex
from typing import Sequence

from nps_active_space.utils.enums import AcousticModel


def build_layer_command_parts(
    environment: str,
    unit: str,
    site: str,
    year: int,
    altitude_m: int,
    ambience_arg: str,
    model: AcousticModel = AcousticModel.NMSIM,
    extra_args: Sequence[str] = (),
) -> list[str]:
    parts: list[str] = [
        "-e", environment,
        "-u", unit,
        "-s", site,
        "-y", year,
        "-l", altitude_m,
        "--model", model,
        "-a", ambience_arg,
    ]
    parts.extend(extra_args)
    return parts


def format_commands_file_line(designator: str, parts: Sequence[str]) -> str:
    return f"{designator}\t" + " ".join(shlex.quote(str(part)) for part in parts)
