"""Legacy NMSim output path resolvers (read-side backward compatibility only)."""

from __future__ import annotations

import glob
import os
import re

from nps_active_space.utils.paths import (
    NMSIM_OUTPUT_SUBDIR,
    NMSIM_PREDICTIONS_SUBDIR,
    NMSIM_SCRATCH_SUBDIR,
)

LEGACY_PREDICTIONS_SUBDIR = "Output_Data/TIG_TIS"
LEGACY_ACTIVESPACES_SUBDIR = "Output_Data/ACTIVESPACES"

NMSIM_ACTIVESPACES_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/ACTIVESPACES"


def _join(site_dir: str, *parts: str) -> str:
    return os.path.join(site_dir, *parts)


def _deployment_id(unit: str, site: str, year) -> str:
    return f"{unit}{site}{year}"


def is_standard_altitude_layer_dir(layer_dir: str, usy: str) -> bool:
    """True for ``DENATRLA2025_1000m``; false for experiment dirs like ``DENATRLA2025_1000m_aam``."""
    name = os.path.basename(layer_dir)
    return re.fullmatch(f"{re.escape(usy)}_\\d+m", name) is not None


def _filter_altitude_layer_dirs(layer_dirs: list[str], usy: str) -> list[str]:
    return sorted(d for d in layer_dirs if is_standard_altitude_layer_dir(d, usy))


def _has_layer_subdirs(activespaces_dir: str) -> bool:
    if not os.path.isdir(activespaces_dir):
        return False
    return bool(glob.glob(_join(activespaces_dir, "*_*m")))


def resolve_nmsim_predictions_dir(site_dir: str, *, for_write: bool = False) -> str:
    """Return predictions cache dir; writes always use the new layout."""
    new_dir = _join(site_dir, NMSIM_PREDICTIONS_SUBDIR)
    if for_write:
        return new_dir
    if os.path.isdir(new_dir) or glob.glob(_join(new_dir, "*.csv")):
        return new_dir
    return _join(site_dir, LEGACY_PREDICTIONS_SUBDIR)


def resolve_nmsim_scratch_dir(site_dir: str, *, for_write: bool = False) -> str:
    """Return NMSim scratch dir (.tis); writes always use the new layout."""
    new_dir = _join(site_dir, NMSIM_SCRATCH_SUBDIR)
    if for_write:
        return new_dir
    if os.path.isdir(new_dir):
        return new_dir
    return _join(site_dir, LEGACY_PREDICTIONS_SUBDIR)


def resolve_nmsim_activespaces_dir(site_dir: str, *, for_write: bool = False) -> str:
    """Return ACTIVESPACES root; writes always use the new layout."""
    new_dir = _join(site_dir, NMSIM_ACTIVESPACES_SUBDIR)
    if for_write or _has_layer_subdirs(new_dir):
        return new_dir
    return _join(site_dir, LEGACY_ACTIVESPACES_SUBDIR)


def _altitude_m_from_layer_dir(layer_dir: str) -> int:
    name = os.path.basename(layer_dir)
    return int(name.rsplit("_", 1)[-1].removesuffix("m"))


def alternate_activespace_layer_dir(layer_dir: str) -> str | None:
    """Return the other NMSim ACTIVESPACES folder for the same altitude, if it exists."""
    layer_dir = os.path.abspath(layer_dir)
    name = os.path.basename(layer_dir)
    activespaces_dir = os.path.dirname(layer_dir)
    parent = os.path.dirname(activespaces_dir)
    if os.path.basename(parent) == "nmsim":
        other = _join(os.path.dirname(parent), "ACTIVESPACES", name)
        return other if os.path.isdir(other) else None
    if os.path.basename(parent) == "Output_Data":
        other = _join(parent, "nmsim", "ACTIVESPACES", name)
        return other if os.path.isdir(other) else None
    return None


def find_layer_geojson(layer_dir: str, gain: float) -> str | None:
    """Find ``*_O_+020.geojson`` in this layer dir, or the sibling NMSim layout."""
    sign = "-" if gain < 0 else "+"
    gain_string = str(abs(int(10 * gain))).zfill(3)
    pattern = f"*_O_{sign}{gain_string}.geojson"
    for directory in (layer_dir, alternate_activespace_layer_dir(layer_dir)):
        if not directory:
            continue
        matches = glob.glob(_join(directory, pattern))
        if matches:
            return matches[0]
    return None


def resolve_activespace_layer_dirs(
    project_dir: str,
    unit: str,
    site: str,
    year,
) -> list[str]:
    """Union new + legacy layer dirs; new layout wins when both exist for an altitude."""
    site_path = _join(project_dir, f"{unit}{site}")
    usy = _deployment_id(unit, site, year)
    pattern = f"{usy}_*m"
    by_alt: dict[int, str] = {}
    for match in _filter_altitude_layer_dirs(
        glob.glob(_join(site_path, LEGACY_ACTIVESPACES_SUBDIR, pattern)), usy,
    ):
        by_alt[_altitude_m_from_layer_dir(match)] = match
    for match in _filter_altitude_layer_dirs(
        glob.glob(_join(site_path, NMSIM_ACTIVESPACES_SUBDIR, pattern)), usy,
    ):
        by_alt[_altitude_m_from_layer_dir(match)] = match
    return [by_alt[alt] for alt in sorted(by_alt)]


def resolve_activespace_geojson(
    project_dir: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    gain_sign: str,
    gain_string: str,
) -> str:
    """Prefer new geojson path when the file exists; otherwise return legacy path."""
    site_path = _join(project_dir, f"{unit}{site}")
    usy = _deployment_id(unit, site, year)
    layer_name = f"{usy}_{altitude_m}m"
    filename = f"{usy}_O_{gain_sign}{gain_string}.geojson"
    new_path = _join(site_path, NMSIM_ACTIVESPACES_SUBDIR, layer_name, filename)
    if os.path.isfile(new_path):
        return new_path
    return _join(site_path, LEGACY_ACTIVESPACES_SUBDIR, layer_name, filename)
