"""Legacy NMSim output path resolvers (read-side backward compatibility only)."""

from __future__ import annotations

import glob
import os
import re

LEGACY_PREDICTIONS_SUBDIR = "Output_Data/TIG_TIS"
LEGACY_ACTIVESPACES_SUBDIR = "Output_Data/ACTIVESPACES"

NMSIM_OUTPUT_SUBDIR = "Output_Data/nmsim"
NMSIM_PREDICTIONS_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/predictions"
NMSIM_SCRATCH_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/scratch"
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


def resolve_activespace_layer_dirs(
    project_dir: str,
    unit: str,
    site: str,
    year,
) -> list[str]:
    """Glob layer dirs under new layout first, then legacy ACTIVESPACES."""
    site_path = _join(project_dir, f"{unit}{site}")
    usy = _deployment_id(unit, site, year)
    pattern = f"{usy}_*m"
    new_matches = glob.glob(_join(site_path, NMSIM_ACTIVESPACES_SUBDIR, pattern))
    if new_matches:
        return _filter_altitude_layer_dirs(new_matches, usy)
    legacy_matches = glob.glob(_join(site_path, LEGACY_ACTIVESPACES_SUBDIR, pattern))
    return _filter_altitude_layer_dirs(legacy_matches, usy)


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
