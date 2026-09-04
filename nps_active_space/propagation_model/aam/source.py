"""Resolve NMSim omni tunings and native NetCDF into cached AAM source files."""

from __future__ import annotations

import math
import re
import shutil
from pathlib import Path

from aam_translator.bands import band_number_for_frequency

AAM_TEMPLATE_NC_FILENAME = "OMNI_200.nc"
AAM_NC_CACHE_REL = Path("Input_Data/aam/NCfiles")
OMNI_STEM_RE = re.compile(r"^O_([+-])(\d{3})$")
AAM_TOKEN_RE = re.compile(r"^[A-Z0-9_]{5}\d+$")

FT_PER_M = 3.280839895
AAM_DEFAULT_REF_M = 304.8
DEFAULT_RADIUS_FT = 1000.0

__all__ = [
    "AAM_NC_CACHE_REL",
    "AAM_TEMPLATE_NC_FILENAME",
    "aam_source_id_from_omni",
    "ensure_aam_nc_for_source",
    "omni_stem_to_aam_token",
    "read_avg_spectrum_db",
    "site_ncfiles_dir",
    "stage_run_ncfiles",
    "write_aam_nc",
]


def site_ncfiles_dir(root_dir: str | Path) -> Path:
    return Path(root_dir).resolve() / AAM_NC_CACHE_REL


def omni_stem_to_aam_token(stem: str) -> str:
    """Map NMSim omni stem ``O_±XXX`` to an AAM 5-char-prefix token."""
    match = OMNI_STEM_RE.fullmatch(stem)
    if match is None:
        raise ValueError(f"not an NMSim omni stem: {stem!r}")
    sign, digits = match.groups()
    if sign == "+":
        return f"OMNI_{digits}"
    return f"OMNIM{digits}"


def aam_source_id_from_omni(omni_source: str) -> str:
    """Return the AAM ``source_id`` token for a tuning path without writing NetCDF."""
    path = Path(omni_source)
    stem = path.stem
    if path.suffix.lower() == ".nc":
        if not AAM_TOKEN_RE.fullmatch(stem):
            raise ValueError(
                f"AAM NetCDF stem must be 5-char prefix + digits, got {stem!r}",
            )
        return stem
    if stem.startswith("OMNI") and AAM_TOKEN_RE.fullmatch(stem):
        return stem
    return omni_stem_to_aam_token(stem)


def read_avg_template(path: str | Path) -> tuple[list[float], float, list[float]]:
    """Parse an NPS ``.avg`` header -> (theta list, phi value, band freqs)."""
    with open(path, "r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    header = lines[0]
    band_freqs = [float(x) for x in header.split()[2:]]
    thetas: list[float] = []
    phi_val = 78.0
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 2 + len(band_freqs):
            continue
        try:
            thetas.append(float(parts[0]))
            phi_val = float(parts[1])
        except ValueError:
            continue
    return thetas, phi_val, band_freqs


def read_avg_spectrum_db(avg_path: str | Path) -> dict[int, float]:
    """Read the first data row of an ``.avg`` file as ANSI band -> dB."""
    _, _, band_freqs = read_avg_template(avg_path)
    levels_cb: list[int] | None = None
    with open(avg_path, errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2 + len(band_freqs):
                continue
            try:
                float(parts[0])
                float(parts[1])
                levels_cb = [int(parts[2 + i]) for i in range(len(band_freqs))]
                break
            except ValueError:
                continue
    if levels_cb is None:
        raise ValueError(f"no spectrum row found in {avg_path}")
    levels_db = [cb / 10.0 for cb in levels_cb]
    return {
        band_number_for_frequency(freq): level_db
        for freq, level_db in zip(band_freqs, levels_db, strict=True)
    }


def write_aam_nc(
    levels_db: dict[int, float],
    template_nc: str | Path,
    out_nc: str | Path,
    *,
    radius_ft: float = DEFAULT_RADIUS_FT,
) -> Path:
    """Clone ``template_nc`` and overwrite ``AMPLITUDE`` (omni) and ``RADIUS``."""
    import netCDF4 as nc
    import numpy as np

    template_nc = Path(template_nc)
    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(template_nc) as tmpl:
        freq = np.asarray(tmpl.variables["FREQUENCY"][:], dtype=float)
        amp_shape = tmpl.variables["AMPLITUDE"].shape
        per_band = np.array(
            [_level_at_band(levels_db, f) for f in freq],
            dtype=float,
        )
        amplitude = np.empty(amp_shape, dtype=float)
        amplitude[:, :, :] = per_band[np.newaxis, np.newaxis, :]

        dst = nc.Dataset(out_nc, "w", format=tmpl.file_format)
        try:
            dst.setncatts({k: tmpl.getncattr(k) for k in tmpl.ncattrs()})
            for dname, dim in tmpl.dimensions.items():
                dst.createDimension(dname, None if dim.isunlimited() else len(dim))
            for vname, var in tmpl.variables.items():
                out = dst.createVariable(vname, var.datatype, var.dimensions)
                out.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
                if vname == "AMPLITUDE":
                    out[:] = amplitude
                elif vname == "RADIUS":
                    out[:] = radius_ft
                else:
                    out[:] = var[:]
        finally:
            dst.close()
    return out_nc


def ensure_aam_nc_for_source(
    source_path: str | Path,
    root_dir: str | Path,
    template_nc: str | Path,
) -> tuple[str, Path]:
    """Stage or generate the NetCDF for ``source_path`` under the site cache."""
    source_path = Path(source_path).resolve()
    template_nc = Path(template_nc)
    if not template_nc.is_file():
        raise FileNotFoundError(
            f"AAM template NetCDF not found: {template_nc}",
        )

    cache_dir = site_ncfiles_dir(root_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    suffix = source_path.suffix.lower()
    if suffix == ".nc":
        token = aam_source_id_from_omni(str(source_path))
        cached = cache_dir / f"{token}.nc"
        if _needs_refresh(source_path, cached):
            shutil.copy2(source_path, cached)
        return token, cached

    if suffix == ".src":
        avg_path = source_path.with_suffix(".avg")
        if not avg_path.is_file():
            raise FileNotFoundError(
                f"AAM source {source_path} needs sibling {avg_path.name} or a .nc path",
            )
        token = aam_source_id_from_omni(str(source_path))
        cached = cache_dir / f"{token}.nc"
        if _needs_refresh(avg_path, cached):
            levels_db = read_avg_spectrum_db(avg_path)
            write_aam_nc(levels_db, template_nc, cached)
        return token, cached

    raise ValueError(
        f"AAM source must be .src (with sibling .avg) or .nc, got {source_path}",
    )


def stage_run_ncfiles(work_dir: str | Path, omni_nc: str | Path) -> Path:
    """Stage one omni NetCDF under ``work_dir/NCfiles`` for a single AAM subprocess.

    AAM loads every sphere under ``ROTOR_NOISE``. Per-gain files cloned from the
    vendor template share PROFILE metadata, so a shared site cache with multiple
    ``OMNI_*.nc`` files triggers "Same Profile data found in multiple spheres".
    """
    work_dir = Path(work_dir)
    omni_nc = Path(omni_nc)
    run_nc_dir = work_dir / "NCfiles"
    run_nc_dir.mkdir(parents=True, exist_ok=True)
    dest = run_nc_dir / omni_nc.name
    shutil.copy2(omni_nc, dest)
    return run_nc_dir


def _needs_refresh(source: Path, cached: Path) -> bool:
    if not cached.is_file():
        return True
    return source.stat().st_mtime > cached.stat().st_mtime


def _level_at_band(levels_db: dict[int, float], freq_hz: float) -> float:
    band = band_number_for_frequency(freq_hz)
    if band in levels_db:
        return levels_db[band]
    nearest = min(levels_db, key=lambda b: abs(b - band))
    return levels_db[nearest]
