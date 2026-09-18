#!/usr/bin/env python
"""Sea-level source mesh + one propagation predict (NMSim or AAM).

Exercises DEM/ELV clearance and a minimal ``_run_propagation_model`` batch.
Run inside the container, e.g.:

  docker/run_activespace.sh docker/validate_sea_level_sources.py --model nmsim
  docker/run_activespace.sh -m aam docker/validate_sea_level_sources.py --model aam
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from shapely.geometry import box

import nps_active_space.utils.config as cfg
from nps_active_space.active_space.active_space_setup import (
    build_active_space_generator,
)
from nps_active_space.propagation_model.protocol import THIRD_OCTAVE_BANDS
from nps_active_space.utils.computation import build_src_point_mesh
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.setup.elevation import (
    GRIDFLOAT_NODATA,
    NODATA_INT16,
    get_project_setup_elevation,
    write_gridfloat,
)
from nps_active_space.utils.helpers import get_deployment, get_omni_sources


def _ensure_nmsim_gridfloat(site_dir: Path) -> None:
    """Create .flt/.hdr siblings from an existing project_setup GeoTIFF if missing."""
    try:
        get_project_setup_elevation(site_dir)
        return
    except FileNotFoundError:
        elev_dir = site_dir / "Input_Data" / "01_ELEVATION"
        tif_paths = sorted(elev_dir.glob("elevation_m_nad83_utm*.tif"))
        if not tif_paths:
            raise
        tif_path = tif_paths[0]
        with rasterio.open(tif_path) as dem:
            band = dem.read(1)
            write_gridfloat(
                tif_path.with_suffix(""),
                band.astype(np.int16),
                dem.transform,
                dem.width,
                dem.height,
                NODATA_INT16,
                GRIDFLOAT_NODATA,
            )
        _log(f"wrote GridFloat siblings for {tif_path.name}")


def _log(msg: str) -> None:
    print(f"[sea-level-sources] {msg}", flush=True)


def _mic_clip(study: gpd.GeoDataFrame, mic, half_width_m: float = 2000.0) -> gpd.GeoDataFrame:
    crs = study.crs
    to_utm = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x_m, y_m = to_utm.transform(mic.lon, mic.lat)
    window = box(
        x_m - half_width_m,
        y_m - half_width_m,
        x_m + half_width_m,
        y_m + half_width_m,
    )
    clip = study.clip(gpd.GeoDataFrame(geometry=[window], crs=crs))
    if clip.empty:
        raise RuntimeError("mic-centered clip did not intersect study area")
    return clip


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        type=AcousticModel,
        choices=list(AcousticModel),
        default=None,
        help="Propagation model (default: ACOUSTIC_MODEL env or nmsim)",
    )
    parser.add_argument("-e", "--environment", default="container")
    parser.add_argument("-u", "--unit", default="GLBA")
    parser.add_argument("-s", "--site", default="LSTL")
    parser.add_argument("-y", "--year", type=int, default=2024)
    args = parser.parse_args()

    model = args.model or AcousticModel.parse(os.environ.get("ACOUSTIC_MODEL", "nmsim"))
    repo = Path(__file__).resolve().parents[1]
    os.chdir(repo)

    cfg.initialize(args.environment)
    project_dir = cfg.read("project", "dir")
    site_dir = Path(project_dir) / f"{args.unit}{args.site}"
    study_path = site_dir / f"{args.unit}{args.site}_study_area.shp"
    if not study_path.is_file():
        _log(f"ERROR: missing study area {study_path}")
        return 1

    _ensure_nmsim_gridfloat(site_dir)

    mic = get_deployment(project_dir, args.unit, args.site, args.year, elevation=False)
    study = gpd.read_file(study_path)
    clip = _mic_clip(study.to_crs(study.estimate_utm_crs()), mic)

    ambience = pd.Series({str(band): 45.0 for band in THIRD_OCTAVE_BANDS})
    generator = build_active_space_generator(
        str(site_dir),
        clip,
        ambience,
        mic,
        model,
    )

    mesh = build_src_point_mesh(clip, density=6, altitude=0)
    if mesh.empty:
        _log("ERROR: empty source mesh")
        return 1

    above, below = generator._determine_underground_pts(mesh)
    _log(f"DEM clearance: {len(above)} above, {len(below)} underground (mesh n={len(mesh)})")
    if above.empty:
        _log("ERROR: no above-ground points after DEM clearance")
        return 1

    sample = above.head(3)
    omni_sources = get_omni_sources(lower=0, upper=0, step_db=5)
    if not omni_sources:
        _log("ERROR: no omni source on ladder at 0 dB")
        return 1
    omni_source = omni_sources[0]
    _log(f"predicting {len(sample)} points with {model.value} omni={Path(omni_source).name}")

    audibility_pts = generator._run_propagation_model(
        "sea_level_smoke",
        sample,
        omni_source,
        altitude_m=0,
        heading=0,
    )
    n_audible = int(audibility_pts["audible"].sum()) if not audibility_pts.empty else 0
    _log(f"audibility rows={len(audibility_pts)} audible={n_audible}")
    if audibility_pts.empty:
        _log("ERROR: propagation returned no audibility rows")
        return 1
    if audibility_pts["audible"].isna().any():
        _log("ERROR: audibility column has NaN")
        return 1

    _log(f"OK: {model.value} sea-level source predict completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
