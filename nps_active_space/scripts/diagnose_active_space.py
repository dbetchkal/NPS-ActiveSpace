"""
Summarize site inputs and recent outputs to explain empty ACTIVESPACES geojsons.

Run from the repository root:

    python -m nps_active_space.scripts.diagnose_active_space -e DENA_example -u DENA -s TRLA -y 2025
"""

from __future__ import annotations

import argparse
import glob
import pickle
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

import nps_active_space.utils.config as cfg
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.utils.computation import NMSIM_bbox_utm, build_src_point_mesh
from nps_active_space.utils.helpers import get_deployment, get_omni_sources, load_studyarea


def _sample_elevation_m(tif_path: Path, xs: pd.Series, ys: pd.Series, crs: str) -> np.ndarray:
    with rasterio.open(tif_path) as dem:
        proj = Transformer.from_crs(crs, dem.crs, always_xy=True)
        px, py = proj.transform(xs, ys)
        pts = np.stack([px, py], axis=1)
        samples = [s[0] for s in dem.sample(pts)]
        nodata = dem.nodata
        return np.array(
            [np.nan if v is None or v == nodata else float(v) for v in samples],
            dtype=float,
        )


def _summarize_tested_points(site_dir: Path, altitude_m: int | None) -> None:
    pattern = str(site_dir / "Output_Data/TESTED_POINTS" / "**" / "*.pkl")
    pkls = sorted(glob.glob(pattern, recursive=True))
    if not pkls:
        print("  No TESTED_POINTS .pkl files found.")
        return

    print(f"  Found {len(pkls)} tested-points pickle(s). Summarizing latest per omni stem:")
    by_stem: dict[str, str] = {}
    for path in pkls:
        stem = Path(path).stem
        by_stem[stem] = path

    for stem, path in sorted(by_stem.items())[:5]:
        with open(path, "rb") as f:
            tested_by_heading = pickle.load(f)
        total = 0
        audible = 0
        for heading_pts in tested_by_heading.values():
            total += len(heading_pts)
            audible += int((heading_pts["audible"] == 1).sum())
        pct = 100.0 * audible / total if total else 0.0
        print(f"    {stem}: {audible}/{total} audible ({pct:.1f}%)")
    if len(by_stem) > 5:
        print(f"    ... and {len(by_stem) - 5} more")


def _summarize_tig_tis(site_dir: Path) -> None:
    tis_dir = site_dir / "Output_Data/TIG_TIS"
    tis_files = sorted(tis_dir.glob("*.tis"))
    csv_files = sorted(tis_dir.glob("*.csv"))
    print(f"  TIG_TIS: {len(tis_files)} .tis, {len(csv_files)} prediction .csv")
    if csv_files:
        latest = csv_files[-1]
        df = pd.read_csv(latest)
        if "A" in df.columns:
            a_col = df["A"].astype(float) / 10.0
            print(
                f"    Latest csv {latest.name}: n={len(df)}, "
                f"A dB min={a_col.min():.1f} max={a_col.max():.1f} mean={a_col.mean():.1f}"
            )


def _summarize_geojsons(site_dir: Path) -> None:
    geojsons = sorted((site_dir / "Output_Data/ACTIVESPACES").glob("**/*.geojson"))
    if not geojsons:
        print("  No ACTIVESPACES geojsons found.")
        return
    nonempty = 0
    for path in geojsons[:10]:
        gdf = gpd.read_file(path)
        has_geom = (
            not gdf.empty
            and gdf.geometry.notna().any()
            and (~gdf.geometry.is_empty).any()
        )
        if has_geom:
            nonempty += 1
        print(f"    {path.name}: rows={len(gdf)} nonempty={has_geom}")
    if len(geojsons) > 10:
        print(f"    ... {len(geojsons)} total geojsons, {nonempty} nonempty in first 10")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-e", "--environment", required=True)
    parser.add_argument("-u", "--unit", required=True)
    parser.add_argument("-s", "--site", required=True)
    parser.add_argument("-y", "--year", type=int, required=True)
    parser.add_argument("-l", "--altitude", type=int, default=None, help="Source altitude (m) to probe.")
    args = parser.parse_args()

    cfg.initialize(environment=args.environment)
    project_dir = Path(cfg.read("project", "dir"))
    site_dir = project_dir / f"{args.unit}{args.site}"
    usy = f"{args.unit}{args.site}{args.year}"

    print(f"=== Active space diagnostics: {usy} ===")
    print(f"Site dir: {site_dir}")

    nmsim = cfg.read("project", "nmsim")
    print(f"NMSIM binary: {nmsim!r} {'(MISSING/empty)' if not nmsim or not Path(nmsim).is_file() else '(ok)'}")

    tif_path: Path | None = None
    try:
        tif_path, flt_path = get_project_setup_elevation(site_dir)
        print(f"Elevation tif: {tif_path.name} ({tif_path.stat().st_size // 1024} KB)")
        print(f"Elevation flt: {flt_path.name} ({flt_path.stat().st_size // 1024} KB)")
        if not flt_path.is_file():
            print("  ERROR: .flt missing — NMSIM cannot run propagation.")
    except FileNotFoundError as exc:
        print(f"  ERROR: {exc}")

    sit_path = site_dir / "Input_Data/05_SITES" / f"{usy}.sit"
    if sit_path.is_file():
        sit_lines = sit_path.read_text().splitlines()
        print(f".sit listener: {sit_lines[2].strip() if len(sit_lines) > 2 else '?'}"
              f"\n.sit dem path: {sit_lines[3].strip() if len(sit_lines) > 3 else '?'}")
        dem_ref = sit_lines[3].strip() if len(sit_lines) > 3 else ""
        if dem_ref:
            candidates = [
                Path(dem_ref),
                site_dir / dem_ref,
                sit_path.parent / dem_ref,
            ]
            if not any(p.is_file() for p in candidates):
                print("  WARNING: .sit dem path does not resolve to an existing file.")

    study_area = load_studyarea(str(project_dir), args.unit, args.site, args.year)
    mic = get_deployment(str(project_dir), args.unit, args.site, args.year, elevation=False)
    crs = NMSIM_bbox_utm(study_area)
    mic_proj = mic.to_crs(crs)
    print(f"Mic UTM: E={mic_proj.x:.1f} N={mic_proj.y:.1f} z_agl={mic_proj.z:.2f} m crs={crs}")

    if tif_path is not None:
        altitude_m = args.altitude or 3000
        grid = build_src_point_mesh(study_area.to_crs(crs), 48, altitude_m)
        elevs = _sample_elevation_m(tif_path, grid.geometry.x, grid.geometry.y, crs)
        valid_elevs = np.isfinite(elevs)
        underground = int((grid.geometry.z.values[valid_elevs] < elevs[valid_elevs]).sum()) if valid_elevs.any() else 0
        nodata = int((~valid_elevs).sum())
        above = len(grid) - underground - nodata
        print(
            f"Source grid @ {altitude_m}m: {len(grid)} pts, "
            f"aboveground={above}, underground={underground}, dem_nodata={nodata}"
        )
        if valid_elevs.any():
            print(
                f"  DEM ground sample: min={np.nanmin(elevs):.0f} max={np.nanmax(elevs):.0f} "
                f"mean={np.nanmean(elevs):.0f} m"
            )
        if above == 0:
            print("  ERROR: All source points classified underground — geojsons will be empty.")
            print("  Check elevation units (meter tif) and source altitude vs terrain.")

    omni_sources = get_omni_sources(lower=-10, upper=40)
    print(f"Omni sources to run: {len(omni_sources)} (default -10..40 dB range)")

    print("\n--- Outputs ---")
    _summarize_geojsons(site_dir)
    _summarize_tested_points(site_dir, args.altitude)
    _summarize_tig_tis(site_dir)

    print("\n--- Likely causes of empty geojsons ---")
    print("  1. All NMSIM-tested points inaudible (ambience threshold too high vs omni gain).")
    print("  2. All source points underground (wrong DEM units or altitude far below terrain).")
    print("  3. NMSIM ran but produced empty/invalid .tis (binary path, permissions, or too many points).")
    print("  4. Contouring found no audible/inaudible boundary (zero audible points in mesh).")
    print("  Inspect TESTED_POINTS pickles: audible count should be > 0 for at least one omni.")


if __name__ == "__main__":
    main()
