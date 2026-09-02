#!/usr/bin/env python
"""Compare ELV pre-filter vs real AAM below-ground rejections on DENATRLA.

Builds a density-10 mesh (same helper as generate_active_space), classifies
each point with ``split_below_aam_terrain``, then runs unfiltered AAM n=1 on a
stratified sample. A 2-point hop probe checks interpolated-track terrain.

Run:
  docker/run_activespace.sh -m aam docker/validate_aam_below_ground.py

Requires staged ``vendor/aam-runtime/`` with ``FLATO200.nc``.
"""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from pyproj import Transformer
from shapely.geometry import Point

from aam_translator import hop_speed_kn, lonlat_to_model_ft, read_poi, read_run_log
from aam_translator.constants import FT_PER_M

from nps_active_space.active_space.aam_propagation_model import (
    AAM_INP_BASENAME,
    AamPropagationModel,
    _pad_single_point_track,
    aam_source_id_from_omni,
)
from nps_active_space.active_space.aam_terrain import (
    AAM_BELOW_SURFACE_TOLERANCE_M,
    _bilinear_sample_grid,
    _elv_grid_values,
    _terrain_surface_elevation_m,
    split_below_aam_terrain,
)
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.utils.computation import NMSIM_bbox_utm, build_src_point_mesh
from nps_active_space.utils.helpers import get_deployment

REPO = Path(os.environ.get("REPO_ROOT", "/repo"))
SITE_DIR = REPO / "example_data" / "site_projects" / "DENATRLA"
PROJECT_DIR = SITE_DIR.parent
AAM_SHIM = Path("/usr/local/bin/aam")
OMNI = Path(
    os.environ.get(
        "OMNI_SOURCE",
        str(REPO / "nps_active_space" / "data" / "tuning" / "O_+000.src"),
    )
)
MESH_DENSITY = int(os.environ.get("MESH_DENSITY", "10"))
ALTITUDE_M = int(os.environ.get("ALTITUDE_M", "1500"))
MAX_AAM_N1 = int(os.environ.get("MAX_AAM_N1", "40"))
N_BELOW_SAMPLE = 10
N_LOW_AGL_ABOVE = 15
N_HIGH_AGL_CONTROL = 5


def log(msg: str) -> None:
    print(f"[aam-below] {msg}", flush=True)


def _aeqd_from_utm(terrain, source_pts: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    to_aeqd = Transformer.from_crs(source_pts.crs, terrain.aeqd_crs, always_xy=True)
    xs = source_pts.geometry.x.to_numpy(dtype=np.float64)
    ys = source_pts.geometry.y.to_numpy(dtype=np.float64)
    if xs.size == 1:
        ax, ay = to_aeqd.transform(float(xs[0]), float(ys[0]))
        return np.asarray([ax], dtype=np.float64), np.asarray([ay], dtype=np.float64)
    ax, ay = to_aeqd.transform(xs, ys)
    return np.asarray(ax, dtype=np.float64), np.asarray(ay, dtype=np.float64)


def _aeqd_from_lonlat(terrain, wgs84_pts: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    spec = terrain.spec
    ax = np.empty(len(wgs84_pts), dtype=np.float64)
    ay = np.empty(len(wgs84_pts), dtype=np.float64)
    for i, geom in enumerate(wgs84_pts.geometry):
        x_ft, y_ft = lonlat_to_model_ft(terrain, float(geom.x), float(geom.y))
        ax[i] = spec.grid_origin_x_m + x_ft / FT_PER_M
        ay[i] = spec.grid_origin_y_m + y_ft / FT_PER_M
    return ax, ay


def _surface_from_lonlat_m(terrain, wgs84_pts: gpd.GeoDataFrame) -> np.ndarray:
    """Bilinear ELV sample at write_inp / lonlat_to_model_ft indices."""
    values = _elv_grid_values(Path(terrain.elv_path))
    spec = terrain.spec
    col_i = np.empty(len(wgs84_pts), dtype=np.float64)
    row_j = np.empty(len(wgs84_pts), dtype=np.float64)
    for i, geom in enumerate(wgs84_pts.geometry):
        x_ft, y_ft = lonlat_to_model_ft(terrain, float(geom.x), float(geom.y))
        col_i[i] = x_ft / (spec.cell_dx_m * FT_PER_M)
        row_j[i] = y_ft / (spec.cell_dy_m * FT_PER_M)
    raw = _bilinear_sample_grid(values, col_i, row_j)
    if terrain.elv_header_feet:
        return raw / FT_PER_M
    return raw


def _classify_aam_log(text: str, *, poi_rows: int, exc: str | None) -> str:
    lower = (text or "").lower()
    if "below ground" in lower:
        return "below_ground"
    if "read error" in lower:
        return "read_error"
    if "fpa" in lower or "forrtl" in lower:
        return "fpa"
    if poi_rows == 0 and exc:
        if "empty" in exc.lower() or "no data rows" in exc.lower():
            return "empty_poi"
        return "exception"
    if poi_rows == 0:
        return "empty_poi"
    return "ok"


def _excerpt(text: str, max_lines: int = 6) -> str:
    if not text:
        return ""
    keys = re.compile(r"below ground|terrain|read error|fpa|forrtl|error:", re.I)
    hits = [ln.strip() for ln in text.splitlines() if keys.search(ln)]
    if not hits:
        nonempty = [ln.strip() for ln in text.splitlines() if ln.strip()]
        hits = nonempty[-max_lines:]
    return " | ".join(hits[:max_lines])[:240]


def _run_aam_unfiltered(
    model: AamPropagationModel,
    site,
    source_pts: gpd.GeoDataFrame,
    omni: str,
    job_name: str,
) -> dict:
    """Stage + launch AAM without the ELV pre-filter."""
    work_dir = model._runs_dir / job_name
    work_dir.mkdir(parents=True, exist_ok=True)
    track = _pad_single_point_track(model._build_track(source_pts))
    pois = model._build_pois(site)
    source_id = aam_source_id_from_omni(omni)
    heading_deg = 0.0
    speed_kn = hop_speed_kn(track, site.terrain)
    inp_path = work_dir / f"{AAM_INP_BASENAME}.inp"
    log_path = work_dir / f"{AAM_INP_BASENAME}.txt"
    poi_path = work_dir / f"{AAM_INP_BASENAME}.POI"
    exc_msg = None
    poi_rows = 0
    try:
        model._stage_run_dir(
            work_dir, site, track, pois, source_id, job_name, heading_deg, speed_kn,
        )
        model._run_aam(inp_path, work_dir)
        if poi_path.is_file():
            histories = read_poi(poi_path)
            poi_rows = 0 if not histories else histories[0].n_samples
        run_log = read_run_log(log_path) if log_path.is_file() else None
        if run_log is not None and not run_log.ok:
            raise RuntimeError(f"AAM run failed: {run_log.read_error}")
        if poi_rows == 0:
            raise RuntimeError(f"empty .POI file: {poi_path}")
        ok = True
    except Exception as exc:
        ok = False
        exc_msg = str(exc)

    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
    stderr_path = work_dir / "aam_stderr.txt"
    if stderr_path.is_file():
        log_text = log_text + "\n" + stderr_path.read_text(encoding="utf-8", errors="replace")
    reason = _classify_aam_log(log_text, poi_rows=poi_rows, exc=exc_msg)
    return {
        "ok": ok,
        "reason": reason if not ok else "ok",
        "poi_rows": poi_rows,
        "excerpt": _excerpt(log_text) or (exc_msg or "")[:240],
        "work_dir": str(work_dir),
    }


def _dem_stats(dem_path: Path) -> tuple[float, float]:
    with rasterio.open(dem_path) as src:
        arr = src.read(1)
        nodata = src.nodata
        if nodata is not None:
            arr = arr[arr != nodata]
        return float(np.nanmin(arr)), float(np.nanmax(arr))


def _build_mesh(study_area: gpd.GeoDataFrame, altitude_m: int) -> gpd.GeoDataFrame:
    crs = NMSIM_bbox_utm(study_area.iloc[[0]])
    area_utm = study_area.to_crs(crs)
    mesh = build_src_point_mesh(area_utm, MESH_DENSITY, altitude_m)
    inside = mesh[mesh.within(area_utm.union_all())].copy()
    inside["pt_id"] = np.arange(len(inside))
    return inside


def _annotate_filter(mesh: gpd.GeoDataFrame, terrain) -> gpd.GeoDataFrame:
    above, below = split_below_aam_terrain(terrain, mesh, job_name="diag")
    out = mesh.copy()
    surface = _terrain_surface_elevation_m(mesh, terrain)
    wgs = mesh.to_crs("EPSG:4326")
    surface_inp = _surface_from_lonlat_m(terrain, wgs)
    ax_f, ay_f = _aeqd_from_utm(terrain, mesh)
    ax_i, ay_i = _aeqd_from_lonlat(terrain, wgs)
    z = mesh.geometry.z.to_numpy(dtype=np.float64)
    out["surface_filter_m"] = surface
    out["surface_inp_m"] = surface_inp
    out["agl_filter_m"] = z - surface
    out["agl_inp_m"] = z - surface_inp
    out["dxy_m"] = np.hypot(ax_f - ax_i, ay_f - ay_i)
    out["dsurf_m"] = surface_inp - surface
    out["filter_below"] = False
    if len(below) > 0:
        out.loc[below.index, "filter_below"] = True
    out["filter_nan"] = np.isnan(surface)
    return out


def _select_n1_sample(classified: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    below = classified[classified["filter_below"] | classified["filter_nan"]]
    above = classified[~classified["filter_below"] & ~classified["filter_nan"]]
    parts: list[gpd.GeoDataFrame] = []
    if len(below) > 0:
        parts.append(below.nsmallest(min(N_BELOW_SAMPLE, len(below)), "agl_filter_m"))
    if len(above) > 0:
        parts.append(above.nsmallest(min(N_LOW_AGL_ABOVE, len(above)), "agl_filter_m"))
        parts.append(above.nlargest(min(N_HIGH_AGL_CONTROL, len(above)), "agl_filter_m"))
    if not parts:
        return classified.iloc[0:0]
    sample = gpd.GeoDataFrame(pd_concat_unique(parts), crs=classified.crs)
    if len(sample) > MAX_AAM_N1:
        sample = sample.iloc[:MAX_AAM_N1]
    return sample


def pd_concat_unique(frames: list[gpd.GeoDataFrame]) -> gpd.GeoDataFrame:
    import pandas as pd

    combined = pd.concat(frames)
    return combined[~combined.index.duplicated(keep="first")]


def _print_row(row, aam: dict | None) -> None:
    filt = "below" if row.filter_below else ("nan" if row.filter_nan else "above")
    aam_s = "skip"
    reason = ""
    excerpt = ""
    if aam is not None:
        aam_s = "ok" if aam["ok"] else "FAIL"
        reason = aam["reason"]
        excerpt = aam["excerpt"]
    log(
        f"  id={int(row.pt_id):3d}  filter={filt:5s}  agl_f={row.agl_filter_m:8.2f}m  "
        f"agl_inp={row.agl_inp_m:8.2f}m  dxy={row.dxy_m:6.3f}m  dsurf={row.dsurf_m:7.3f}m  "
        f"aam={aam_s:4s}  {reason}  {excerpt}"
    )


def main() -> int:
    if not AAM_SHIM.is_file():
        log(f"ERROR: shim missing at {AAM_SHIM}")
        return 1
    if not SITE_DIR.is_dir():
        log(f"ERROR: site dir missing: {SITE_DIR}")
        return 1

    dem_tif, _ = get_project_setup_elevation(SITE_DIR)
    dem_min, dem_max = _dem_stats(dem_tif)
    log(f"DEM {dem_tif.name}: min={dem_min:.1f}m max={dem_max:.1f}m")

    study_area = gpd.read_file(SITE_DIR / "DENATRLA_study_area.shp")
    mic = get_deployment(str(PROJECT_DIR), "DENA", "TRLA", 2025, elevation=False)
    log(f"mic {mic.name} lon={mic.lon:.5f} lat={mic.lat:.5f} z_agl={mic.z:.2f}m")

    altitudes = [ALTITUDE_M]
    if not (dem_min < ALTITUDE_M < dem_max):
        mid = int(round((dem_min + dem_max) / 2))
        log(f"altitude {ALTITUDE_M}m is outside DEM range; also probing {mid}m")
        altitudes.append(mid)

    with tempfile.TemporaryDirectory(prefix="aam_below_ground_") as tmp:
        root = Path(tmp) / "DENATRLA"
        (root / "Input_Data").mkdir(parents=True)
        model = AamPropagationModel(str(root), aam_shim=str(AAM_SHIM))
        site = model.prepare_site(
            str(dem_tif),
            study_area.iloc[[0]],
            mic,
            project_dem=False,
            suffix=f"_{mic.name}",
        )
        log(
            f"ELV {site.terrain.elv_path} header_feet={site.terrain.elv_header_feet} "
            f"tol={AAM_BELOW_SURFACE_TOLERANCE_M}m"
        )

        for altitude_m in altitudes:
            mesh = _build_mesh(study_area, altitude_m)
            log(f"--- mesh density={MESH_DENSITY} altitude={altitude_m}m n={len(mesh)} ---")
            classified = _annotate_filter(mesh, site.terrain)
            n_below = int(classified["filter_below"].sum())
            n_nan = int(classified["filter_nan"].sum())
            n_above = len(classified) - n_below
            log(
                f"filter: above={n_above} below={n_below} nan={n_nan}  "
                f"median |dxy|={classified['dxy_m'].median():.4f}m  "
                f"median |dsurf|={classified['dsurf_m'].abs().median():.4f}m  "
                f"max |dsurf|={classified['dsurf_m'].abs().max():.4f}m"
            )

            sample = _select_n1_sample(classified)
            log(f"AAM n=1 sample: {len(sample)} points")
            aam_by_id: dict[int, dict] = {}
            for _, row in sample.iterrows():
                pts = classified.loc[[row.name]]
                job = f"n1_a{altitude_m}_p{int(row.pt_id)}"
                result = _run_aam_unfiltered(model, site, pts, str(OMNI), job)
                aam_by_id[int(row.pt_id)] = result
                _print_row(row, result)

            above_fail = []
            below_ok = []
            for _, row in sample.iterrows():
                res = aam_by_id[int(row.pt_id)]
                if (not row.filter_below) and (not row.filter_nan) and (not res["ok"]):
                    above_fail.append((row, res))
                if (row.filter_below or row.filter_nan) and res["ok"]:
                    below_ok.append((row, res))

            log(
                f"mismatches: filter-above/AAM-fail={len(above_fail)}  "
                f"filter-below/AAM-ok={len(below_ok)}"
            )
            if above_fail:
                log("worst filter-above/AAM-fail (smallest AGL):")
                above_fail.sort(key=lambda t: t[0].agl_filter_m)
                for row, res in above_fail[:3]:
                    log(
                        f"  id={int(row.pt_id)} agl_f={row.agl_filter_m:.2f}m "
                        f"agl_inp={row.agl_inp_m:.2f}m dxy={row.dxy_m:.3f}m "
                        f"reason={res['reason']} dir={res['work_dir']}"
                    )
                    log(f"    excerpt: {res['excerpt']}")

            successes = [
                classified.loc[idx]
                for idx in sample.index
                if aam_by_id[int(classified.loc[idx].pt_id)]["ok"]
                and not classified.loc[idx].filter_below
            ]
            if len(successes) >= 2:
                a = successes[0]
                rest = successes[1:]
                dists = [
                    Point(a.geometry.x, a.geometry.y).distance(Point(b.geometry.x, b.geometry.y))
                    for b in rest
                ]
                b = rest[int(np.argmin(dists))]
                pair = classified.loc[[a.name, b.name]]
                hop = _run_aam_unfiltered(
                    model, site, pair, str(OMNI), f"hop_a{altitude_m}_{int(a.pt_id)}_{int(b.pt_id)}",
                )
                log(
                    f"hop probe: n=1 ids {int(a.pt_id)}+{int(b.pt_id)} both ok; "
                    f"2-pt track {'ok' if hop['ok'] else 'FAIL'} ({hop['reason']})  {hop['excerpt']}"
                )
            else:
                log("hop probe: skipped (need two AAM-ok filter-above points)")

            if n_below > 0 or altitude_m != altitudes[0]:
                break
            if n_below == 0 and len(altitudes) == 1:
                mid = int(round((dem_min + dem_max) / 2))
                if mid != altitude_m:
                    log(f"no filter-below at {altitude_m}m; adding altitude {mid}m")
                    altitudes.append(mid)

    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
