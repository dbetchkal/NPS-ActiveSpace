#!/usr/bin/env python
"""Experiment: can two AAM Wine processes run concurrently?

generate_active_space.py uses mp.Pool over omni workers; each job gets its own
cwd / run dir while NCfiles and the Wine prefix are shared.

Uses the two-point ridge COMPUTEPOI adapter (same path as generate), not
noisecon.inp — that deck needs CH146.nc which is not in the currently staged
NCfiles (FLATO200 / OMNI_200 only).

Run:
  docker/run_activespace.sh -m aam docker/validate_aam_parallel.py
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import tempfile
import time
from pathlib import Path

import geopandas as gpd
from rasterio import open as rio_open
from shapely.geometry import Point, box

from nps_active_space.propagation_model.aam.model import AamPropagationModel
from nps_active_space.utils.models import Microphone

AAM_SHIM = Path("/usr/local/bin/aam")
FIXTURE_DIR = Path(
    os.environ.get(
        "AAM_FIXTURE_DIR",
        "/repo/tests/active_space/fixtures/two_point_ridge",
    ),
)
N_WORKERS = int(os.environ.get("AAM_PARALLEL_N", "2"))
DBA_TOL = 0.05


def log(msg: str) -> None:
    print(f"[aam-parallel] {msg}", flush=True)


_WORKER_STATE: dict = {}


def _build_adapter_site(root: Path) -> tuple[AamPropagationModel, object, gpd.GeoDataFrame, str]:
    meta = json.loads((FIXTURE_DIR / "case_meta.json").read_text())
    dem_path = FIXTURE_DIR / "parent_dem_utm.tif"
    (root / "Input_Data").mkdir(parents=True)

    rx_lon, rx_lat = meta["receiver_lonlat"]
    mic = Microphone(name="Receiver", lat=rx_lat, lon=rx_lon, z=4.92)
    rows = meta["source_points_utm"]
    crs = "EPSG:32606"
    source_pts = gpd.GeoDataFrame(
        {"label": [r["label"] for r in rows]},
        geometry=[Point(r["x"], r["y"], r["z"]) for r in rows],
        crs=crs,
    )
    with rio_open(dem_path) as ds:
        bounds = ds.bounds
    aoi = box(bounds.left, bounds.bottom, bounds.right, bounds.top)

    model = AamPropagationModel(str(root), aam_shim=str(AAM_SHIM))
    site = model.prepare_site(
        str(dem_path),
        gpd.GeoDataFrame(geometry=[aoi], crs=crs),
        mic,
        project_dem=False,
    )
    omni = os.environ.get(
        "OMNI_SOURCE",
        "/repo/nps_active_space/propagation_model/nmsim/data/tuning/O_+200.avg",
    )
    return model, site, source_pts, omni


def _adapter_predict(job_name: str) -> dict:
    model = _WORKER_STATE["model"]
    site = _WORKER_STATE["site"]
    source_pts = _WORKER_STATE["source_pts"]
    omni = _WORKER_STATE["omni"]
    start = time.perf_counter()
    preds = model.predict(
        site, source_pts, omni, altitude_m=400, job_name=job_name, heading=90,
    )
    return {
        "job_name": job_name,
        "A": [float(v) for v in preds["A"].tolist()],
        "n": len(preds),
        "elapsed_s": time.perf_counter() - start,
        "error": None,
    }


def _adapter_predict_guarded(job_name: str) -> dict:
    try:
        return _adapter_predict(job_name)
    except Exception as exc:
        return {
            "job_name": job_name,
            "A": [],
            "n": 0,
            "elapsed_s": 0.0,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _load_worker_state(root: Path) -> None:
    model, site, source_pts, omni = _build_adapter_site(root)
    _WORKER_STATE["model"] = model
    _WORKER_STATE["site"] = site
    _WORKER_STATE["source_pts"] = source_pts
    _WORKER_STATE["omni"] = omni


def run_adapter_phase(n: int) -> int:
    log(f"--- adapter COMPUTEPOI (two-point ridge, n={n}) ---")
    with tempfile.TemporaryDirectory(prefix="aam_par_adapter_") as tmp:
        seq_root = Path(tmp) / "seq"
        par_root = Path(tmp) / "par"

        _load_worker_state(seq_root)
        start = time.perf_counter()
        seq = [_adapter_predict_guarded(f"seq_{i}") for i in range(n)]
        seq_wall = time.perf_counter() - start

        _load_worker_state(par_root)
        start = time.perf_counter()
        ctx = mp.get_context("fork")
        with ctx.Pool(n) as pool:
            par = pool.map(_adapter_predict_guarded, [f"par_{i}" for i in range(n)])
        par_wall = time.perf_counter() - start

        seq_cpu = sum(r["elapsed_s"] for r in seq)
        par_cpu = sum(r["elapsed_s"] for r in par)
        log(
            f"adapter sequential: wall={seq_wall:.1f}s cpu={seq_cpu:.1f}s"
        )
        log(
            f"adapter parallel:   wall={par_wall:.1f}s cpu={par_cpu:.1f}s "
            f"speedup={seq_wall / par_wall if par_wall else 0:.2f}x vs sequential wall"
        )

        failed = False
        for kind, rows in ("seq", seq), ("par", par):
            for r in rows:
                if r["error"]:
                    log(f"ERROR: {kind} {r['job_name']}: {r['error']}")
                    failed = True
                    continue
                log(
                    f"  {kind} {r['job_name']}: n={r['n']} "
                    f"A={['%.2f' % a for a in r['A']]} {r['elapsed_s']:.1f}s"
                )
                if r["n"] != 2:
                    log(f"ERROR: {r['job_name']} expected 2 rows, got {r['n']}")
                    failed = True

        if failed:
            return 1

        baseline = seq[0]["A"]
        for r in seq[1:] + par:
            for a, b in zip(r["A"], baseline, strict=True):
                if abs(a - b) > DBA_TOL:
                    log(
                        f"ERROR: {r['job_name']} A={r['A']} differs from "
                        f"sequential baseline {baseline} by >{DBA_TOL} dB"
                    )
                    return 1
        log(f"adapter spectra match sequential within {DBA_TOL} dB")
        if par_wall >= seq_wall * 0.9:
            log(
                "NOTE: little/no wall-time speedup — wineserver may still serialize, "
                "but outputs did not clobber"
            )
        return 0


def main() -> int:
    if not AAM_SHIM.is_file():
        log(f"ERROR: shim not found at {AAM_SHIM}")
        return 1
    n = N_WORKERS
    if n < 2:
        log("ERROR: AAM_PARALLEL_N must be >= 2")
        return 1

    log(f"WINEPREFIX={os.environ.get('WINEPREFIX', '')} n={n}")
    log("skipping noisecon.inp (needs CH146.nc; staged NCfiles are FLATO200/OMNI_200)")
    return run_adapter_phase(n)


if __name__ == "__main__":
    raise SystemExit(main())
