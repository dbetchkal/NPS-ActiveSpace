#!/usr/bin/env python
"""Docker integration test: AAM adapter two-point reciprocal ridge (tier-4 parity).

Run:
  docker/stage_aam_runtime.sh ~/dev/nmsim-aam-experiments/activespace-experiments/runs/tier4_reciprocal_two_point
  docker/run_activespace.sh -m aam docker/validate_aam_propagation_model.py

Requires ``FLATO200.nc`` in staged ``vendor/aam-runtime/NCfiles/`` (tier-4 case includes it).
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import geopandas as gpd
from rasterio import open as rio_open
from shapely.geometry import box, Point

from nps_active_space.active_space.aam_propagation_model import AamPropagationModel
from nps_active_space.utils.models import Microphone

TIER4 = Path(
    os.environ.get(
        "TIER4_FIXTURE_DIR",
        "/repo/tests/active_space/fixtures/tier4_two_point",
    ),
)
AAM_SHIM = Path("/usr/local/bin/aam")


def log(msg: str) -> None:
    print(f"[aam-adapter] {msg}", flush=True)


def main() -> int:
    if not AAM_SHIM.is_file():
        log(f"ERROR: {AAM_SHIM} missing (rebuild image: docker/build.sh)")
        return 1

    meta = json.loads((TIER4 / "case_meta.json").read_text())
    dem_path = TIER4 / "parent_dem_utm.tif"
    if not dem_path.is_file():
        dem_path = Path(meta.get("dem_utm", ""))
    if not dem_path.is_file():
        log(f"ERROR: DEM not found at {dem_path}")
        return 1

    with tempfile.TemporaryDirectory(prefix="aam_adapter_") as tmp:
        root = Path(tmp) / "site"
        (root / "Input_Data").mkdir(parents=True)

        rx_lon, rx_lat = meta["receiver_lonlat"]
        mic = Microphone(name="Receiver", lat=rx_lat, lon=rx_lon, z=4.92)

        rows = meta["source_points_utm"]
        crs = "EPSG:32606"
        geoms = [Point(r["x"], r["y"], r["z"]) for r in rows]
        source_pts = gpd.GeoDataFrame(
            {"label": [r["label"] for r in rows]},
            geometry=geoms,
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
            "/repo/nps_active_space/data/tuning/O_+200.avg",
        )
        preds = model.predict(
            site, source_pts, omni, altitude_m=400, job_name="tier4", heading=90,
        )

        if len(preds) != 2:
            log(f"ERROR: expected 2 predictions, got {len(preds)}")
            return 1
        if preds["A"].isna().any():
            log("ERROR: missing A-weighted levels")
            return 1

        log(f"OK: west={preds['A'].iloc[0]:.1f} dBA, east={preds['A'].iloc[1]:.1f} dBA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
