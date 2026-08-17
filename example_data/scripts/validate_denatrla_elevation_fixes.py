"""Validate elevation + set_dem fixes against DENATRLA example data.

Run from repo root:
    python example_data/scripts/validate_denatrla_elevation_fixes.py
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds

import nps_active_space.utils.config as cfg
from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
from nps_active_space.setup import setup_site
from nps_active_space.setup.elevation import get_project_setup_elevation
from nps_active_space.utils.enums import SourceElevationUnits
from nps_active_space.utils.helpers import get_deployment, load_studyarea

REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_DIR = REPO_ROOT / "example_data/site_projects"
SITE_DIR = PROJECT_DIR / "DENATRLA"
SOURCE_DEM = REPO_ROOT / "example_data/source_dems/DENATRLA_trla_m.tif"


def main() -> None:
    checks: list[tuple[str, bool]] = []

    tif, flt = get_project_setup_elevation(SITE_DIR)
    checks.append(
        ("get_project_setup_elevation", tif.name == "elevation_m_nad83_utm6.tif" and flt.is_file())
    )

    with rasterio.open(tif) as ds:
        ref = ds.read(1)
        nodata = np.int16(ds.nodata)
    with open(flt, "rb") as f:
        flt_data = np.frombuffer(f.read(), dtype=np.float32).reshape(ref.shape)
    valid = ref != nodata
    checks.append(
        (
            "tif_flt_values_match",
            bool(np.all(flt_data[valid] == ref[valid].astype(np.float32))),
        )
    )

    study_area = load_studyarea(str(PROJECT_DIR), "DENA", "TRLA", 2025)
    mic = get_deployment(str(PROJECT_DIR), "DENA", "TRLA", 2025, elevation=False)
    stub = REPO_ROOT / ".validation_nmsim_stub.exe"
    stub.write_text("stub")
    try:
        gen = ActiveSpaceGenerator(
            NMSIM=str(stub),
            study_area=study_area,
            root_dir=str(SITE_DIR),
            dem_src=str(SOURCE_DEM),
            ambience=pd.Series({"1000": 40.0}),
        )
        gen.set_dem(mic)
        checks.append(
            ("set_dem_paths", gen._dem_file.endswith(".tif") and gen._flt_file.endswith(".flt"))
        )
        checks.append(("set_dem_not_gdal_elevation.tif", "elevation.tif" not in Path(gen._dem_file).name))

        cfg.initialize(environment="DENA_example")
        cfg_dem = Path(cfg.read("data", "dem"))
        checks.append(("config_dem_exists", cfg_dem.is_file()))
        checks.append(("config_dem_is_meters_path", cfg_dem.name == "DENATRLA_trla_m.tif"))
        checks.append(("config_units_meters", cfg.read_dem_elevation_units() == SourceElevationUnits.METERS))

        tmpdir = Path(tempfile.mkdtemp(prefix="denatrla_validate_"))
        try:
            tmp_site = tmpdir / "DENATRLA"
            setup_site(
                site_dir=tmp_site,
                unit="DENA",
                site="TRLA",
                year=2025,
                mic_coord=(mic.lon, mic.lat),
                studyarea_sw=(-149.27687011178102, 63.49743904384526),
                studyarea_ne=(-148.4660630742463, 63.78496721672992),
                source_dem=cfg_dem,
                source_elevation_units=SourceElevationUnits.METERS,
            )
            out_tif = tmp_site / "Input_Data/01_ELEVATION/elevation_m_nad83_utm6.tif"
            out_flt = out_tif.with_suffix(".flt")
            checks.append(("project_setup_writes_flt", out_flt.is_file()))
            with rasterio.open(tif) as r, rasterio.open(out_tif) as o:
                win = from_bounds(*o.bounds, transform=r.transform)
                ref_clip = r.read(1)[
                    int(win.row_off):int(win.row_off + win.height),
                    int(win.col_off):int(win.col_off + win.width),
                ]
                out_data = o.read(1)
                rn, on = np.int16(r.nodata), np.int16(o.nodata)
                overlap_valid = (ref_clip != rn) & (out_data != on)
                checks.append(
                    (
                        "project_setup_roundtrip_exact",
                        bool(overlap_valid.any()) and bool(np.all(ref_clip[overlap_valid] == out_data[overlap_valid])),
                    )
                )
        finally:
            shutil.rmtree(tmpdir)
    finally:
        stub.unlink(missing_ok=True)

    print("DENATRLA validation results:")
    failed: list[str] = []
    for name, ok in checks:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}")
        if not ok:
            failed.append(name)

    if failed:
        raise SystemExit(f"{len(failed)} check(s) failed: {failed}")
    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
