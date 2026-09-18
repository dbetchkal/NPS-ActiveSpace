"""Shared fixtures for AAM propagation model tests."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import pytest
import rasterio
from pyproj import Transformer
from rasterio.transform import xy
from shapely.geometry import Point, box

pytest.importorskip("aam_translator")

from aam_translator import write_terrain

from nps_active_space.propagation_model.aam.model import AamPropagationModel

TWO_POINT_RIDGE_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "two_point_ridge"
RIDGE_CRS = "EPSG:32606"


def fixture_dem_path() -> Path:
    meta = json.loads((TWO_POINT_RIDGE_FIXTURES / "case_meta.json").read_text())
    return TWO_POINT_RIDGE_FIXTURES / meta["dem_utm"]


def dem_center_point_msl(dem_path: Path) -> tuple[float, float, float]:
    """Return UTM x, y and MSL elevation at the fixture DEM center."""
    with rasterio.open(dem_path) as dem:
        row_i, col_i = dem.height // 2, dem.width // 2
        x_m, y_m = xy(dem.transform, row_i, col_i)
        z_m = float(dem.read(1)[row_i, col_i])
    return x_m, y_m, z_m


def aoi_for_dem(dem_path: Path) -> box:
    with rasterio.open(dem_path) as dem:
        to_wgs84 = Transformer.from_crs(dem.crs, "EPSG:4326", always_xy=True)
        lons, lats = to_wgs84.transform(
            [float(dem.bounds.left), float(dem.bounds.right)],
            [float(dem.bounds.bottom), float(dem.bounds.top)],
        )
        lon_min, lon_max = lons
        lat_min, lat_max = lats
    inset_deg = 0.0005
    return box(
        lon_min + inset_deg,
        lat_min + inset_deg,
        lon_max - inset_deg,
        lat_max - inset_deg,
    )


@pytest.fixture
def case_meta() -> dict:
    return json.loads((TWO_POINT_RIDGE_FIXTURES / "case_meta.json").read_text())


@pytest.fixture
def ridge_terrain(tmp_path: Path):
    dem_path = fixture_dem_path()
    return write_terrain(
        dem_path,
        aoi_for_dem(dem_path),
        tmp_path / "terrain",
        crs_in="EPSG:4326",
    )


@pytest.fixture
def center_utm() -> tuple[float, float, float]:
    return dem_center_point_msl(fixture_dem_path())


@pytest.fixture
def ridge_source_pts(case_meta: dict) -> gpd.GeoDataFrame:
    rows = case_meta["source_points_utm"]
    geoms = [Point(r["x"], r["y"], r["z"]) for r in rows]
    return gpd.GeoDataFrame(
        {"label": [r["label"] for r in rows]},
        geometry=geoms,
        crs=RIDGE_CRS,
    )


@pytest.fixture
def aam_predict_harness(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """AAM model with terrain filter and track splitting disabled for predict() tests."""

    def passthrough(self, site, source_pts, job_name=""):
        return source_pts, source_pts.iloc[0:0]

    monkeypatch.setattr(AamPropagationModel, "filter_below_terrain", passthrough)
    monkeypatch.setattr(
        "nps_active_space.propagation_model.aam.model.split_safe_aam_track_runs",
        lambda terrain, pts, job_name="": [pts],
    )
    model = AamPropagationModel(str(tmp_path))
    site = SimpleNamespace(terrain=None)
    return model, site
