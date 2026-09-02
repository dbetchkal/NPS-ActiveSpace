"""Tests for AAM terrain below-surface filtering."""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from pyproj import Transformer
from rasterio.transform import xy
from shapely.geometry import Point, box

pytest.importorskip("aam_translator")

from aam_translator import write_terrain

from nps_active_space.active_space.aam_terrain import (
    AAM_BELOW_SURFACE_TOLERANCE_M,
    _hop_segment_below_terrain,
    _terrain_surface_elevation_m,
    split_below_aam_terrain,
    split_safe_aam_track_runs,
)

FIXTURES = Path(__file__).parent / "fixtures" / "two_point_ridge"


def _fixture_dem_path() -> Path:
    meta = json.loads((FIXTURES / "case_meta.json").read_text())
    return FIXTURES / meta["dem_utm"]


def _dem_center_point_msl(dem_path: Path) -> tuple[float, float, float]:
    """Return UTM x, y and MSL elevation at the fixture DEM center."""
    with rasterio.open(dem_path) as dem:
        row_i, col_i = dem.height // 2, dem.width // 2
        x_m, y_m = xy(dem.transform, row_i, col_i)
        z_m = float(dem.read(1)[row_i, col_i])
    return x_m, y_m, z_m


def _aoi_for_dem(dem_path: Path) -> box:
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
def terrain(tmp_path: Path):
    dem_path = _fixture_dem_path()
    return write_terrain(
        dem_path,
        _aoi_for_dem(dem_path),
        tmp_path / "terrain",
        crs_in="EPSG:4326",
    )


@pytest.fixture
def center_utm() -> tuple[float, float, float]:
    return _dem_center_point_msl(_fixture_dem_path())


class TestSplitBelowAamTerrain:

    def test_keeps_points_above_elv_surface(
        self,
        terrain,
        center_utm: tuple[float, float, float],
    ) -> None:
        x_m, y_m, z_m = center_utm
        source_pts = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(x_m, y_m, z_m + 100.0)],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 0

    def test_filters_points_below_elv_surface(
        self,
        terrain,
        center_utm: tuple[float, float, float],
    ) -> None:
        x_m, y_m, z_m = center_utm
        source_pts = gpd.GeoDataFrame(
            {"id": [0, 1]},
            geometry=[
                Point(x_m, y_m, z_m + 50.0),
                Point(x_m, y_m, z_m - 50.0),
            ],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 1
        assert above["id"].iloc[0] == 0
        assert below["id"].iloc[0] == 1

    def test_just_above_surface_passes(
        self,
        terrain,
        center_utm: tuple[float, float, float],
    ) -> None:
        x_m, y_m, _ = center_utm
        probe = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(x_m, y_m, 0.0)],
            crs="EPSG:26906",
        )
        surface_m = float(_terrain_surface_elevation_m(probe, terrain)[0])
        z_above_surface = surface_m + AAM_BELOW_SURFACE_TOLERANCE_M + 0.05
        source_pts = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(x_m, y_m, z_above_surface)],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 0


class TestSplitSafeAamTrackRuns:
    def test_keeps_one_run_when_hops_are_clear(
        self,
        terrain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2]},
            geometry=[Point(0, 0, 100), Point(10, 0, 100), Point(20, 0, 100)],
            crs="EPSG:32606",
        )
        monkeypatch.setattr(
            "nps_active_space.active_space.aam_terrain._hop_segment_below_terrain",
            lambda *args, **kwargs: False,
        )
        runs = split_safe_aam_track_runs(terrain, pts)
        assert len(runs) == 1
        assert runs[0]["id"].tolist() == [0, 1, 2]

    def test_splits_when_a_hop_clips_terrain(
        self,
        terrain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2]},
            geometry=[Point(0, 0, 100), Point(10, 0, 100), Point(20, 0, 100)],
            crs="EPSG:32606",
        )

        def fake_hop(terrain_ctx, start, end, source_crs, to_aeqd, from_aeqd) -> bool:
            return float(start.x) == 10.0

        monkeypatch.setattr(
            "nps_active_space.active_space.aam_terrain._hop_segment_below_terrain",
            fake_hop,
        )
        runs = split_safe_aam_track_runs(terrain, pts)
        assert [run["id"].tolist() for run in runs] == [[0, 1], [2]]

    def test_hop_interior_below_surface_is_detected(
        self,
        terrain,
        center_utm: tuple[float, float, float],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        x_m, y_m, z_m = center_utm
        start = Point(x_m - 200.0, y_m, z_m + 50.0)
        end = Point(x_m + 200.0, y_m, z_m + 50.0)
        to_aeqd = Transformer.from_crs("EPSG:26906", terrain.aeqd_crs, always_xy=True)
        from_aeqd = Transformer.from_crs(terrain.aeqd_crs, "EPSG:26906", always_xy=True)

        def ridge_at_midpoint(samples, terr):
            surface = np.full(len(samples), z_m, dtype=float)
            surface[len(samples) // 2] = z_m + 200.0
            return surface

        monkeypatch.setattr(
            "nps_active_space.active_space.aam_terrain._terrain_surface_elevation_m",
            ridge_at_midpoint,
        )
        assert _hop_segment_below_terrain(
            terrain, start, end, "EPSG:26906", to_aeqd, from_aeqd,
        ) is True

        monkeypatch.setattr(
            "nps_active_space.active_space.aam_terrain._terrain_surface_elevation_m",
            lambda samples, terr: np.full(len(samples), z_m, dtype=float),
        )
        assert _hop_segment_below_terrain(
            terrain, start, end, "EPSG:26906", to_aeqd, from_aeqd,
        ) is False
