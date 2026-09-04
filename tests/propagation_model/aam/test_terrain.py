"""Tests for AAM terrain below-surface filtering."""

from __future__ import annotations

import json
import os
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
from aam_translator.constants import FT_PER_M

from nps_active_space.propagation_model.aam.terrain import (
    AAM_BELOW_SURFACE_TOLERANCE_M,
    _bilinear_sample_grid,
    _elv_grid_values,
    _hop_segment_below_terrain,
    _northup_row_from_model_j,
    _terrain_surface_elevation_m,
    split_below_aam_terrain,
    split_safe_aam_track_runs,
    terrain_dir_for_site,
)

FIXTURES = Path(__file__).resolve().parents[2] / "active_space" / "fixtures" / "two_point_ridge"


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


def _utm_probe_at_elv_ij(terrain, col: float, row_south: float) -> gpd.GeoDataFrame:
    spec = terrain.spec
    aeqd_x_m = spec.grid_origin_x_m + col * spec.cell_dx_m
    aeqd_y_m = spec.grid_origin_y_m + row_south * spec.cell_dy_m
    from_aeqd = Transformer.from_crs(terrain.aeqd_crs, "EPSG:26906", always_xy=True)
    x_m, y_m = from_aeqd.transform(aeqd_x_m, aeqd_y_m)
    return gpd.GeoDataFrame(geometry=[Point(x_m, y_m, 0.0)], crs="EPSG:26906")


class TestElvNorthUpIndexing:
    def test_model_j_zero_is_south_array_row(self) -> None:
        assert float(_northup_row_from_model_j(np.array([0.0]), 873)[0]) == 872.0
        assert float(_northup_row_from_model_j(np.array([872.0]), 873)[0]) == 0.0

    def test_south_probe_matches_south_elv_not_north_row(
        self, terrain, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        values = _elv_grid_values(Path(terrain.elv_path))
        scale = FT_PER_M if terrain.elv_header_feet else 1.0
        nrows, ncols = values.shape
        # North-up gradient: row 0 (north) = 1000 m, last row (south) = 0 m.
        fake = np.zeros_like(values, dtype=np.float64)
        for i in range(nrows):
            fake[i, :] = 1000.0 * (1.0 - i / (nrows - 1)) * scale
        monkeypatch.setattr(
            "nps_active_space.propagation_model.aam.terrain._elv_grid_values",
            lambda _path: fake,
        )
        col = ncols / 2.0
        row_south = 2.0
        sampled_m = float(
            _terrain_surface_elevation_m(_utm_probe_at_elv_ij(terrain, col, row_south), terrain)[0]
        )
        row_north = _northup_row_from_model_j(np.array([row_south]), nrows)
        expected_raw = float(_bilinear_sample_grid(fake, np.array([col]), row_north)[0])
        wrong_raw = float(_bilinear_sample_grid(fake, np.array([col]), np.array([row_south]))[0])
        expected_m = expected_raw / scale
        wrong_m = wrong_raw / scale
        assert sampled_m == pytest.approx(expected_m, abs=0.05)
        assert abs(sampled_m - wrong_m) > 100.0

    def test_filters_against_south_surface_not_flipped_north_row(
        self, terrain, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        values = _elv_grid_values(Path(terrain.elv_path))
        scale = FT_PER_M if terrain.elv_header_feet else 1.0
        nrows, ncols = values.shape
        fake = np.zeros_like(values, dtype=np.float64)
        for i in range(nrows):
            fake[i, :] = 1000.0 * (1.0 - i / (nrows - 1)) * scale
        monkeypatch.setattr(
            "nps_active_space.propagation_model.aam.terrain._elv_grid_values",
            lambda _path: fake,
        )
        col = ncols / 2.0
        row_south = 2.0
        probe = _utm_probe_at_elv_ij(terrain, col, row_south)
        south_m = float(_terrain_surface_elevation_m(probe, terrain)[0])
        north_row_m = float(
            _bilinear_sample_grid(fake, np.array([col]), np.array([row_south]))[0]
        ) / scale
        # Midway: above the south cell, below the unflipped north-row reading.
        z_m = (south_m + north_row_m) / 2.0
        assert z_m > south_m
        assert z_m < north_row_m
        source_pts = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(probe.geometry.iloc[0].x, probe.geometry.iloc[0].y, z_m)],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 0

    def test_ridge_south_of_center_differs_from_unflipped_row(self, terrain) -> None:
        values = _elv_grid_values(Path(terrain.elv_path))
        spec = terrain.spec
        col = spec.cell_count_x / 2.0
        row_south = spec.cell_count_y / 2.0 - 20.0
        sampled_m = float(
            _terrain_surface_elevation_m(_utm_probe_at_elv_ij(terrain, col, row_south), terrain)[0]
        )
        row_north = _northup_row_from_model_j(np.array([row_south]), values.shape[0])
        expected_raw = float(
            _bilinear_sample_grid(values, np.array([col]), row_north)[0]
        )
        wrong_raw = float(
            _bilinear_sample_grid(values, np.array([col]), np.array([row_south]))[0]
        )
        expected_m = expected_raw / FT_PER_M if terrain.elv_header_feet else expected_raw
        wrong_m = wrong_raw / FT_PER_M if terrain.elv_header_feet else wrong_raw
        assert sampled_m == pytest.approx(expected_m, abs=0.05)
        assert abs(sampled_m - wrong_m) > 10.0


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
            "nps_active_space.propagation_model.aam.terrain._hop_segment_below_terrain",
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
            "nps_active_space.propagation_model.aam.terrain._hop_segment_below_terrain",
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
            "nps_active_space.propagation_model.aam.terrain._terrain_surface_elevation_m",
            ridge_at_midpoint,
        )
        assert _hop_segment_below_terrain(
            terrain, start, end, "EPSG:26906", to_aeqd, from_aeqd,
        ) is True

        monkeypatch.setattr(
            "nps_active_space.propagation_model.aam.terrain._terrain_surface_elevation_m",
            lambda samples, terr: np.full(len(samples), z_m, dtype=float),
        )
        assert _hop_segment_below_terrain(
            terrain, start, end, "EPSG:26906", to_aeqd, from_aeqd,
        ) is False


class TestTerrainDirForSite:
    def test_falls_back_to_legacy_flat_terrain_suffix(self, tmp_path: Path) -> None:
        legacy = tmp_path / "Input_Data" / "AAM" / "terrain_mic1"
        legacy.mkdir(parents=True)

        assert terrain_dir_for_site(tmp_path, "_mic1") == legacy

    def test_creates_canonical_when_missing(self, tmp_path: Path) -> None:
        expected = tmp_path / "Input_Data" / "aam" / "terrain" / "mic1"
        assert terrain_dir_for_site(tmp_path, "_mic1") == expected
        assert expected.is_dir()

    @pytest.mark.skipif(
        os.name == "nt" or os.uname().sysname == "Darwin",
        reason="case-insensitive filesystem collapses Input_Data/aam and Input_Data/AAM",
    )
    def test_prefers_canonical_lowercase_path(self, tmp_path: Path) -> None:
        canonical = tmp_path / "Input_Data" / "aam" / "terrain" / "mic1"
        canonical.mkdir(parents=True)
        legacy = tmp_path / "Input_Data" / "AAM" / "terrain" / "mic1"
        legacy.mkdir(parents=True)

        assert terrain_dir_for_site(tmp_path, "_mic1") == canonical

    @pytest.mark.skipif(
        os.name == "nt" or os.uname().sysname == "Darwin",
        reason="case-insensitive filesystem collapses Input_Data/aam and Input_Data/AAM",
    )
    def test_falls_back_to_legacy_uppercase_terrain_dir(self, tmp_path: Path) -> None:
        legacy = tmp_path / "Input_Data" / "AAM" / "terrain" / "mic1"
        legacy.mkdir(parents=True)

        assert terrain_dir_for_site(tmp_path, "_mic1") == legacy
