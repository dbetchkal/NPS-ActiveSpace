"""Tests for AAM terrain below-surface filtering."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from pyproj import Transformer
from shapely.geometry import Point

pytest.importorskip("aam_translator")

from aam_translator.constants import FT_PER_M

from nps_active_space.propagation_model.aam.terrain import (
    AAM_BELOW_SURFACE_TOLERANCE_M,
    _bilinear_sample_grid,
    _elv_grid_values,
    _hop_segment_below_terrain,
    _northup_row_from_model_j,
    _split_sequential_hop_runs,
    _terrain_surface_elevation_m,
    split_below_aam_terrain,
    split_safe_aam_track_runs,
    terrain_dir_for_site,
)

_TERRAIN_SAMPLING = "nps_active_space.propagation_model.aam.terrain_sampling"


def _patch_northup_gradient_grid(
    monkeypatch: pytest.MonkeyPatch,
    ridge_terrain,
) -> tuple[np.ndarray, float, int, int]:
    values = _elv_grid_values(Path(ridge_terrain.elv_path))
    scale = FT_PER_M if ridge_terrain.elv_header_feet else 1.0
    nrows, ncols = values.shape
    fake = np.zeros_like(values, dtype=np.float64)
    for i in range(nrows):
        fake[i, :] = 1000.0 * (1.0 - i / (nrows - 1)) * scale
    monkeypatch.setattr(
        f"{_TERRAIN_SAMPLING}._elv_grid_values",
        lambda _path: fake,
    )
    return fake, scale, nrows, ncols


class TestSplitBelowAamTerrain:
    def test_keeps_points_above_elv_surface(
        self,
        ridge_terrain,
        center_utm: tuple[float, float, float],
    ) -> None:
        x_m, y_m, z_m = center_utm
        source_pts = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(x_m, y_m, z_m + 100.0)],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(ridge_terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 0

    def test_filters_points_below_elv_surface(
        self,
        ridge_terrain,
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
        above, below = split_below_aam_terrain(ridge_terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 1
        assert above["id"].iloc[0] == 0
        assert below["id"].iloc[0] == 1

    def test_just_above_surface_passes(
        self,
        ridge_terrain,
        center_utm: tuple[float, float, float],
    ) -> None:
        x_m, y_m, _ = center_utm
        probe = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(x_m, y_m, 0.0)],
            crs="EPSG:26906",
        )
        surface_m = float(_terrain_surface_elevation_m(probe, ridge_terrain)[0])
        z_above_surface = surface_m + AAM_BELOW_SURFACE_TOLERANCE_M + 0.05
        source_pts = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(x_m, y_m, z_above_surface)],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(ridge_terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 0


def _utm_probe_at_elv_ij(ridge_terrain, col: float, row_south: float) -> gpd.GeoDataFrame:
    spec = ridge_terrain.spec
    aeqd_x_m = spec.grid_origin_x_m + col * spec.cell_dx_m
    aeqd_y_m = spec.grid_origin_y_m + row_south * spec.cell_dy_m
    from_aeqd = Transformer.from_crs(ridge_terrain.aeqd_crs, "EPSG:26906", always_xy=True)
    x_m, y_m = from_aeqd.transform(aeqd_x_m, aeqd_y_m)
    return gpd.GeoDataFrame(geometry=[Point(x_m, y_m, 0.0)], crs="EPSG:26906")


class TestElvNorthUpIndexing:
    def test_model_j_zero_is_south_array_row(self) -> None:
        assert float(_northup_row_from_model_j(np.array([0.0]), 873)[0]) == 872.0
        assert float(_northup_row_from_model_j(np.array([872.0]), 873)[0]) == 0.0

    def test_south_probe_matches_south_elv_not_north_row(
        self,
        ridge_terrain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        fake, scale, nrows, ncols = _patch_northup_gradient_grid(monkeypatch, ridge_terrain)
        col = ncols / 2.0
        row_south = 2.0
        sampled_m = float(
            _terrain_surface_elevation_m(
                _utm_probe_at_elv_ij(ridge_terrain, col, row_south),
                ridge_terrain,
            )[0]
        )
        row_north = _northup_row_from_model_j(np.array([row_south]), nrows)
        expected_raw = float(_bilinear_sample_grid(fake, np.array([col]), row_north)[0])
        wrong_raw = float(_bilinear_sample_grid(fake, np.array([col]), np.array([row_south]))[0])
        expected_m = expected_raw / scale
        wrong_m = wrong_raw / scale
        assert sampled_m == pytest.approx(expected_m, abs=0.05)
        assert abs(sampled_m - wrong_m) > 100.0

    def test_filters_against_south_surface_not_flipped_north_row(
        self,
        ridge_terrain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        fake, scale, nrows, ncols = _patch_northup_gradient_grid(monkeypatch, ridge_terrain)
        col = ncols / 2.0
        row_south = 2.0
        probe = _utm_probe_at_elv_ij(ridge_terrain, col, row_south)
        south_m = float(_terrain_surface_elevation_m(probe, ridge_terrain)[0])
        north_row_m = float(
            _bilinear_sample_grid(fake, np.array([col]), np.array([row_south]))[0]
        ) / scale
        z_m = (south_m + north_row_m) / 2.0
        assert z_m > south_m
        assert z_m < north_row_m
        source_pts = gpd.GeoDataFrame(
            {"id": [0]},
            geometry=[Point(probe.geometry.iloc[0].x, probe.geometry.iloc[0].y, z_m)],
            crs="EPSG:26906",
        )
        above, below = split_below_aam_terrain(ridge_terrain, source_pts)
        assert len(above) == 1
        assert len(below) == 0

    def test_ridge_south_of_center_differs_from_unflipped_row(self, ridge_terrain) -> None:
        values = _elv_grid_values(Path(ridge_terrain.elv_path))
        spec = ridge_terrain.spec
        col = spec.cell_count_x / 2.0
        row_south = spec.cell_count_y / 2.0 - 20.0
        sampled_m = float(
            _terrain_surface_elevation_m(
                _utm_probe_at_elv_ij(ridge_terrain, col, row_south),
                ridge_terrain,
            )[0]
        )
        row_north = _northup_row_from_model_j(np.array([row_south]), values.shape[0])
        expected_raw = float(
            _bilinear_sample_grid(values, np.array([col]), row_north)[0]
        )
        wrong_raw = float(
            _bilinear_sample_grid(values, np.array([col]), np.array([row_south]))[0]
        )
        expected_m = expected_raw / FT_PER_M if ridge_terrain.elv_header_feet else expected_raw
        wrong_m = wrong_raw / FT_PER_M if ridge_terrain.elv_header_feet else wrong_raw
        assert sampled_m == pytest.approx(expected_m, abs=0.05)
        assert abs(sampled_m - wrong_m) > 10.0


class TestSplitSafeAamTrackRuns:
    def test_keeps_one_run_when_hops_are_clear(
        self,
        ridge_terrain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2]},
            geometry=[Point(0, 0, 100), Point(10, 0, 100), Point(20, 0, 100)],
            crs="EPSG:32606",
        )
        monkeypatch.setattr(
            f"{_TERRAIN_SAMPLING}._hop_segment_below_terrain",
            lambda *args, **kwargs: False,
        )
        runs = split_safe_aam_track_runs(ridge_terrain, pts)
        assert len(runs) == 1
        assert runs[0]["id"].tolist() == [0, 1, 2]

    def test_splits_when_a_hop_clips_terrain(
        self,
        ridge_terrain,
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
            f"{_TERRAIN_SAMPLING}._hop_segment_below_terrain",
            fake_hop,
        )
        runs = split_safe_aam_track_runs(ridge_terrain, pts)
        assert [run["id"].tolist() for run in runs] == [[0, 1], [2]]

    def test_reconnects_around_a_clipping_snake_gap(
        self,
        ridge_terrain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        pts = gpd.GeoDataFrame(
            {"id": [0, 1, 2, 3]},
            geometry=[
                Point(0, 0, 100),
                Point(1, 0, 100),
                Point(2, 0, 100),
                Point(1, 1, 100),
            ],
            crs="EPSG:32606",
        )

        def fake_hop(terrain_ctx, start, end, source_crs, to_aeqd, from_aeqd) -> bool:
            xs = sorted((float(start.x), float(end.x)))
            ys = sorted((float(start.y), float(end.y)))
            return xs[0] < 1.5 < xs[1] and max(ys) < 0.5

        monkeypatch.setattr(
            f"{_TERRAIN_SAMPLING}._hop_segment_below_terrain",
            fake_hop,
        )
        sequential = _split_sequential_hop_runs(ridge_terrain, pts)
        packed = [run["id"].tolist() for run in split_safe_aam_track_runs(ridge_terrain, pts)]
        assert sequential == [[0, 1], [2, 3]]
        assert packed == [[0, 1, 3, 2]]

    def test_hop_interior_below_surface_is_detected(
        self,
        ridge_terrain,
        center_utm: tuple[float, float, float],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        x_m, y_m, z_m = center_utm
        start = Point(x_m - 200.0, y_m, z_m + 50.0)
        end = Point(x_m + 200.0, y_m, z_m + 50.0)
        to_aeqd = Transformer.from_crs("EPSG:26906", ridge_terrain.aeqd_crs, always_xy=True)
        from_aeqd = Transformer.from_crs(ridge_terrain.aeqd_crs, "EPSG:26906", always_xy=True)

        def ridge_at_midpoint(samples, terr):
            surface = np.full(len(samples), z_m, dtype=float)
            surface[len(samples) // 2] = z_m + 200.0
            return surface

        monkeypatch.setattr(
            f"{_TERRAIN_SAMPLING}._terrain_surface_elevation_m",
            ridge_at_midpoint,
        )
        assert _hop_segment_below_terrain(
            ridge_terrain, start, end, "EPSG:26906", to_aeqd, from_aeqd,
        ) is True

        monkeypatch.setattr(
            f"{_TERRAIN_SAMPLING}._terrain_surface_elevation_m",
            lambda samples, terr: np.full(len(samples), z_m, dtype=float),
        )
        assert _hop_segment_below_terrain(
            ridge_terrain, start, end, "EPSG:26906", to_aeqd, from_aeqd,
        ) is False


class TestTerrainDirForSite:
    def test_creates_canonical_when_missing(self, tmp_path: Path) -> None:
        expected = tmp_path / "Input_Data" / "aam" / "terrain" / "mic1"
        assert terrain_dir_for_site(tmp_path, "_mic1") == expected
        assert expected.is_dir()

    def test_returns_existing_canonical(self, tmp_path: Path) -> None:
        expected = tmp_path / "Input_Data" / "aam" / "terrain" / "mic1"
        expected.mkdir(parents=True)
        (expected / "scenario.elv").write_text("")

        assert terrain_dir_for_site(tmp_path, "_mic1") == expected
