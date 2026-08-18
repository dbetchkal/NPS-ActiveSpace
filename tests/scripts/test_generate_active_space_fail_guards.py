"""Tests for generate_active_space.py fail-guard exit paths without running NMSIM."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString, Polygon

import nps_active_space.scripts.generate_active_space as generate_active_space

MODULE = "nps_active_space.scripts.generate_active_space"


def _empty_active_space() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(geometry=[Polygon()], crs="EPSG:4326")


def _nonempty_active_space() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs="EPSG:4326",
    )


def _study_area() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"Unit": ["DENA"], "Site": ["TRLA"], "Year": [2025]},
        geometry=[Polygon([(-136.2, 58.4), (-135.8, 58.4), (-135.8, 58.8), (-136.2, 58.8)])],
        crs="EPSG:4269",
    )


def _annotations() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"audible": [True]},
        geometry=[LineString([(0, 0, 1000), (1, 1, 1000)])],
        crs="EPSG:4326",
    )


class _FakeTqdm:
    def __init__(self, iterable=None, **kwargs):
        self._iterable = iterable if iterable is not None else []

    def __iter__(self):
        return iter(self._iterable)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def update(self, n=1):
        pass


def _pool_factory(output):
    class _FakeAsyncResult:
        def get(self):
            return output

    class _FakePool:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def apply_async(self, func, callback=None, error_callback=None, kwds=None):
            if callback is not None:
                callback(None)
            return _FakeAsyncResult()

    return _FakePool


def _patch_main_dependencies(
    mod,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    pool_output,
    select_optimal_return: tuple[str | None, float, float, float, pd.DataFrame] | None = None,
) -> None:
    ambience_pkl = tmp_path / "ambience.pkl"
    pd.Series({"1000": 40.0}).to_pickle(ambience_pkl)
    project_dir = tmp_path / "project"
    site_dir = project_dir / "DENATRLA"
    site_dir.mkdir(parents=True)
    (site_dir / "DENATRLA_study_area.shp").write_text("placeholder")
    nmsim_exe = tmp_path / "nmsim.exe"
    nmsim_exe.write_text("stub")
    dem_path = tmp_path / "dem.tif"
    dem_path.write_text("stub")
    omni_source = tmp_path / "O_+010.src"
    omni_source.write_text("stub")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_active_space.py",
            "-e",
            "test_env",
            "-u",
            "DENA",
            "-s",
            "TRLA",
            "-y",
            "2025",
            "-l",
            "3000",
            "--ambience",
            str(ambience_pkl),
        ],
    )
    monkeypatch.setattr(mod.cfg, "initialize", lambda environment: None)
    monkeypatch.setattr(
        mod.cfg,
        "read",
        lambda section, key: {
            ("project", "dir"): str(project_dir),
            ("project", "nmsim"): str(nmsim_exe),
            ("data", "dem"): str(dem_path),
        }[(section, key)],
    )
    monkeypatch.setattr(mod, "get_omni_sources", lambda **kwargs: [str(omni_source)])
    monkeypatch.setattr(mod, "get_logger", lambda name: MagicMock())
    monkeypatch.setattr(mod, "get_deployment", lambda *args, **kwargs: MagicMock(name="DENATRLA2025"))
    monkeypatch.setattr(mod.gpd, "read_file", lambda path: _study_area())
    monkeypatch.setattr(
        mod.glob,
        "glob",
        MagicMock(return_value=[str(site_dir / "DENATRLA_study_area.shp")]),
    )
    monkeypatch.setattr(mod, "load_annotations", lambda *args, **kwargs: _annotations())
    monkeypatch.setattr(
        mod,
        "normalize_point_density",
        lambda valid_points, study_area, random_seed: valid_points,
    )
    monkeypatch.setattr(mod, "ActiveSpaceGenerator", lambda **kwargs: MagicMock())
    monkeypatch.setattr(mod.mp, "cpu_count", lambda: 2)
    monkeypatch.setattr(mod.mp, "Pool", _pool_factory(pool_output))
    monkeypatch.setattr(mod, "tqdm", _FakeTqdm)
    if select_optimal_return is not None:
        monkeypatch.setattr(mod, "select_optimal", lambda **kwargs: select_optimal_return)


def _run_generate_active_space_main(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    pool_output,
    select_optimal_return: tuple[str | None, float, float, float, pd.DataFrame] | None = None,
) -> None:
    mod = importlib.import_module(MODULE)
    _patch_main_dependencies(
        mod,
        tmp_path,
        monkeypatch,
        pool_output=pool_output,
        select_optimal_return=select_optimal_return,
    )
    source = Path(mod.__file__).read_text()
    main_block = source[source.index("if __name__ == '__main__':"):]
    namespace = mod.__dict__.copy()
    namespace["__name__"] = "__main__"
    exec(compile(main_block, mod.__file__, "exec"), namespace)


class TestNonemptyActiveSpaceCount:
    def test_counts_only_layers_with_geometry(self):
        results = [
            ("O_+010", _nonempty_active_space()),
            ("O_+020", gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")),
            ("O_+030", _empty_active_space()),
        ]

        assert generate_active_space._nonempty_active_space_count(results) == 1


class TestFailActiveSpaceGeneration:
    def test_exits_with_status_one(self, monkeypatch):
        generate_active_space.__dict__["logger"] = MagicMock()

        def capture_exit(code: int = 0) -> None:
            raise SystemExit(code)

        monkeypatch.setattr(generate_active_space.sys, "exit", capture_exit)

        with pytest.raises(SystemExit) as exc_info:
            generate_active_space._fail_active_space_generation("test failure")

        assert exc_info.value.code == 1


class TestGenerateActiveSpaceFailGuards:
    def test_exits_when_no_results_generated(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        with pytest.raises(SystemExit) as exc_info:
            _run_generate_active_space_main(tmp_path, monkeypatch, pool_output=None)

        assert exc_info.value.code == 1

    def test_exits_when_all_geojson_layers_are_empty(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        pool_output = ("O_+010", _empty_active_space(), {})

        with pytest.raises(SystemExit) as exc_info:
            _run_generate_active_space_main(tmp_path, monkeypatch, pool_output=pool_output)

        assert exc_info.value.code == 1

    def test_exits_when_f1_is_zero_with_audible_annotations(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        pool_output = ("O_+010", _nonempty_active_space(), {})
        select_optimal_return = ("O_+010", 0.0, 0.0, 0.0, pd.DataFrame())

        with pytest.raises(SystemExit) as exc_info:
            _run_generate_active_space_main(
                tmp_path,
                monkeypatch,
                pool_output=pool_output,
                select_optimal_return=select_optimal_return,
            )

        assert exc_info.value.code == 1
