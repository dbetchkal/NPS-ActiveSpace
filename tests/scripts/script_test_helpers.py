"""Importable helpers for generate_active_space script tests."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString, Polygon

MODULE = "nps_active_space.scripts.generate_active_space"


def make_stub_generate_active_space_script(
    *,
    altitude_m: int | None = None,
    segments: int = 5,
    kde_reduction: str = "12.5%",
    octave_gain: float = 12.5,
    f1: float = 0.91,
) -> str:
    """Return stub script source that writes JSON results for ``run_deployment`` tests."""
    altitude_lines = ""
    mean_altitude = "3000"
    if altitude_m is not None:
        altitude_lines = 'parser.add_argument("-l", "--altitude", type=int, required=True)\n'
        mean_altitude = "args.altitude"

    return f"""
import argparse
import json

parser = argparse.ArgumentParser()
{altitude_lines}parser.add_argument("--results-out", required=True)
args, _ = parser.parse_known_args()

with open(args.results_out, "w") as results_file:
    json.dump({{
        "Number of valid annotated segments": {segments},
        "Mean altitude": {mean_altitude},
        "KDE reduction (%)": "{kde_reduction}",
        "1/3rd Octave Gain (F1)": {octave_gain},
        "F1": {f1},
    }}, results_file)
"""


def stub_generate_active_space_cmd(
    tmp_path: Path,
    *,
    altitude_m: int | None = None,
) -> list[str]:
    stub_path = tmp_path / "stub_generate_active_space.py"
    stub_path.write_text(make_stub_generate_active_space_script(altitude_m=altitude_m))
    return [sys.executable, str(stub_path)]


def exec_script_main(module) -> None:
    source = Path(module.__file__).read_text()
    main_marker = 'if __name__ == "__main__":'
    if main_marker not in source:
        main_marker = "if __name__ == '__main__':"
    main_block = source[source.index(main_marker):]
    namespace = module.__dict__.copy()
    namespace["__name__"] = "__main__"
    exec(compile(main_block, module.__file__, "exec"), namespace)


def empty_active_space() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(geometry=[Polygon()], crs="EPSG:4326")


def nonempty_active_space() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs="EPSG:4326",
    )


def study_area() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"Unit": ["DENA"], "Site": ["TRLA"], "Year": [2025]},
        geometry=[Polygon([(-136.2, 58.4), (-135.8, 58.4), (-135.8, 58.8), (-136.2, 58.8)])],
        crs="EPSG:4269",
    )


def annotations() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {"audible": [True]},
        geometry=[LineString([(0, 0, 1000), (1, 1, 1000)])],
        crs="EPSG:4326",
    )


class FakeTqdm:
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


def pool_factory(output):
    class FakeAsyncResult:
        def get(self):
            return output

    class FakePool:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def apply_async(self, func, callback=None, error_callback=None, kwds=None):
            if callback is not None:
                callback(None)
            return FakeAsyncResult()

    return FakePool


def patch_generate_active_space_main(
    mod,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    pool_output,
    select_optimal_return: tuple[str | None, float, float, float, pd.DataFrame] | None = None,
    extra_argv: list[str] | None = None,
) -> None:
    ambience_pkl = tmp_path / "ambience.pkl"
    pd.Series({"1000": 40.0}).to_pickle(ambience_pkl)
    project_dir = tmp_path / "project"
    site_dir = project_dir / "DENATRLA"
    site_dir.mkdir(parents=True)
    (site_dir / "DENATRLA_study_area.shp").write_text("placeholder")
    nmsim_exe = tmp_path / "nmsim.exe"
    nmsim_exe.write_text("stub")
    omni_source = tmp_path / "O_+010.src"
    omni_source.write_text("stub")

    argv = [
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
    ]
    if extra_argv is not None:
        argv.extend(extra_argv)

    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(mod.cfg, "initialize", lambda environment: None)
    monkeypatch.setattr(
        mod.cfg,
        "read",
        lambda section, key: {
            ("project", "dir"): str(project_dir),
            ("project", "nmsim"): str(nmsim_exe),
        }[(section, key)],
    )
    monkeypatch.setattr(mod, "get_omni_sources", lambda **kwargs: [str(omni_source)])
    monkeypatch.setattr(mod, "get_logger", lambda name: MagicMock())
    monkeypatch.setattr(mod, "get_deployment", lambda *args, **kwargs: MagicMock(name="DENATRLA2025"))
    monkeypatch.setattr(mod.gpd, "read_file", lambda path: study_area())
    monkeypatch.setattr(
        mod.glob,
        "glob",
        MagicMock(return_value=[str(site_dir / "DENATRLA_study_area.shp")]),
    )
    monkeypatch.setattr(mod, "load_annotations", lambda *args, **kwargs: annotations())
    monkeypatch.setattr(
        mod,
        "normalize_point_density",
        lambda valid_points, study_area, random_seed: valid_points,
    )
    monkeypatch.setattr(mod, "ActiveSpaceGenerator", lambda **kwargs: MagicMock())
    monkeypatch.setattr(mod.mp, "cpu_count", lambda: 2)
    monkeypatch.setattr(mod.mp, "Pool", pool_factory(pool_output))
    monkeypatch.setattr(mod, "tqdm", FakeTqdm)
    if select_optimal_return is not None:
        monkeypatch.setattr(mod, "select_optimal", lambda **kwargs: select_optimal_return)


def run_patched_generate_main(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    pool_output,
    select_optimal_return: tuple[str | None, float, float, float, pd.DataFrame] | None = None,
    extra_argv: list[str] | None = None,
) -> None:
    mod = importlib.import_module(MODULE)
    patch_generate_active_space_main(
        mod,
        tmp_path,
        monkeypatch,
        pool_output=pool_output,
        select_optimal_return=select_optimal_return,
        extra_argv=extra_argv,
    )
    exec_script_main(mod)
