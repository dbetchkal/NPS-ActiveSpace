"""Tests for generate_active_space.py fail-guard exit paths without running NMSIM."""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest

import nps_active_space.scripts.generate_active_space as generate_active_space
from script_test_helpers import empty_active_space, nonempty_active_space
from nps_active_space.scripts.generate_active_space_batch import (
    REQUIRED_RESULT_KEYS,
    read_results_file,
)


class TestNonemptyActiveSpaceCount:
    def test_counts_only_layers_with_geometry(self):
        results = [
            ("O_+010", nonempty_active_space()),
            ("O_+020", gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")),
            ("O_+030", empty_active_space()),
        ]

        assert generate_active_space._nonempty_active_space_count(results) == 1


class TestGenerateActiveSpaceFailGuards:
    @pytest.mark.parametrize(
        ("pool_output", "select_optimal_return"),
        [
            (None, None),
            (("O_+010", empty_active_space(), {}), None),
            (
                ("O_+010", nonempty_active_space(), {}),
                ("O_+010", 0.0, 0.0, 0.0, pd.DataFrame()),
            ),
        ],
    )
    def test_exits_on_generation_failure(
        self,
        patched_generate_main,
        pool_output,
        select_optimal_return,
    ):
        with pytest.raises(SystemExit) as exc_info:
            patched_generate_main(
                pool_output=pool_output,
                select_optimal_return=select_optimal_return,
            )

        assert exc_info.value.code == 1


class TestGenerateActiveSpaceResultsOut:
    def test_success_writes_json_consumed_by_read_results_file(
        self,
        tmp_path: Path,
        patched_generate_main,
    ):
        results_path = tmp_path / "results.json"
        pool_output = ("O_+010", nonempty_active_space(), {})
        select_optimal_return = ("O_+010", 0.87, 0.0, 0.0, pd.DataFrame())

        patched_generate_main(
            pool_output=pool_output,
            select_optimal_return=select_optimal_return,
            extra_argv=["--results-out", str(results_path)],
        )

        data = json.loads(results_path.read_text())
        assert all(key in data for key in REQUIRED_RESULT_KEYS)
        assert data["F1"] == 0.87
        assert data["Mean altitude"] == 3000

        series, error_message = read_results_file(results_path, "DENATRLA20253000m")
        assert error_message is None
        assert series is not None
        assert series["Designator"] == "DENATRLA20253000m"
        assert series["F1"] == 0.87
