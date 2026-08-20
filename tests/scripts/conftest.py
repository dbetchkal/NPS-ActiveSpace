"""Pytest fixtures for generate_active_space script tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from script_test_helpers import run_patched_generate_main


@pytest.fixture
def patched_generate_main(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    def _run(
        pool_output,
        select_optimal_return: tuple[str | None, float, float, pd.DataFrame] | None = None,
        extra_argv: list[str] | None = None,
    ) -> None:
        run_patched_generate_main(
            tmp_path,
            monkeypatch,
            pool_output=pool_output,
            select_optimal_return=select_optimal_return,
            extra_argv=extra_argv,
        )

    return _run
