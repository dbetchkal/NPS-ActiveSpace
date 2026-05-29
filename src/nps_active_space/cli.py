"""Console entry points for pipeline scripts (see pyproject.toml [project.scripts])."""

from __future__ import annotations

import runpy
import sys
import warnings
from typing import NoReturn

# Parsed here so script argparse never sees them. Default matches historic ``-W ignore``.
_SHOW_WARNINGS = "--show-warnings"
_IGNORE_WARNINGS = "--ignore-warnings"
def _apply_console_argv(argv: list[str] | None = None) -> bool:
    """Strip console-only flags from argv; return whether to show warnings."""
    target = sys.argv if argv is None else argv

    show_warnings = False
    script_argv = [target[0]]
    for arg in target[1:]:
        if arg == _SHOW_WARNINGS:
            show_warnings = True
        elif arg == _IGNORE_WARNINGS:
            show_warnings = False
        else:
            script_argv.append(arg)

    target[:] = script_argv
    return show_warnings


def _configure_runtime(*, show_warnings: bool) -> None:
    if not show_warnings:
        warnings.filterwarnings("ignore")
    for stream in (sys.stdout, sys.stderr):
        if stream is not None and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(line_buffering=True)
            except (OSError, ValueError):
                pass


def _run_script(module: str) -> NoReturn:
    show_warnings = _apply_console_argv()
    _configure_runtime(show_warnings=show_warnings)
    runpy.run_module(module, run_name="__main__")
    raise SystemExit(0)


def run_ground_truthing() -> NoReturn:
    _run_script("nps_active_space.scripts.run_ground_truthing")


def run_audible_transits() -> NoReturn:
    _run_script("nps_active_space.scripts.run_audible_transits")


def generate_active_space() -> NoReturn:
    _run_script("nps_active_space.scripts.generate_active_space")


def generate_active_space_batch() -> NoReturn:
    _run_script("nps_active_space.scripts.generate_active_space_batch")


def generate_active_space_mesh() -> NoReturn:
    _run_script("nps_active_space.scripts.generate_active_space_mesh")


def generate_3d_active_space() -> NoReturn:
    _run_script("nps_active_space.scripts.generate_3d_active_space")


def fit_3d_active_space() -> NoReturn:
    _run_script("nps_active_space.scripts.fit_3d_active_space")


def get_acoustic_metrics() -> NoReturn:
    _run_script("nps_active_space.scripts.get_acoustic_metrics")


def get_geographic_metrics() -> NoReturn:
    _run_script("nps_active_space.scripts.get_geographic_metrics")


def plot_altitudes() -> NoReturn:
    _run_script("nps_active_space.scripts.plot_altitudes")


def viz() -> NoReturn:
    _run_script("nps_active_space.scripts.viz")


def check_study_duration_robustness() -> NoReturn:
    _run_script("nps_active_space.scripts.check_study_duration_robustness")


def run_init_config() -> NoReturn:
    _run_script("nps_active_space.scripts.init_config")
