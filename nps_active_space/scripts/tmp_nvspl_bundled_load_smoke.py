"""
Temporary smoke test: load bundled example NVSPL files and print a short report.

Run from the repository root (Mac/Linux/Windows):

    python -m nps_active_space.scripts.tmp_nvspl_bundled_load_smoke

Delete this module after NVSPL read optimization (#108) is validated on Windows hosts.
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import time
from pathlib import Path

import numpy as np

# Disable tqdm progress bars for cleaner timing output.
os.environ.setdefault("TQDM_DISABLE", "1")

BUNDLED_NVSPL_DIRS = (
    (
        "GLBALSTL 2024-05-24",
        Path(
            "example_data/nvspl_archive/2024 GLBALSTL Lower South Tidal Inlet"
            "/01 DATA/NVSPL"
        ),
    ),
    (
        "DENATRLA 2025-06-23",
        Path(
            "example_data/nvspl_archive/2025 DENATRLA Triple Lakes"
            "/01 DATA/NVSPL"
        ),
    ),
)

NUMERIC_SPOT_CHECKS = ("12.5", "1000", "20000", "dbA", "WindSpeed", "WindDir")
DROPPED_IF_EMPTY = (
    "INVID",
    "INSID",
    "GChar2",
    "AdjustmentsApplied",
    "CalibrationAdjustment",
    "GPSTimeAdjustment",
)


def find_repo_root(start: Path) -> Path | None:
    for candidate in [start, *start.parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "nps_active_space").is_dir():
            return candidate
    return None


def collect_nvspl_files(repo_root: Path) -> list[tuple[str, Path]]:
    entries: list[tuple[str, Path]] = []
    for label, rel_dir in BUNDLED_NVSPL_DIRS:
        nvspl_dir = repo_root / rel_dir
        if not nvspl_dir.is_dir():
            raise FileNotFoundError(f"Missing bundled NVSPL directory: {nvspl_dir}")
        files = sorted(nvspl_dir.glob("*.txt"))
        if not files:
            raise FileNotFoundError(f"No NVSPL .txt files in: {nvspl_dir}")
        entries.extend((label, path) for path in files)
    return entries


def format_mb(bytes_value: float) -> str:
    return f"{bytes_value / 1e6:.1f} MB"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Load bundled example NVSPL files and print a smoke-test report."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="Repository root (default: auto-detect from cwd or this file).",
    )
    args = parser.parse_args()

    repo_root = args.repo_root
    if repo_root is None:
        repo_root = find_repo_root(Path.cwd())
        if repo_root is None:
            repo_root = find_repo_root(Path(__file__).resolve().parents[2])
    if repo_root is None:
        print("ERROR: Could not find repository root. Run from repo root or pass --repo-root.")
        return 1

    repo_root = repo_root.resolve()
    os.chdir(repo_root)

    from nps_active_space.utils.models import Nvspl

    print("NVSPL bundled load smoke test")
    print(f"repo_root: {repo_root}")
    print(f"python: {sys.version.split()[0]}")
    print(f"platform: {sys.platform}")
    print()

    try:
        labeled_files = collect_nvspl_files(repo_root)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}")
        return 1

    file_paths = [str(path) for _, path in labeled_files]
    labels = {str(path): label for label, path in labeled_files}

    gc.collect()
    start = time.perf_counter()
    nvspl = Nvspl(file_paths)
    elapsed_s = time.perf_counter() - start
    mem_bytes = nvspl.memory_usage(deep=True).sum()

    print("Load summary")
    print(f"  files: {len(file_paths)}")
    print(f"  rows: {len(nvspl)}")
    print(f"  columns: {len(nvspl.columns)}")
    print(f"  time: {elapsed_s:.3f} s")
    print(f"  memory: {format_mb(mem_bytes)}")
    print(f"  index: {nvspl.index.min()} -> {nvspl.index.max()}")
    print(f"  monotonic index: {nvspl.index.is_monotonic_increasing}")
    print()

    print("Per-bundle file counts")
    for label in {label for label, _ in labeled_files}:
        count = sum(1 for path in file_paths if labels[path] == label)
        print(f"  {label}: {count} hourly files")
    print()

    print("Numeric dtypes (spot check)")
    for col in NUMERIC_SPOT_CHECKS:
        if col not in nvspl.columns:
            print(f"  {col}: MISSING")
        else:
            print(f"  {col}: {nvspl[col].dtype}")
    print()

    expected_numeric = set(Nvspl._NUMERIC_COLS) & set(nvspl.columns)
    non_float32 = sorted(
        col for col in expected_numeric if nvspl[col].dtype != np.float32
    )
    if non_float32:
        print(f"WARNING: expected float32 but got other dtypes: {non_float32}")
    else:
        print("All present _NUMERIC_COLS are float32")
    print()

    dropped = [col for col in DROPPED_IF_EMPTY if col not in nvspl.columns]
    present_empty_candidates = [col for col in DROPPED_IF_EMPTY if col in nvspl.columns]
    print(f"Dropped empty metadata columns: {dropped or 'none'}")
    if present_empty_candidates:
        print(f"Still present (had values): {present_empty_candidates}")
    print()

    print("dbA sanity")
    print(f"  min: {float(nvspl['dbA'].min()):.2f}")
    print(f"  max: {float(nvspl['dbA'].max()):.2f}")
    print(f"  mean: {float(nvspl['dbA'].mean()):.2f}")
    print()

    print("OK: bundled NVSPL load completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
