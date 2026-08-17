"""
Temporary smoke test: load bundled example NVSPL files and print a short report.

Works on main and perf/nvspl-read-v2 — checkout either branch and run the same command.

From repo root (Mac/Linux/Windows):

    python -m nps_active_space.scripts.tmp_nvspl_bundled_load_smoke

Optional: average multiple timed runs (default 1):

    python -m nps_active_space.scripts.tmp_nvspl_bundled_load_smoke --runs 3

Workflow to compare main vs branch (single checkout):

    git checkout main
    python nps_active_space/scripts/tmp_nvspl_bundled_load_smoke.py --runs 3

    git checkout perf/nvspl-read-v2
    python nps_active_space/scripts/tmp_nvspl_bundled_load_smoke.py --runs 3

Or copy this file once to the repo root and run after each checkout:

    python tmp_nvspl_bundled_load_smoke.py --runs 3

Each run prints the current git branch/commit in the header.

Delete after NVSPL read optimization (#108) is validated on Windows hosts.
"""

from __future__ import annotations

import argparse
import gc
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

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

# Local copy — do not import from Nvspl (main does not define _NUMERIC_COLS).
NUMERIC_COLS = frozenset({
    "12.5", "15.8", "20", "25", "31.5", "40", "50", "63", "80", "100",
    "125", "160", "200", "250", "315", "400", "500", "630", "800", "1000",
    "1250", "1600", "2000", "2500", "3150", "4000", "5000", "6300", "8000",
    "10000", "12500", "16000", "20000", "dbA", "dbC", "dbF",
    "Voltage", "WindSpeed", "WindDir", "TempIns", "TempOut", "Humidity",
})

NUMERIC_SPOT_CHECKS = ("12.5", "1000", "20000", "dbA", "WindSpeed", "WindDir")
DROPPED_IF_EMPTY = (
    "INVID",
    "INSID",
    "GChar2",
    "AdjustmentsApplied",
    "CalibrationAdjustment",
    "GPSTimeAdjustment",
)


@dataclass
class LoadReport:
    files: int
    rows: int
    columns: int
    time_s: float
    mem_mb: float
    index_min: str
    index_max: str
    monotonic_index: bool
    dba_dtype: str
    dba_min: float
    dba_max: float
    dba_mean: float
    dropped_metadata: list[str]
    present_metadata: list[str]
    non_float32_numeric: list[str]
    spot_dtypes: dict[str, str]


def find_repo_root(start: Path) -> Path | None:
    for candidate in [start, *start.parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "nps_active_space").is_dir():
            return candidate
    return None


def git_label(repo_root: Path) -> str:
    try:
        branch = subprocess.run(
            ["git", "branch", "--show-current"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        return f"{branch} ({commit})" if branch else f"detached@{commit}"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


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


def load_once(repo_root: Path, file_paths: list[str]) -> LoadReport:
    for key in list(sys.modules):
        if key.startswith("nps_active_space"):
            del sys.modules[key]
    repo_str = str(repo_root)
    if repo_str in sys.path:
        sys.path.remove(repo_str)
    sys.path.insert(0, repo_str)

    from nps_active_space.utils.models import Nvspl

    gc.collect()
    start = time.perf_counter()
    nvspl = Nvspl(file_paths)
    elapsed_s = time.perf_counter() - start
    mem_bytes = nvspl.memory_usage(deep=True).sum()

    expected_numeric = NUMERIC_COLS & set(nvspl.columns)
    non_float32 = sorted(
        col for col in expected_numeric if nvspl[col].dtype != np.float32
    )
    dropped = [col for col in DROPPED_IF_EMPTY if col not in nvspl.columns]
    present_metadata = [col for col in DROPPED_IF_EMPTY if col in nvspl.columns]

    return LoadReport(
        files=len(file_paths),
        rows=len(nvspl),
        columns=len(nvspl.columns),
        time_s=elapsed_s,
        mem_mb=mem_bytes / 1e6,
        index_min=str(nvspl.index.min()),
        index_max=str(nvspl.index.max()),
        monotonic_index=bool(nvspl.index.is_monotonic_increasing),
        dba_dtype=str(nvspl["dbA"].dtype),
        dba_min=float(nvspl["dbA"].min()),
        dba_max=float(nvspl["dbA"].max()),
        dba_mean=float(nvspl["dbA"].mean()),
        dropped_metadata=dropped,
        present_metadata=present_metadata,
        non_float32_numeric=non_float32,
        spot_dtypes={
            col: str(nvspl[col].dtype) if col in nvspl.columns else "MISSING"
            for col in NUMERIC_SPOT_CHECKS
        },
    )


def print_report(
    repo_root: Path,
    labeled_files: list[tuple[str, Path]],
    file_paths: list[str],
    report: LoadReport,
    runs: int,
    time_avg_s: float | None,
    mem_avg_mb: float | None,
) -> None:
    labels = {str(path): label for label, path in labeled_files}

    print("NVSPL bundled load smoke test")
    print(f"git:       {git_label(repo_root)}")
    print(f"repo_root: {repo_root}")
    print(f"python:    {sys.version.split()[0]}")
    print(f"platform:  {sys.platform}")
    print()

    print("Load summary")
    print(f"  files:   {report.files}")
    print(f"  rows:    {report.rows}")
    print(f"  columns: {report.columns}")
    if runs > 1 and time_avg_s is not None:
        print(f"  time:    {time_avg_s:.3f} s avg ({runs} runs; last {report.time_s:.3f} s)")
        print(f"  memory:  {format_mb(mem_avg_mb * 1e6)} avg (last {format_mb(report.mem_mb * 1e6)})")
    else:
        print(f"  time:    {report.time_s:.3f} s")
        print(f"  memory:  {format_mb(report.mem_mb * 1e6)}")
    print(f"  index:   {report.index_min} -> {report.index_max}")
    print(f"  monotonic index: {report.monotonic_index}")
    print()

    print("Per-bundle file counts")
    for label in sorted({label for label, _ in labeled_files}):
        count = sum(1 for path in file_paths if labels[path] == label)
        print(f"  {label}: {count} hourly files")
    print()

    print("Numeric dtypes (spot check)")
    for col in NUMERIC_SPOT_CHECKS:
        print(f"  {col}: {report.spot_dtypes[col]}")
    print()

    if report.non_float32_numeric:
        print(f"Numeric columns not float32: {report.non_float32_numeric}")
    else:
        print("All present NUMERIC_COLS are float32")
    print()

    print(f"Dropped empty metadata columns: {report.dropped_metadata or 'none'}")
    if report.present_metadata:
        print(f"Still present (had values): {report.present_metadata}")
    print()

    print("dbA sanity")
    print(f"  min:  {report.dba_min:.2f}")
    print(f"  max:  {report.dba_max:.2f}")
    print(f"  mean: {report.dba_mean:.2f}")
    print()

    print("OK: bundled NVSPL load completed successfully.")


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
    parser.add_argument(
        "--runs",
        type=int,
        default=1,
        metavar="N",
        help="Number of timed loads; report average time and memory (default: 1).",
    )
    args = parser.parse_args()

    if args.runs < 1:
        print("ERROR: --runs must be >= 1")
        return 1

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

    try:
        labeled_files = collect_nvspl_files(repo_root)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}")
        return 1

    file_paths = [str(path) for _, path in labeled_files]

    times: list[float] = []
    mems: list[float] = []
    report: LoadReport | None = None
    for _ in range(args.runs):
        report = load_once(repo_root, file_paths)
        times.append(report.time_s)
        mems.append(report.mem_mb)

    time_avg = sum(times) / len(times) if args.runs > 1 else None
    mem_avg = sum(mems) / len(mems) if args.runs > 1 else None

    print_report(
        repo_root,
        labeled_files,
        file_paths,
        report,
        args.runs,
        time_avg,
        mem_avg,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
