"""Append-only site log for AAM runs (``Output_Data/aam/active_space.log``)."""

from __future__ import annotations

import hashlib
import os
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

from nps_active_space.active_space.propagation_model import (
    AAM_OUTPUT_SUBDIR,
    AAM_RUN_LOG_FILENAME,
)

try:
    import fcntl
except ImportError:
    fcntl = None

_configured_root: Path | None = None
_log_path: Path | None = None


def aam_run_log_path(root_dir: str | Path) -> Path:
    """Return the site-relative path to the AAM run log file."""
    return Path(root_dir) / AAM_OUTPUT_SUBDIR / AAM_RUN_LOG_FILENAME


def configure_aam_run_log(root_dir: str | Path) -> Path:
    """Open (append) the site log and write a session header."""
    global _configured_root, _log_path

    root = Path(root_dir).resolve()
    if _configured_root == root and _log_path is not None:
        return _log_path

    _configured_root = root
    _log_path = aam_run_log_path(root)
    _log_path.parent.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    _append_line(f"=== session {timestamp} ===")
    return _log_path


def site_relative_path(root_dir: Path, path: str | Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(root_dir.resolve()))
    except ValueError:
        return str(path)


def summarize_aam_cli_output(text: str) -> str:
    """Keep the diagnostic line; drop Fortran/Wine stack frames."""
    if not text:
        return ""
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("Image") and "PC" in line:
            continue
        if line.endswith(".exe") or ".exe " in line[:40]:
            continue
        if line.startswith(("kernel32.dll", "ntdll.dll", "wine:")):
            continue
        if line.count("Unknown") >= 2:
            continue
        return line[:300]
    return text.splitlines()[0].strip()[:300]


def summarize_aam_error(message: str) -> str:
    """Short console/site-log reason for an AAM batch failure."""
    text = summarize_aam_cli_output(str(message))
    lower = text.lower()
    if "filename" in lower and ("length" in lower or "substring" in lower):
        return "AAM path too long (FILENAME 140)"
    if "fpa" in lower:
        return "AAM FPA bounds"
    if "forrtl" in lower:
        return text[:120]
    if "empty .poi" in lower or "no data rows" in lower:
        return "empty POI"
    if "read error" in lower:
        return "AAM READ ERROR"
    if "below aam terrain" in lower:
        return "below terrain"
    return text[:160]


def short_aam_work_dir_name(job_name: str) -> str:
    """Short ``runs/`` subdirectory name so ``ROTOR_NOISE`` paths stay under AAM's 140-char FILENAME cap."""
    digest = hashlib.sha1(job_name.encode()).hexdigest()[:12]
    return f"x{digest}"


def aam_log(category: str, msg: str, *, to_console: bool = True) -> None:
    """Write one line to the site log, and optionally stdout."""
    line = f"[aam-{category}] {msg}"
    if to_console:
        print(line, flush=True)
    _append_line(line)


def log_run_batch(
    root_dir: Path,
    *,
    job_name: str,
    n_track: int,
    source_id: str,
    heading_deg: float,
    speed_kn: float,
    elapsed_s: float,
    inp_path: Path,
    aam_log_path: Path | None = None,
    ok: bool = True,
    error: str | None = None,
    to_console: bool | None = None,
) -> None:
    """Log one AAM propagation batch with pointers to on-disk artifacts.

    Successful runs default to site log only; failures print to console unless
    ``to_console`` is set explicitly. Set ``AAM_VERBOSE=1`` for per-run console lines.
    """
    if to_console is None:
        verbose = os.environ.get("AAM_VERBOSE", "").strip().lower() in {
            "1",
            "true",
            "yes",
        }
        to_console = (not ok) or verbose
    rel_inp = site_relative_path(root_dir, inp_path)
    parts = [
        job_name,
        f"n={n_track}",
        f"source={source_id}",
        f"heading={heading_deg:g}",
        f"speed_kn={speed_kn:g}",
        f"{elapsed_s:.1f}s",
        f"inp={rel_inp}",
    ]
    if ok and aam_log_path is not None:
        parts.append(f"aam_log={site_relative_path(root_dir, aam_log_path)}")
    if not ok:
        parts.insert(0, "skip")
        if error:
            parts.append(f"reason={summarize_aam_error(error)}")
    aam_log("run", "  ".join(parts), to_console=to_console)


def append_aam_run_summary(lines: Iterable[str]) -> None:
    """Append a summary block (e.g. per-gain results from validate/generate scripts)."""
    _append_line("")
    _append_line("=== summary ===")
    for line in lines:
        _append_line(line)


def _append_line(line: str) -> None:
    if _log_path is None:
        return
    try:
        with _log_path.open("a", encoding="utf-8") as handle:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                handle.write(line + "\n")
            finally:
                if fcntl is not None:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    except OSError as exc:
        print(
            f"[aam-log] could not append to {_log_path} ({exc}); line dropped",
            flush=True,
        )
