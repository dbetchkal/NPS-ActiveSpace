#!/usr/bin/env python
"""NMSim Wine smoke — Nord2000batch.exe launches through the container shim.

Does not compute an active space. Pipeline check:
  docker/run_activespace.sh nps_active_space/scripts/generate_active_space.py ...
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

SHIM = Path("/usr/local/bin/nord2000")
TIMEOUT_S = 60
WINE_FAIL_MARKERS = (
    "err:module:import_dll",
    "err:module:loader_init",
    "wine: could not load",
    "wine: cannot find",
)


def log(msg: str) -> None:
    print(f"[nmsim-smoke] {msg}", flush=True)


def main() -> int:
    home = Path(os.environ.get("NMSIM_HOME", "/opt/nmsim"))
    exe = home / "Nord2000batch.exe"
    rnd = home / "RND" / "directories.ini"

    if not SHIM.is_file():
        log(f"ERROR: shim not found at {SHIM} (rebuild image: docker/build.sh)")
        return 1
    if not exe.is_file():
        log(f"ERROR: {exe} missing")
        log("Stage: docker/stage_nmsim_runtime.sh /path/to/NMSim")
        return 1
    if not rnd.is_file():
        log(f"ERROR: {rnd} missing (RND/ is required at runtime)")
        return 1

    # Missing batch on purpose: we only need Wine to load the PE.
    proc = subprocess.run(
        [str(SHIM), "/tmp/nps_nmsim_smoke_missing.batch"],
        capture_output=True,
        text=True,
        timeout=TIMEOUT_S,
    )
    out = (proc.stdout or "") + (proc.stderr or "")
    if out:
        print(out, end="")

    lowered = out.lower()
    if any(marker in lowered for marker in WINE_FAIL_MARKERS):
        log("ERROR: Wine failed to launch Nord2000batch.exe")
        return 1

    log(f"OK: Nord2000batch.exe ran via Wine (exit {proc.returncode})")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.TimeoutExpired:
        log(f"ERROR: Nord2000batch.exe did not exit within {TIMEOUT_S}s")
        sys.exit(1)
