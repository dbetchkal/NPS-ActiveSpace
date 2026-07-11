#!/usr/bin/env python
"""AAM Wine smoke test — proves the binary runs in nps-activespace:linux.

Uses noisecon.inp from the staged runtime (vendor/aam-runtime/, mounted at /opt/aam).
Run: docker/run_activespace.sh -m aam docker/validate_aam_smoke.py
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

AAM_SHIM = Path("/usr/local/bin/aam")
TIMEOUT_S = 300


def log(msg: str) -> None:
    print(f"[aam-smoke] {msg}", flush=True)


def main() -> int:
    aam_home = Path(os.environ.get("AAM_HOME", "/opt/aam"))
    aam_exe = os.environ.get("AAM_EXE", "AAM_3.0.0.exe")
    noisecon = aam_home / "noisecon.inp"

    if not (aam_home / aam_exe).is_file():
        log(f"ERROR: {aam_home / aam_exe} missing")
        log("Stage: docker/stage_aam_runtime.sh /path/to/AAM  (dir with noisecon.inp)")
        log("Run:   docker/run_activespace.sh -m aam docker/validate_aam_smoke.py")
        return 1
    if not (aam_home / "NCfiles").is_dir():
        log("ERROR: NCfiles/ missing under AAM_HOME")
        return 1
    if not AAM_SHIM.is_file():
        log(f"ERROR: shim not found at {AAM_SHIM} (rebuild image: docker/build.sh)")
        return 1
    if not noisecon.is_file():
        log(f"ERROR: {noisecon} missing — re-stage from a dir that includes noisecon.inp")
        return 1

    with tempfile.TemporaryDirectory(prefix="aam_smoke_") as tmp:
        work = Path(tmp)
        shutil.copy(noisecon, work / "noisecon.inp")

        log(f"running aam noisecon.inp in {work} ...")
        proc = subprocess.run(
            [str(AAM_SHIM), "noisecon.inp"],
            cwd=work,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
        )
        if proc.returncode != 0:
            if proc.stdout:
                print(proc.stdout, end="", file=sys.stderr)
            if proc.stderr:
                print(proc.stderr, end="", file=sys.stderr)
            log(f"ERROR: aam exited {proc.returncode}")
            return proc.returncode or 1

        grd = work / "noisecon.GRD"
        txt = work / "noisecon.txt"
        if not grd.is_file() or grd.stat().st_size < 100:
            log(f"ERROR: expected noisecon.GRD (got {grd.stat().st_size if grd.is_file() else 'missing'} bytes)")
            return 1

        body = txt.read_text(errors="replace")
        if "ADVANCED ACOUSTIC MODEL" not in body:
            log("ERROR: noisecon.txt missing expected AAM log header")
            return 1
        if "ch146" not in body.lower():
            log("ERROR: noisecon.txt does not reference ch146 NC files (NC env vars?)")
            return 1

        log(f"OK: {grd.stat().st_size} byte GRD, log {txt.stat().st_size} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
