"""Docker integration: sea-level mesh clearance + one NMSim/AAM predict."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
VENDOR_ROOT = REPO.parent.parent if REPO.parent.name == ".worktrees" else REPO
VALIDATE = "docker/validate_sea_level_sources.py"
RUN = REPO / "docker" / "run_activespace.sh"


def _docker_available() -> bool:
    if shutil.which("docker") is None:
        return False
    probe = subprocess.run(
        ["docker", "info"],
        capture_output=True,
        timeout=30,
    )
    return probe.returncode == 0


def _vendor_runtime_ready(model: str, vendor_root: Path) -> bool:
    if model == "nmsim":
        runtime = vendor_root / "vendor" / "nmsim-runtime"
        return (
            (runtime / "Nord2000batch.exe").is_file()
            and (runtime / "RND" / "directories.ini").is_file()
        )
    runtime = vendor_root / "vendor" / "aam-runtime"
    aam_exe = os.environ.get("AAM_EXE", "AAM_3.0.0.exe")
    return (runtime / aam_exe).is_file() and (runtime / "NCfiles").is_dir()


@pytest.mark.integration
@pytest.mark.parametrize("model", ["nmsim", "aam"])
def test_sea_level_sources_in_docker(model: str) -> None:
    if not _docker_available():
        pytest.skip("docker not available")
    if not (REPO / VALIDATE).is_file() or not RUN.is_file():
        pytest.skip("docker validate script missing")
    if not _vendor_runtime_ready(model, VENDOR_ROOT):
        pytest.skip(
            f"vendor {model} runtime not staged under {VENDOR_ROOT / 'vendor'} "
            "(see docker/README.md)"
        )

    env = os.environ.copy()
    env.setdefault("NMSIM_RUNTIME", str(VENDOR_ROOT / "vendor" / "nmsim-runtime"))
    env.setdefault("AAM_RUNTIME", str(VENDOR_ROOT / "vendor" / "aam-runtime"))

    cmd = [str(RUN), "-m", model, str(VALIDATE), "--model", model]
    proc = subprocess.run(
        cmd,
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
    )
    combined = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        pytest.fail(
            f"docker sea-level {model} failed (exit {proc.returncode}):\n{combined[-8000:]}"
        )
    assert "[sea-level-sources] OK:" in combined
