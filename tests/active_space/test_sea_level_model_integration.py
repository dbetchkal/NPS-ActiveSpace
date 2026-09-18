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
    return shutil.which("docker") is not None


@pytest.mark.integration
@pytest.mark.parametrize("model", ["nmsim", "aam"])
def test_sea_level_sources_in_docker(model: str) -> None:
    if not _docker_available():
        pytest.skip("docker not available")
    if not (REPO / VALIDATE).is_file() or not RUN.is_file():
        pytest.skip("docker validate script missing")

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
