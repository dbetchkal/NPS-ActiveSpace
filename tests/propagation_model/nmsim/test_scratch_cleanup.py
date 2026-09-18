"""Tests for NMSim scratch-file cleanup after a successful predict."""

from __future__ import annotations

from pathlib import Path

import pytest

from nps_active_space.propagation_model.nmsim.model import NmsimPropagationModel
from nps_active_space.utils.paths import NMSIM_SCRATCH_SUBDIR


def _model(tmp_path: Path) -> NmsimPropagationModel:
    nmsim_exe = tmp_path / "Nord2000batch.exe"
    nmsim_exe.write_text("")
    return NmsimPropagationModel(str(nmsim_exe), str(tmp_path))


def _job_files(tmp_path: Path, job: str = "TRLA_1000m_mesh1") -> tuple[Path, Path]:
    trj_dir = tmp_path / "Input_Data" / "03_TRAJECTORY"
    trj_dir.mkdir(parents=True)
    tis_dir = tmp_path / NMSIM_SCRATCH_SUBDIR
    tis_dir.mkdir(parents=True)
    return trj_dir / f"{job}.trj", tis_dir / f"{job}.tis"


class TestNmsimPredictScratchCleanup:
    def test_control_and_batch_paths_use_trajectory_stem(self, tmp_path: Path) -> None:
        model = _model(tmp_path)
        trj, _ = _job_files(tmp_path)
        control, batch = model._control_and_batch_paths(str(trj))
        scratch = tmp_path / NMSIM_SCRATCH_SUBDIR
        assert Path(control) == scratch / "control_TRLA_1000m_mesh1.nms"
        assert Path(batch) == scratch / "batch_TRLA_1000m_mesh1.txt"

    def test_cleanup_removes_trj_tis_control_and_batch(self, tmp_path: Path) -> None:
        model = _model(tmp_path)
        trj, tis = _job_files(tmp_path)
        trj.write_text("trajectory")
        tis.write_text("tis")
        model._create_instruction_files(
            "elev.flt", "site.sit", str(trj), "O_+000.src",
        )
        control, batch = model._control_and_batch_paths(str(trj))
        leftover = tmp_path / "keep_me.txt"
        leftover.write_text("other")

        model._cleanup_predict_scratch(str(trj), str(tis))

        assert not trj.exists()
        assert not tis.exists()
        assert not Path(control).exists()
        assert not Path(batch).exists()
        assert leftover.exists()

    def test_failed_read_leaves_scratch_files(self, tmp_path: Path) -> None:
        model = _model(tmp_path)
        trj, tis = _job_files(tmp_path)
        trj.write_text("not a trajectory")
        model._create_instruction_files(
            "elev.flt", "site.sit", str(trj), "O_+000.src",
        )
        control, batch = model._control_and_batch_paths(str(trj))

        with pytest.raises(Exception):
            model._read_trj_tis(str(trj), str(tis))

        assert trj.exists()
        assert Path(control).exists()
        assert Path(batch).exists()
