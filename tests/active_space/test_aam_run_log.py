"""Tests for AAM site run log."""

from __future__ import annotations

import multiprocessing as mp
import re
from pathlib import Path

from nps_active_space.active_space.aam_run_log import (
    aam_log,
    aam_run_log_path,
    append_aam_run_summary,
    configure_aam_run_log,
    log_run_batch,
    short_aam_work_dir_name,
    summarize_aam_cli_output,
    summarize_aam_error,
)
from nps_active_space.active_space.propagation_model import AAM_RUN_LOG_FILENAME

_LOG_LINE_PREFIX = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z ")


def _strip_log_timestamp(line: str) -> str:
    return _LOG_LINE_PREFIX.sub("", line, count=1)


class TestAamRunLog:
    def test_configure_creates_log_with_session_header(self, tmp_path: Path) -> None:
        log_path = configure_aam_run_log(tmp_path)
        assert log_path == aam_run_log_path(tmp_path)
        assert log_path.name == AAM_RUN_LOG_FILENAME
        text = log_path.read_text(encoding="utf-8")
        assert text.startswith("=== session ")
        assert log_path.parent.name == "aam"

    def test_aam_log_appends_tagged_line(self, tmp_path: Path) -> None:
        configure_aam_run_log(tmp_path)
        aam_log("terrain", "grid 10×10 cells")
        lines = aam_run_log_path(tmp_path).read_text(encoding="utf-8").splitlines()
        assert _strip_log_timestamp(lines[-1]) == "[aam-terrain] grid 10×10 cells"
        assert _LOG_LINE_PREFIX.match(lines[-1])

    def test_log_run_batch_records_artifact_paths(self, tmp_path: Path) -> None:
        configure_aam_run_log(tmp_path)
        runs_dir = tmp_path / "Output_Data" / "aam" / "runs" / "job1"
        runs_dir.mkdir(parents=True)
        inp_path = runs_dir / "scenario.inp"
        inp_path.write_text("COMPUTEPOI\n", encoding="ascii")
        log_path = runs_dir / "scenario.txt"
        log_path.write_text("ok\n", encoding="utf-8")

        log_run_batch(
            tmp_path,
            job_name="job1",
            n_track=12,
            source_id="FLATO200",
            heading_deg=0.0,
            speed_kn=100.0,
            elapsed_s=6.2,
            inp_path=inp_path,
            aam_log_path=log_path,
        )

        last = _strip_log_timestamp(
            aam_run_log_path(tmp_path).read_text(encoding="utf-8").splitlines()[-1]
        )
        assert "[aam-run]" in last
        assert "job1" in last
        assert "n=12" in last
        assert "inp=Output_Data/aam/runs/job1/scenario.inp" in last
        assert "aam_log=Output_Data/aam/runs/job1/scenario.txt" in last

    def test_log_run_batch_failure_is_one_line_reason(self, tmp_path: Path) -> None:
        configure_aam_run_log(tmp_path)
        inp_path = tmp_path / "scenario.inp"
        inp_path.write_text("COMPUTEPOI\n", encoding="ascii")
        traceback = (
            "AAM exited 152: forrtl: severe (408): fort: (11): "
            "Subscript #2 of the array FPA has value 0\n"
            "Image              PC                Routine            Line        Source\n"
            "AAM_3.0.0.exe      00000001403D194C  Unknown               Unknown  Unknown\n"
        )
        log_run_batch(
            tmp_path,
            job_name="job_bad",
            n_track=1,
            source_id="FLATO200",
            heading_deg=0.0,
            speed_kn=0.0,
            elapsed_s=6.0,
            inp_path=inp_path,
            ok=False,
            error=traceback,
            to_console=False,
        )
        last = _strip_log_timestamp(
            aam_run_log_path(tmp_path).read_text(encoding="utf-8").splitlines()[-1]
        )
        assert "skip" in last
        assert "reason=AAM FPA bounds" in last
        assert "Unknown" not in last
        assert "FAILED" not in last

    def test_summarize_drops_fortran_stack(self) -> None:
        blob = (
            "forrtl: severe (408): fort: (11): Subscript #2 of the array FPA "
            "has value 0 which is less than the lower bound of 1\n\n"
            "Image              PC                Routine            Line        Source\n"
            "AAM_3.0.0.exe      00000001403D194C  Unknown               Unknown  Unknown\n"
            "kernel32.dll       000000007B627E49  Unknown               Unknown  Unknown\n"
        )
        assert "forrtl" in summarize_aam_cli_output(blob)
        assert "Unknown" not in summarize_aam_cli_output(blob)
        assert summarize_aam_error(blob) == "AAM FPA bounds"
        assert summarize_aam_error("empty .POI file: /tmp/scenario.POI") == "empty POI"

    def test_summarize_filename_path_length(self) -> None:
        blob = (
            "forrtl: severe (408): fort: (12): Variable FILENAME has substring "
            "ending point 164 which is greater than the variable length of 140"
        )
        assert summarize_aam_error(blob) == "AAM path too long (FILENAME 140)"

    def test_short_work_dir_name_is_stable_and_short(self) -> None:
        long_job = "DENASUSH2026O_+000_2100m_mesh1_r000_a0_a1_a2_a3"
        name = short_aam_work_dir_name(long_job)
        assert name == short_aam_work_dir_name(long_job)
        assert name.startswith("x")
        assert len(name) == 13

    def test_append_summary_block(self, tmp_path: Path) -> None:
        configure_aam_run_log(tmp_path)
        append_aam_run_summary(["gain=0.0 tested=100 audible=40"])
        text = aam_run_log_path(tmp_path).read_text(encoding="utf-8")
        assert "=== summary ===" in text
        assert "gain=0.0 tested=100 audible=40" in text
        summary_line = next(
            line for line in text.splitlines() if "gain=0.0 tested=100 audible=40" in line
        )
        assert _LOG_LINE_PREFIX.match(summary_line)


def _append_log_lines(root: str, tag: str, n: int) -> None:
    configure_aam_run_log(root)
    for i in range(n):
        aam_log("run", f"{tag}-{i:03d}")


class TestAamRunLogConcurrency:
    def test_concurrent_appends_keep_lines_intact(self, tmp_path: Path) -> None:
        configure_aam_run_log(tmp_path)
        n_per_worker = 40
        ctx = mp.get_context("spawn")
        with ctx.Pool(2) as pool:
            pool.starmap(
                _append_log_lines,
                [(str(tmp_path), "a", n_per_worker), (str(tmp_path), "b", n_per_worker)],
            )
        run_lines = [
            _strip_log_timestamp(line)
            for line in aam_run_log_path(tmp_path).read_text(encoding="utf-8").splitlines()
            if "[aam-run]" in line
        ]
        assert len(run_lines) == n_per_worker * 2
        assert {line.split("] ", 1)[1] for line in run_lines} == {
            f"{tag}-{i:03d}" for tag in ("a", "b") for i in range(n_per_worker)
        }
