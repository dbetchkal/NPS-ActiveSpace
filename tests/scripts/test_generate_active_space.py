"""Tests for generate_active_space CLI."""

from __future__ import annotations

from nps_active_space.scripts.generate_active_space import (
    build_parser,
    resolve_pool_n_workers,
)
from nps_active_space.utils.enums import AcousticModel


class TestGenerateActiveSpaceParser:
    def test_model_flag(self):
        args = build_parser().parse_args([
            "-e", "container",
            "-u", "DENA",
            "-s", "TRLA",
            "-y", "2025",
            "--model", "aam",
        ])
        assert args.model is AcousticModel.AAM

    def test_model_defaults_none_for_env_resolution(self):
        args = build_parser().parse_args([
            "-e", "container",
            "-u", "DENA",
            "-s", "TRLA",
            "-y", "2025",
        ])
        assert args.model is None

    def test_density_optional(self):
        args = build_parser().parse_args([
            "-e", "container",
            "-u", "DENA",
            "-s", "TRLA",
            "-y", "2025",
            "--density", "12",
        ])
        assert args.density == 12

    def test_source_flag(self):
        args = build_parser().parse_args([
            "-e", "container",
            "-u", "DENA",
            "-s", "TRLA",
            "-y", "2025",
            "--source", "/path/O_+000.src",
            "--source", "/path/custom.nc",
        ])
        assert args.sources == ["/path/O_+000.src", "/path/custom.nc"]


class TestResolvePoolNWorkers:
    def test_aam_matches_nmsim_without_docker_cap(self, monkeypatch):
        monkeypatch.delenv("AAM_PARALLEL_N", raising=False)
        nmsim_n = resolve_pool_n_workers(AcousticModel.NMSIM)
        aam_n = resolve_pool_n_workers(AcousticModel.AAM)
        assert aam_n == nmsim_n

    def test_docker_script_can_cap_aam(self, monkeypatch):
        monkeypatch.setenv("AAM_PARALLEL_N", "2")
        assert resolve_pool_n_workers(AcousticModel.AAM) <= 2
