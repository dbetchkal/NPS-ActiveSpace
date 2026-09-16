"""Tests for 3D active-space command-line builders."""

from __future__ import annotations

from nps_active_space.scripts.generate_3d_commands import (
    build_layer_command_parts,
    format_commands_file_line,
)
from nps_active_space.utils.enums import AcousticModel


class TestBuildLayerCommandParts:
    def test_includes_model_and_ambience(self):
        parts = build_layer_command_parts(
            "container",
            "DENA",
            "TRLA",
            2025,
            1000,
            "/tmp/ambience.pkl",
            model=AcousticModel.AAM,
            extra_args=["--omni-min", "0", "--omni-max", "2"],
        )
        assert parts == [
            "-e", "container",
            "-u", "DENA",
            "-s", "TRLA",
            "-y", 2025,
            "-l", 1000,
            "--model", "aam",
            "-a", "/tmp/ambience.pkl",
            "--omni-min", "0",
            "--omni-max", "2",
        ]

    def test_format_commands_file_line_quotes_paths(self):
        line = format_commands_file_line(
            "DENATRLA2025_1000m",
            ["-a", "/tmp/with spaces.pkl"],
        )
        assert line.startswith("DENATRLA2025_1000m\t")
        assert "'/tmp/with spaces.pkl'" in line
