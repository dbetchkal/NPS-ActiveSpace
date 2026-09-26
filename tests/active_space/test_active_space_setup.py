"""Tests for shared active-space generator setup."""

from __future__ import annotations

import os

import pandas as pd
import pytest

from nps_active_space.active_space.active_space_setup import (
    DEFAULT_AAM_SHIM,
    omni_stem_to_gain_db,
    precision_recall_plot_path,
    resolve_3d_fit_gain,
    resolve_aam_exe,
    resolve_acoustic_model,
    upsert_project_fit,
)
from nps_active_space.utils.enums import AcousticModel


class TestResolveAcousticModel:
    def test_defaults_to_nmsim(self):
        assert resolve_acoustic_model(None) is AcousticModel.NMSIM

    def test_cli_overrides_env(self, monkeypatch):
        monkeypatch.setenv("ACOUSTIC_MODEL", "aam")
        assert resolve_acoustic_model(AcousticModel.NMSIM) is AcousticModel.NMSIM

    def test_env_aam(self, monkeypatch):
        monkeypatch.delenv("ACOUSTIC_MODEL", raising=False)
        monkeypatch.setenv("ACOUSTIC_MODEL", "aam")
        assert resolve_acoustic_model(None) is AcousticModel.AAM

    def test_parses_string_cli(self):
        assert resolve_acoustic_model("aam") is AcousticModel.AAM

    def test_invalid_model_raises(self):
        with pytest.raises(ValueError, match="Unknown acoustic model"):
            resolve_acoustic_model("invalid")


class TestResolveAamExe:
    def test_explicit_override(self):
        assert resolve_aam_exe(aam_exe=r"C:\AAM\AAM_3.0.0.exe") == r"C:\AAM\AAM_3.0.0.exe"

    def test_reads_project_aam(self, monkeypatch):
        monkeypatch.setattr(
            "nps_active_space.active_space.active_space_setup.cfg.read",
            lambda section, option=None: r"C:\AAM\AAM_3.0.0.exe",
        )
        assert resolve_aam_exe() == r"C:\AAM\AAM_3.0.0.exe"

    def test_blank_config_uses_docker_shim(self, monkeypatch):
        monkeypatch.setattr(
            "nps_active_space.active_space.active_space_setup.cfg.read",
            lambda section, option=None: "",
        )
        assert resolve_aam_exe() == DEFAULT_AAM_SHIM


class TestOmniStemToGainDb:
    def test_positive_gain(self):
        assert omni_stem_to_gain_db("O_+125") == 12.5

    def test_negative_gain(self):
        assert omni_stem_to_gain_db("O_-050") == -5.0


class TestUpsertProjectFit:
    def test_upsert_by_model(self, tmp_path):
        project_dir = tmp_path / "site_projects"
        project_dir.mkdir()
        upsert_project_fit(
            str(project_dir), "DENATRLA2025", AcousticModel.NMSIM,
            {"1/3rd Octave Gain (F1)": 2.0, "F1": 0.5},
        )
        upsert_project_fit(
            str(project_dir), "DENATRLA2025", AcousticModel.AAM,
            {"1/3rd Octave Gain (F1)": 0.5, "F1": 0.6},
        )
        df = pd.read_csv(project_dir / "fits.csv")
        assert len(df) == 2
        assert set(df["Model"]) == {"nmsim", "aam"}


class TestResolve3dFitGain:
    def test_reads_model_specific_row(self, tmp_path):
        project_dir = tmp_path / "site_projects"
        project_dir.mkdir()
        upsert_project_fit(
            str(project_dir), "DENATRLA2025", AcousticModel.AAM,
            {"1/3rd Octave Gain (F1)": 0.5, "F1": 0.6},
        )
        assert resolve_3d_fit_gain(
            str(project_dir), "DENA", "TRLA", 2025, AcousticModel.AAM,
        ) == 0.5
        assert resolve_3d_fit_gain(
            str(project_dir), "DENA", "TRLA", 2025, AcousticModel.NMSIM,
        ) is None


class TestPrecisionRecallPlotPath:
    def test_model_scoped_output(self, tmp_path):
        site_dir = tmp_path / "DENATRLA"
        site_dir.mkdir()
        path = precision_recall_plot_path(
            str(site_dir), "DENA", "TRLA", 2025, 1000, 1.0, AcousticModel.AAM,
        )
        assert path.endswith(
            "Output_Data/aam/PRECISION_RECALL/PrecisionRecallPlot_DENATRLA2025_1000m_1p0.png"
        )
