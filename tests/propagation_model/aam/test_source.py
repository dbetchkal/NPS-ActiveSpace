"""Tests for AAM source resolution and NetCDF staging."""

from __future__ import annotations

import os
import shutil
import time
from pathlib import Path

import pytest

pytest.importorskip("aam_translator")

from aam_translator.bands import band_number_for_frequency

from nps_active_space.propagation_model.aam.source import (
    AAM_TEMPLATE_NC_FILENAME,
    aam_source_id_from_omni,
    aam_subprocess_env,
    ensure_aam_nc_for_source,
    omni_stem_to_aam_token,
    read_avg_spectrum_db,
    resolve_aam_template_ncfiles_dir,
    site_ncfiles_dir,
    stage_run_ncfiles,
    write_aam_nc,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_TEMPLATE = REPO_ROOT / "vendor/aam-runtime/NCfiles/OMNI_200.nc"
TUNING_O_000 = REPO_ROOT / "nps_active_space/propagation_model/nmsim/data/tuning/O_+000.avg"


@pytest.fixture
def template_nc() -> Path:
    if not VENDOR_TEMPLATE.is_file():
        pytest.skip(f"vendor template missing: {VENDOR_TEMPLATE}")
    return VENDOR_TEMPLATE


@pytest.fixture
def o_plus_000_src(tmp_path: Path) -> Path:
    if not TUNING_O_000.is_file():
        pytest.skip("tuning O_+000.avg missing")
    src_dir = tmp_path / "tuning"
    src_dir.mkdir()
    shutil.copy2(TUNING_O_000, src_dir / "O_+000.avg")
    src = src_dir / "O_+000.src"
    src.write_text("placeholder\n")
    return src


class TestOmniStemToAamToken:
    @pytest.mark.parametrize(
        ("stem", "expected"),
        [
            ("O_+000", "OMNI_000"),
            ("O_+200", "OMNI_200"),
            ("O_-100", "OMNIM100"),
        ],
    )
    def test_stem_maps_to_aam_token(self, stem: str, expected: str) -> None:
        assert omni_stem_to_aam_token(stem) == expected


class TestAamSourceIdFromOmni:
    @pytest.mark.parametrize(
        ("path", "expected"),
        [
            ("/data/tuning/O_+000.src", "OMNI_000"),
            ("/data/tuning/O_+200.avg", "OMNI_200"),
            ("/path/OMNI_000.nc", "OMNI_000"),
        ],
    )
    def test_path_maps_to_omni_token(self, path: str, expected: str) -> None:
        assert aam_source_id_from_omni(path) == expected


class TestReadAvgSpectrumDb:
    def test_o_plus_000_first_row(self) -> None:
        if not TUNING_O_000.is_file():
            pytest.skip("tuning O_+000.avg missing")
        levels = read_avg_spectrum_db(TUNING_O_000)
        band_25 = band_number_for_frequency(25.0)
        assert levels[band_25] == 56.0


class TestEnsureAamNcForSource:
    def test_src_avg_writes_cache(
        self,
        tmp_path: Path,
        template_nc: Path,
        o_plus_000_src: Path,
    ) -> None:
        token, cached = ensure_aam_nc_for_source(o_plus_000_src, tmp_path, template_nc)
        assert token == "OMNI_000"
        assert cached == site_ncfiles_dir(tmp_path) / "OMNI_000.nc"
        assert cached.is_file()

    def test_second_ensure_skips_when_mtime_unchanged(
        self,
        tmp_path: Path,
        template_nc: Path,
        o_plus_000_src: Path,
    ) -> None:
        ensure_aam_nc_for_source(o_plus_000_src, tmp_path, template_nc)
        cached = site_ncfiles_dir(tmp_path) / "OMNI_000.nc"
        first_mtime = cached.stat().st_mtime
        time.sleep(0.01)
        ensure_aam_nc_for_source(o_plus_000_src, tmp_path, template_nc)
        assert cached.stat().st_mtime == first_mtime

    def test_nc_pass_through_copies(self, tmp_path: Path, template_nc: Path) -> None:
        src_nc = tmp_path / "OMNI_042.nc"
        shutil.copy2(template_nc, src_nc)
        token, cached = ensure_aam_nc_for_source(src_nc, tmp_path, template_nc)
        assert token == "OMNI_042"
        assert cached.is_file()


class TestStageRunNcfiles:
    def test_stages_single_omni_only(self, tmp_path: Path, template_nc: Path) -> None:
        cache_dir = site_ncfiles_dir(tmp_path)
        cache_dir.mkdir(parents=True)
        omni_000 = cache_dir / "OMNI_000.nc"
        omni_005 = cache_dir / "OMNI_005.nc"
        shutil.copy2(template_nc, omni_000)
        shutil.copy2(template_nc, omni_005)

        work_dir = tmp_path / "runs" / "job_r000"
        staged = stage_run_ncfiles(work_dir, omni_005)

        assert staged == work_dir / "NCfiles"
        staged_names = sorted(p.name for p in staged.iterdir())
        assert staged_names == ["OMNI_005.nc"]


class TestAamSubprocessEnv:
    def test_sets_site_ncfiles_for_shim(self, tmp_path: Path) -> None:
        site_root = tmp_path / "site"
        nc_dir = site_ncfiles_dir(site_root)
        nc_dir.mkdir(parents=True)
        shim = tmp_path / "aam"
        shim.write_text("#!/bin/sh\n")
        env = aam_subprocess_env(shim, nc_dir)
        expected = str(nc_dir.resolve()) + os.sep
        assert env["ROTOR_NOISE"] == expected
        assert env["AAM_NC"] == str(nc_dir.resolve())

    def test_exe_sets_noise_paths_from_site_ncfiles(self, tmp_path: Path) -> None:
        site_root = tmp_path / "site"
        nc_dir = site_ncfiles_dir(site_root)
        nc_dir.mkdir(parents=True)
        exe = tmp_path / "AAM_3.0.0.exe"
        exe.write_bytes(b"")
        env = aam_subprocess_env(exe, nc_dir)
        expected = str(nc_dir.resolve()) + os.sep
        assert env["ROTOR_NOISE"] == expected
        assert env["FWING_NOISE"] == expected
        assert env["QUARRY_NOISE"] == expected
        assert env["AAM_NC"] == str(nc_dir.resolve())

    def test_template_resolution_prefers_parent_when_bin_stub_lacks_template(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.delenv("AAM_NC", raising=False)
        bin_dir = tmp_path / "Bin"
        bin_dir.mkdir()
        exe = bin_dir / "AAM_3.0.0.exe"
        exe.write_bytes(b"")
        stub = bin_dir / "NCfiles"
        stub.mkdir()
        nc = tmp_path / "NCfiles"
        nc.mkdir()
        (nc / AAM_TEMPLATE_NC_FILENAME).write_bytes(b"")
        resolved = resolve_aam_template_ncfiles_dir(exe)
        assert resolved == nc

    def test_aam_home_template_resolution(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("AAM_NC", raising=False)
        aam_home = tmp_path / "opt" / "aam"
        nc = aam_home / "NCfiles"
        nc.mkdir(parents=True)
        (nc / AAM_TEMPLATE_NC_FILENAME).write_bytes(b"")
        shim = tmp_path / "usr" / "local" / "bin" / "aam"
        shim.parent.mkdir(parents=True)
        shim.write_text("#!/bin/sh\n")
        monkeypatch.setenv("AAM_HOME", str(aam_home))
        resolved = resolve_aam_template_ncfiles_dir(shim)
        assert resolved == nc

    def test_aam_nc_override_for_template(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        exe = tmp_path / "Bin" / "AAM_3.0.0.exe"
        exe.parent.mkdir(parents=True)
        exe.write_bytes(b"")
        nc = tmp_path / "custom" / "NCfiles"
        nc.mkdir(parents=True)
        (nc / AAM_TEMPLATE_NC_FILENAME).write_bytes(b"")
        monkeypatch.setenv("AAM_NC", str(nc))
        resolved = resolve_aam_template_ncfiles_dir(exe)
        assert resolved == nc


class TestWriteAamNc:
    def test_writes_radius(self, tmp_path: Path, template_nc: Path) -> None:
        pytest.importorskip("netCDF4")
        from netCDF4 import Dataset

        levels = {14: 56.0}
        out = tmp_path / "OMNI_000.nc"
        write_aam_nc(levels, template_nc, out, radius_ft=1000.0)
        with Dataset(out) as ds:
            assert float(ds.variables["RADIUS"][:]) == 1000.0
