"""Tests for AAM source resolution and NetCDF staging."""

from __future__ import annotations

import shutil
import time
from pathlib import Path

import pytest

pytest.importorskip("aam_translator")
pytest.importorskip("netCDF4")

from aam_translator.bands import band_number_for_frequency

from nps_active_space.active_space.aam_source import (
    aam_source_id_from_omni,
    ensure_aam_nc_for_source,
    omni_stem_to_aam_token,
    read_avg_spectrum_db,
    site_ncfiles_dir,
    write_aam_nc,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
VENDOR_TEMPLATE = REPO_ROOT / "vendor/aam-runtime/NCfiles/OMNI_200.nc"
TUNING_O_000 = REPO_ROOT / "nps_active_space/data/tuning/O_+000.avg"


@pytest.fixture
def template_nc() -> Path:
    if not VENDOR_TEMPLATE.is_file():
        pytest.skip(f"vendor template missing: {VENDOR_TEMPLATE}")
    return VENDOR_TEMPLATE


class TestOmniStemToAamToken:
    def test_o_plus_000(self) -> None:
        assert omni_stem_to_aam_token("O_+000") == "OMNI_000"

    def test_o_plus_200(self) -> None:
        assert omni_stem_to_aam_token("O_+200") == "OMNI_200"

    def test_o_minus_100(self) -> None:
        assert omni_stem_to_aam_token("O_-100") == "OMNIM100"


class TestAamSourceIdFromOmni:
    def test_src_path_maps_to_omni_token(self) -> None:
        assert aam_source_id_from_omni("/data/tuning/O_+000.src") == "OMNI_000"
        assert aam_source_id_from_omni("/data/tuning/O_+200.avg") == "OMNI_200"

    def test_nc_path_uses_stem(self) -> None:
        assert aam_source_id_from_omni("/path/OMNI_000.nc") == "OMNI_000"


class TestReadAvgSpectrumDb:
    def test_o_plus_000_first_row(self) -> None:
        if not TUNING_O_000.is_file():
            pytest.skip("tuning O_+000.avg missing")
        levels = read_avg_spectrum_db(TUNING_O_000)
        band_25 = band_number_for_frequency(25.0)
        assert levels[band_25] == 56.0


class TestEnsureAamNcForSource:
    def test_src_avg_writes_cache(self, tmp_path: Path, template_nc: Path) -> None:
        if not TUNING_O_000.is_file():
            pytest.skip("tuning O_+000.avg missing")
        src_dir = tmp_path / "tuning"
        src_dir.mkdir()
        avg = src_dir / "O_+000.avg"
        shutil.copy2(TUNING_O_000, avg)
        src = src_dir / "O_+000.src"
        src.write_text("placeholder\n")

        token, cached = ensure_aam_nc_for_source(src, tmp_path, template_nc)
        assert token == "OMNI_000"
        assert cached == site_ncfiles_dir(tmp_path) / "OMNI_000.nc"
        assert cached.is_file()

    def test_second_ensure_skips_when_mtime_unchanged(
        self,
        tmp_path: Path,
        template_nc: Path,
    ) -> None:
        if not TUNING_O_000.is_file():
            pytest.skip("tuning O_+000.avg missing")
        src_dir = tmp_path / "tuning"
        src_dir.mkdir()
        avg = src_dir / "O_+000.avg"
        shutil.copy2(TUNING_O_000, avg)
        src = src_dir / "O_+000.src"
        src.write_text("placeholder\n")

        ensure_aam_nc_for_source(src, tmp_path, template_nc)
        cached = site_ncfiles_dir(tmp_path) / "OMNI_000.nc"
        first_mtime = cached.stat().st_mtime
        time.sleep(0.01)
        ensure_aam_nc_for_source(src, tmp_path, template_nc)
        assert cached.stat().st_mtime == first_mtime

    def test_nc_pass_through_copies(self, tmp_path: Path, template_nc: Path) -> None:
        src_nc = tmp_path / "OMNI_042.nc"
        shutil.copy2(template_nc, src_nc)
        token, cached = ensure_aam_nc_for_source(src_nc, tmp_path, template_nc)
        assert token == "OMNI_042"
        assert cached.is_file()


class TestWriteAamNc:
    def test_writes_radius(self, tmp_path: Path, template_nc: Path) -> None:
        from netCDF4 import Dataset

        levels = {14: 56.0}
        out = tmp_path / "OMNI_000.nc"
        write_aam_nc(levels, template_nc, out, radius_ft=1000.0)
        with Dataset(out) as ds:
            assert float(ds.variables["RADIUS"][:]) == 1000.0
