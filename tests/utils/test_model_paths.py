"""Tests for model-scoped output paths and NMSim legacy resolvers."""

from __future__ import annotations

from pathlib import Path

import pytest

from nps_active_space.utils import legacy_nmsim_paths as legacy
from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel


@pytest.fixture
def site_tree(tmp_path: Path) -> tuple[Path, str, str, int]:
    project_dir = tmp_path / "projects"
    unit, site, year = "DENA", "TRLA", 2025
    site_dir = project_dir / f"{unit}{site}"
    site_dir.mkdir(parents=True)
    return project_dir, unit, site, year


class TestLegacyNmsimResolvers:
    @pytest.mark.parametrize(
        ("resolver", "new_subpath", "legacy_subpath"),
        [
            (legacy.resolve_nmsim_predictions_dir, ("nmsim", "predictions"), ("TIG_TIS",)),
            (legacy.resolve_nmsim_scratch_dir, ("nmsim", "scratch"), ("TIG_TIS",)),
        ],
    )
    def test_prefers_new_when_dir_exists(self, site_tree, resolver, new_subpath, legacy_subpath) -> None:
        project_dir, unit, site, _year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        new_dir = site_dir.joinpath("Output_Data", *new_subpath)
        new_dir.mkdir(parents=True)

        assert resolver(str(site_dir)) == str(new_dir)

    @pytest.mark.parametrize(
        ("resolver", "legacy_subpath"),
        [
            (legacy.resolve_nmsim_predictions_dir, ("TIG_TIS",)),
            (legacy.resolve_nmsim_scratch_dir, ("TIG_TIS",)),
        ],
    )
    def test_falls_back_to_legacy(self, site_tree, resolver, legacy_subpath) -> None:
        project_dir, unit, site, _year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        legacy_dir = site_dir.joinpath("Output_Data", *legacy_subpath)
        legacy_dir.mkdir(parents=True)

        assert resolver(str(site_dir)) == str(legacy_dir)

    def test_predictions_write_always_uses_new(self, site_tree) -> None:
        project_dir, unit, site, _year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        expected = site_dir / "Output_Data" / "nmsim" / "predictions"

        assert legacy.resolve_nmsim_predictions_dir(str(site_dir), for_write=True) == str(expected)

    def test_activespaces_prefers_new_with_layers(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        layer_dir = site_dir / "Output_Data" / "nmsim" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        layer_dir.mkdir(parents=True)

        expected = site_dir / "Output_Data" / "nmsim" / "ACTIVESPACES"
        assert legacy.resolve_nmsim_activespaces_dir(str(site_dir)) == str(expected)

    def test_activespaces_falls_back_to_legacy(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        layer_dir = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        layer_dir.mkdir(parents=True)

        expected = site_dir / "Output_Data" / "ACTIVESPACES"
        assert legacy.resolve_nmsim_activespaces_dir(str(site_dir)) == str(expected)

    def test_layer_dirs_merges_new_and_legacy_by_altitude(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        new_layer = site_dir / "Output_Data" / "nmsim" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        new_layer.mkdir(parents=True)
        legacy_layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_500m"
        legacy_layer.mkdir(parents=True)

        matches = legacy.resolve_activespace_layer_dirs(str(project_dir), unit, site, year)
        assert matches == [str(legacy_layer), str(new_layer)]

    def test_layer_dirs_same_altitude_prefers_new(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        new_layer = site_dir / "Output_Data" / "nmsim" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        new_layer.mkdir(parents=True)
        legacy_layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        legacy_layer.mkdir(parents=True)

        matches = legacy.resolve_activespace_layer_dirs(str(project_dir), unit, site, year)
        assert matches == [str(new_layer)]

    def test_find_layer_geojson_falls_back_to_legacy(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        usy = f"{unit}{site}{year}"
        new_layer = site_dir / "Output_Data" / "nmsim" / "ACTIVESPACES" / f"{usy}_1500m"
        new_layer.mkdir(parents=True)
        (new_layer / f"{usy}_O_+000.geojson").write_text("{}")
        legacy_layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{usy}_1500m"
        legacy_layer.mkdir(parents=True)
        legacy_file = legacy_layer / f"{usy}_O_+020.geojson"
        legacy_file.write_text("{}")

        found = legacy.find_layer_geojson(str(new_layer), 2.0)
        assert found == str(legacy_file)

    def test_layer_dirs_falls_back_to_legacy(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        legacy_layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_500m"
        legacy_layer.mkdir(parents=True)

        matches = legacy.resolve_activespace_layer_dirs(str(project_dir), unit, site, year)
        assert matches == [str(legacy_layer)]

    def test_layer_dirs_excludes_experiment_aam_suffix_dirs(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        nmsim_layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        nmsim_layer.mkdir(parents=True)
        experiment = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m_aam"
        experiment.mkdir(parents=True)

        matches = legacy.resolve_activespace_layer_dirs(str(project_dir), unit, site, year)
        assert matches == [str(nmsim_layer)]

    def test_is_standard_altitude_layer_dir(self) -> None:
        usy = "DENATRLA2025"
        assert legacy.is_standard_altitude_layer_dir(f"/tmp/{usy}_1000m", usy)
        assert not legacy.is_standard_altitude_layer_dir(f"/tmp/{usy}_1000m_aam", usy)

    def test_geojson_prefers_existing_new_file(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        layer = site_dir / "Output_Data" / "nmsim" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        layer.mkdir(parents=True)
        geojson = layer / f"{unit}{site}{year}_O_+005.geojson"
        geojson.write_text("{}")

        path = legacy.resolve_activespace_geojson(
            str(project_dir), unit, site, year, 1000, "+", "005",
        )
        assert path == str(geojson)

    def test_geojson_returns_legacy_when_new_missing(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        layer.mkdir(parents=True)
        expected = layer / f"{unit}{site}{year}_O_+005.geojson"

        path = legacy.resolve_activespace_geojson(
            str(project_dir), unit, site, year, 1000, "+", "005",
        )
        assert path == str(expected)


class TestDisplayPath:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            (r"C:\projects\DENATRLA\Output_Data\nmsim", "C:/projects/DENATRLA/Output_Data/nmsim"),
            ("/data/DENATRLA/Output_Data/aam", "/data/DENATRLA/Output_Data/aam"),
            (
                r"C:\projects\DENATRLA\Output_Data/nmsim/scratch",
                "C:/projects/DENATRLA/Output_Data/nmsim/scratch",
            ),
        ],
    )
    def test_normalizes_separators(self, raw: str, expected: str) -> None:
        assert p.display_path(raw) == expected


class TestModelAwarePaths:
    @pytest.mark.parametrize(
        "model",
        [AcousticModel.NMSIM, AcousticModel.AAM],
    )
    def test_model_output_dirs(self, site_tree, model: AcousticModel) -> None:
        project_dir, unit, site, _year = site_tree
        site_dir = str(project_dir / f"{unit}{site}")
        suffix = "nmsim" if model is AcousticModel.NMSIM else "aam"
        assert p.model_output_dir(site_dir, model).endswith(f"Output_Data/{suffix}")

    @pytest.mark.parametrize(
        "model",
        [AcousticModel.NMSIM, AcousticModel.AAM],
    )
    def test_model_activespaces_dirs(self, site_tree, model: AcousticModel) -> None:
        project_dir, unit, site, _year = site_tree
        site_dir = str(project_dir / f"{unit}{site}")
        suffix = "nmsim" if model is AcousticModel.NMSIM else "aam"
        assert p.model_activespaces_dir(site_dir, model).endswith(
            f"Output_Data/{suffix}/ACTIVESPACES"
        )

    def test_aam_layer_dirs_only_new(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = project_dir / f"{unit}{site}"
        new_layer = site_dir / "Output_Data" / "aam" / "ACTIVESPACES" / f"{unit}{site}{year}_1000m"
        new_layer.mkdir(parents=True)
        legacy_layer = site_dir / "Output_Data" / "ACTIVESPACES" / f"{unit}{site}{year}_500m"
        legacy_layer.mkdir(parents=True)

        matches = p.activespace_layer_dirs(
            str(project_dir), unit, site, year, model=AcousticModel.AAM,
        )
        assert matches == [str(new_layer)]

    def test_aam_geojson_builds_new_path_only(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        path = p.activespace_geojson(
            str(project_dir), unit, site, year, 1000, "+", "005",
            model=AcousticModel.AAM,
        )
        assert "Output_Data/aam/ACTIVESPACES" in path
        assert path.endswith(f"{unit}{site}{year}_O_+005.geojson")

    def test_precision_recall_and_tested_points(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        site_dir = str(project_dir / f"{unit}{site}")

        assert p.precision_recall_dir(site_dir, AcousticModel.NMSIM).endswith("Output_Data/nmsim/PRECISION_RECALL")
        tested = p.tested_points_dir(site_dir, unit, site, year, 1000, AcousticModel.AAM)
        assert tested.endswith(f"Output_Data/aam/TESTED_POINTS/{unit}{site}{year}_1000m")

    def test_site_model_paths_reuses_layout(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        layout = p.SiteModelPaths.from_project(
            str(project_dir), unit, site, year, AcousticModel.AAM,
        )
        assert layout.layer_dir(1500) == p.activespace_layer_dir(
            layout.site_dir, unit, site, year, 1500, AcousticModel.AAM,
        )
        assert layout.tested_points_dir(1500) == p.tested_points_dir(
            layout.site_dir, unit, site, year, 1500, AcousticModel.AAM,
        )
        assert "Input_Data/aam" in layout.failure_hint()
        assert "TIG_TIS" not in layout.failure_hint()
        nmsim = p.SiteModelPaths.for_site(layout.site_dir, AcousticModel.NMSIM)
        assert nmsim.nmsim_scratch_glob_patterns()
        assert not layout.nmsim_scratch_glob_patterns()

    def test_precision_recall_3d_plot_model_scoped(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        path = p.precision_recall_3d_plot(
            str(project_dir), unit, site, year, model=AcousticModel.AAM, beta=1.0,
        )
        assert path.endswith(
            f"Output_Data/aam/PRECISION_RECALL/PrecisionRecallPlot_3d_{unit}{site}{year}_f1p0.png"
        )

    def test_site_model_paths_omni_aware_layer_skip(self, site_tree) -> None:
        project_dir, unit, site, year = site_tree
        layout = p.SiteModelPaths.from_project(
            str(project_dir), unit, site, year, AcousticModel.AAM,
        )
        layer = Path(layout.layer_dir(1500))
        layer.mkdir(parents=True)
        (layer / f"{layout.usy}_O_+000.geojson").write_text("{}")
        assert layout.has_layer_outputs(1500)
        assert not layout.has_layer_outputs(1500, 0.0, 2.0)


class TestLayerOutputChecks:
    def test_layer_has_activespace_outputs(self, tmp_path: Path) -> None:
        layer_dir = tmp_path / "DENATRLA2025_1500m"
        assert not p.layer_has_activespace_outputs(layer_dir)
        layer_dir.mkdir()
        assert not p.layer_has_activespace_outputs(layer_dir)
        (layer_dir / "DENATRLA2025_O_+000.geojson").write_text("{}")
        assert p.layer_has_activespace_outputs(layer_dir)

    def test_layer_has_required_omni_outputs_partial_gain(self, tmp_path: Path) -> None:
        layer_dir = tmp_path / "DENASUSH2026_1500m"
        layer_dir.mkdir()
        (layer_dir / "DENASUSH2026_O_+000.geojson").write_text("{}")
        assert not p.layer_has_required_omni_outputs(
            layer_dir, "DENASUSH2026", 0.0, 2.0,
        )
        for stem in ("O_+000", "O_+005", "O_+010", "O_+015", "O_+020"):
            (layer_dir / f"DENASUSH2026_{stem}.geojson").write_text("{}")
        assert p.layer_has_required_omni_outputs(
            layer_dir, "DENASUSH2026", 0.0, 2.0,
        )
