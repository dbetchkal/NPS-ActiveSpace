"""Cross-platform path helpers for the NPS-ActiveSpace project directory layout.

All filesystem paths in the pipeline should be built through these helpers (or
``os.path.join`` directly) rather than hard-coded backslashes or ``f"{a}/{b}"``
strings. Forward slashes often work on Windows too, but ``os.path.join`` is
explicit and keeps globs working on Linux/Mac.

Prefer :class:`SiteModelPaths` when a caller already knows site + model (+
deployment). The module-level functions remain for one-off lookups.

Use :func:`display_path` when printing or logging filesystem paths so output
uses forward slashes on all platforms.
"""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from pathlib import Path

from nps_active_space.utils.enums import AcousticModel

NMSIM_OUTPUT_SUBDIR = "Output_Data/nmsim"
NMSIM_PREDICTIONS_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/predictions"
NMSIM_SCRATCH_SUBDIR = f"{NMSIM_OUTPUT_SUBDIR}/scratch"
AAM_MODEL_SLUG = "aam"
AAM_INPUT_SUBDIR = f"Input_Data/{AAM_MODEL_SLUG}"
AAM_INPUT_SUBDIR_LEGACY = "Input_Data/AAM"
AAM_TERRAIN_SUBDIR = f"{AAM_INPUT_SUBDIR}/terrain"
AAM_OUTPUT_SUBDIR = f"Output_Data/{AAM_MODEL_SLUG}"
AAM_PREDICTIONS_SUBDIR = f"{AAM_OUTPUT_SUBDIR}/predictions"
AAM_RUNS_SUBDIR = f"{AAM_OUTPUT_SUBDIR}/runs"
AAM_RUN_LOG_FILENAME = "active_space.log"
from nps_active_space.utils.legacy_nmsim_paths import (
    LEGACY_PREDICTIONS_SUBDIR,
    NMSIM_ACTIVESPACES_SUBDIR,
    is_standard_altitude_layer_dir,
    resolve_activespace_geojson as _resolve_activespace_geojson,
    resolve_activespace_layer_dirs as _resolve_activespace_layer_dirs,
    resolve_nmsim_activespaces_dir,
)


def join(*parts: str) -> str:
    return os.path.join(*parts)


def display_path(path: str | Path) -> str:
    """Format a filesystem path for logs and user-facing messages.

    Uses forward slashes so paths stay copy-pasteable in Explorer, terminals,
    and docs regardless of host OS. Do not use for ``open()`` or other I/O.
    """
    return os.fspath(Path(path).expanduser()).replace("\\", "/")


def deployment_id(unit: str, site: str, year) -> str:
    return f"{unit}{site}{year}"


def site_dir(project_dir: str, unit: str, site: str) -> str:
    return join(project_dir, f"{unit}{site}")


def site_path(project_dir: str, unit: str, site: str, *parts: str) -> str:
    return join(project_dir, f"{unit}{site}", *parts)


def input_data_dir(project_dir: str, unit: str, site: str) -> str:
    return site_path(project_dir, unit, site, "Input_Data")


def output_data_dir(project_dir: str, unit: str, site: str) -> str:
    return site_path(project_dir, unit, site, "Output_Data")


def model_output_dir(site_dir_path: str, model: AcousticModel) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.NMSIM:
            subdir = NMSIM_OUTPUT_SUBDIR
        case AcousticModel.AAM:
            subdir = AAM_OUTPUT_SUBDIR
    return join(site_dir_path, subdir)


def model_activespaces_dir(site_dir_path: str, model: AcousticModel) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            return join(site_dir_path, AAM_OUTPUT_SUBDIR, "ACTIVESPACES")
        case AcousticModel.NMSIM:
            return join(site_dir_path, NMSIM_ACTIVESPACES_SUBDIR)


def layer_has_activespace_outputs(layer_dir: str | Path) -> bool:
    """True when a layer folder already has omni geojson outputs."""
    layer_path = Path(layer_dir)
    if not layer_path.is_dir():
        return False
    return any(layer_path.glob("*_O_*.geojson"))


def layer_has_required_omni_outputs(
    layer_dir: str | Path,
    usy: str,
    omni_min: float,
    omni_max: float,
    step_db: float = 0.5,
) -> bool:
    """True when every omni in ``[omni_min, omni_max]`` has a layer geojson."""
    from nps_active_space.utils.helpers import get_omni_sources

    layer_path = Path(layer_dir)
    if not layer_path.is_dir():
        return False
    for src in get_omni_sources(omni_min, omni_max, step_db):
        stem = Path(src).stem
        if not (layer_path / f"{usy}_{stem}.geojson").is_file():
            return False
    return True


@dataclass(frozen=True)
class SiteModelPaths:
    """Model-scoped site layout. Build once per generate/batch/fit call."""

    site_dir: str
    model: AcousticModel
    unit: str
    site: str
    year: int | str

    @classmethod
    def from_project(
        cls,
        project_dir: str,
        unit: str,
        site: str,
        year: int | str,
        model: AcousticModel | str,
    ) -> SiteModelPaths:
        return cls(
            site_dir=site_dir(project_dir, unit, site),
            model=AcousticModel.parse(model),
            unit=unit,
            site=site,
            year=year,
        )

    @classmethod
    def for_site(cls, site_dir_path: str, model: AcousticModel | str) -> SiteModelPaths:
        """Layout for site-level dirs that do not need a deployment id."""
        return cls(
            site_dir=site_dir_path,
            model=AcousticModel.parse(model),
            unit="",
            site="",
            year="",
        )

    @property
    def usy(self) -> str:
        return deployment_id(self.unit, self.site, self.year)

    @property
    def output_dir(self) -> str:
        return model_output_dir(self.site_dir, self.model)

    @property
    def activespaces_dir(self) -> str:
        return model_activespaces_dir(self.site_dir, self.model)

    @property
    def predictions_dir(self) -> str:
        match self.model:
            case AcousticModel.AAM:
                return join(self.site_dir, AAM_PREDICTIONS_SUBDIR)
            case AcousticModel.NMSIM:
                return join(self.site_dir, NMSIM_PREDICTIONS_SUBDIR)

    @property
    def scratch_dir(self) -> str:
        match self.model:
            case AcousticModel.AAM:
                return join(self.site_dir, AAM_RUNS_SUBDIR)
            case AcousticModel.NMSIM:
                return join(self.site_dir, NMSIM_SCRATCH_SUBDIR)

    @property
    def precision_recall_dir(self) -> str:
        return join(self.output_dir, "PRECISION_RECALL")

    @property
    def ambience_dir(self) -> str:
        return join(self.site_dir, "Output_Data", "AMBIENCE")

    @property
    def aam_input_dir(self) -> str:
        return join(self.site_dir, AAM_INPUT_SUBDIR)

    @property
    def trajectory_dir(self) -> str:
        return join(self.site_dir, "Input_Data", "03_TRAJECTORY")

    def layer_dir(self, altitude_m: int) -> str:
        return join(self.activespaces_dir, f"{self.usy}_{altitude_m}m")

    def tested_points_dir(self, altitude_m: int) -> str:
        return join(self.output_dir, "TESTED_POINTS", f"{self.usy}_{altitude_m}m")

    def precision_recall_plot(self, altitude_m: int, beta: float) -> str:
        beta_str = str(beta).replace(".", "p")
        return join(
            self.precision_recall_dir,
            f"PrecisionRecallPlot_{self.usy}_{altitude_m}m_{beta_str}.png",
        )

    def has_layer_outputs(
        self,
        altitude_m: int,
        omni_min: float | None = None,
        omni_max: float | None = None,
        omni_step_db: float = 0.5,
    ) -> bool:
        if omni_min is not None and omni_max is not None:
            return layer_has_required_omni_outputs(
                self.layer_dir(altitude_m),
                self.usy,
                omni_min,
                omni_max,
                omni_step_db,
            )
        return layer_has_activespace_outputs(self.layer_dir(altitude_m))

    def failure_hint(self) -> str:
        """Model-specific paths to inspect after a failed layer run."""
        match self.model:
            case AcousticModel.AAM:
                return (
                    f"check {display_path(self.aam_input_dir)}, "
                    f"{display_path(self.scratch_dir)}, "
                    f"and {AAM_RUN_LOG_FILENAME}"
                )
            case AcousticModel.NMSIM:
                return (
                    f"check {display_path(self.trajectory_dir)}, "
                    f"{display_path(self.predictions_dir)}, "
                    f"and {display_path(self.scratch_dir)}"
                )

    def nmsim_scratch_glob_patterns(self) -> list[str]:
        """Glob patterns for ``--cleanup-nmsim-scratch``. Empty for AAM."""
        if self.model is not AcousticModel.NMSIM:
            return []
        return [
            join(self.site_dir, "control*"),
            join(self.site_dir, "batch*"),
            join(self.trajectory_dir, "*.trj"),
            join(self.site_dir, LEGACY_PREDICTIONS_SUBDIR, "*.tis"),
            join(self.site_dir, NMSIM_SCRATCH_SUBDIR, "*.tis"),
        ]


def activespaces_dir(project_dir: str, unit: str, site: str) -> str:
    """Legacy-compatible ACTIVESPACES root (prefers new NMSim layout when populated)."""
    return resolve_nmsim_activespaces_dir(site_dir(project_dir, unit, site))


def activespace_layer_dir(
    site_dir_path: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    model: AcousticModel,
) -> str:
    return SiteModelPaths(
        site_dir_path, AcousticModel.parse(model), unit, site, year,
    ).layer_dir(altitude_m)


def activespace_layer_dirs(
    project_dir: str,
    unit: str,
    site: str,
    year,
    model: AcousticModel = AcousticModel.NMSIM,
) -> list[str]:
    """Glob paths like ``.../ACTIVESPACES/DENATRLA2025_1000m``."""
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            usy = deployment_id(unit, site, year)
            pattern = join(
                model_activespaces_dir(
                    site_dir(project_dir, unit, site), AcousticModel.AAM,
                ),
                f"{usy}_*m",
            )
            matches = glob.glob(pattern)
            return sorted(
                d for d in matches if is_standard_altitude_layer_dir(d, usy)
            )
        case AcousticModel.NMSIM:
            return _resolve_activespace_layer_dirs(project_dir, unit, site, year)


def activespace_geojson(
    project_dir: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    gain_sign: str,
    gain_string: str,
    model: AcousticModel = AcousticModel.NMSIM,
) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            usy = deployment_id(unit, site, year)
            return join(
                model_activespaces_dir(
                    site_dir(project_dir, unit, site), AcousticModel.AAM,
                ),
                f"{usy}_{altitude_m}m",
                f"{usy}_O_{gain_sign}{gain_string}.geojson",
            )
        case AcousticModel.NMSIM:
            return _resolve_activespace_geojson(
                project_dir, unit, site, year, altitude_m, gain_sign, gain_string,
            )


def precision_recall_dir(site_dir_path: str, model: AcousticModel) -> str:
    return join(model_output_dir(site_dir_path, model), "PRECISION_RECALL")


def site_fits_csv(site_dir_path: str) -> str:
    return join(site_dir_path, "fits.csv")


def tested_points_dir(
    site_dir_path: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    model: AcousticModel,
) -> str:
    return SiteModelPaths(
        site_dir_path, AcousticModel.parse(model), unit, site, year,
    ).tested_points_dir(altitude_m)


def annotation_files(project_dir: str, unit: str, site: str, year) -> list[str]:
    """Annotation geojson files saved in the site directory root (ground-truthing output)."""
    usy = deployment_id(unit, site, year)
    pattern = join(site_dir(project_dir, unit, site), f"{usy}*saved_annotations*.geojson")
    return glob.glob(pattern)


def study_area_shapefile(project_dir: str, unit: str, site: str) -> str:
    site = site_dir(project_dir, unit, site)
    matches = glob.glob(join(site, "*study*.shp"))
    if not matches:
        matches = glob.glob(join(site, f"{unit}{site}*study*area*.shp"))
    if not matches:
        raise FileNotFoundError(f"No study area shapefile under {display_path(site)}")
    return matches[0]


def dem_raster(project_dir: str, unit: str, site: str) -> str:
    elevation_dir = join(input_data_dir(project_dir, unit, site), "01_ELEVATION")
    matches = glob.glob(join(elevation_dir, "elevation_m_nad83_utm*.tif"))
    if not matches:
        raise FileNotFoundError(
            f"No elevation_m_nad83_utm*.tif under {display_path(elevation_dir)}"
        )
    return matches[0]


def site_file(project_dir: str, unit: str, site: str, year) -> str:
    return join(input_data_dir(project_dir, unit, site), "05_SITES",
                f"{unit}{site}{year}.sit")


def fits_csv(project_dir: str) -> str:
    return join(project_dir, "fits.csv")


def precision_recall_3d_plot(
    project_dir: str,
    unit: str,
    site: str,
    year,
    model: AcousticModel = AcousticModel.NMSIM,
    beta: float = 1.0,
) -> str:
    usy = deployment_id(unit, site, year)
    beta_str = str(beta).replace(".", "p")
    site_dir_path = site_dir(project_dir, unit, site)
    return join(
        model_output_dir(site_dir_path, model),
        "PRECISION_RECALL",
        f"PrecisionRecallPlot_3d_{usy}_f{beta_str}.png",
    )


def altitude_histogram_plot(project_dir: str, unit: str, site: str, year) -> str:
    usy = deployment_id(unit, site, year)
    return join(site_dir(project_dir, unit, site), f"Altitude_Histogram_{usy}.png")


def all_annotation_files(project_dir: str) -> list[str]:
    """Find annotation geojson files under each site directory in ``project_dir``."""
    return glob.glob(join(project_dir, "*", "*saved_annotations*.geojson"))
