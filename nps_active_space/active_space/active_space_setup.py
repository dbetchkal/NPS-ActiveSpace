"""Shared active-space generator setup for NMSim and AAM propagation models."""

from __future__ import annotations

import glob
import json
import os
import time
from typing import TYPE_CHECKING

import geopandas as gpd
import pandas as pd

import nps_active_space.utils.config as cfg
from nps_active_space.active_space import ActiveSpaceGenerator
from nps_active_space.active_space.propagation_model import NMSIM_SCRATCH_SUBDIR
from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.utils.helpers import omni_to_gain
from nps_active_space.utils.legacy_nmsim_paths import LEGACY_PREDICTIONS_SUBDIR

if TYPE_CHECKING:
    from nps_active_space.utils.models import Microphone

DEFAULT_AAM_SHIM = "/usr/local/bin/aam"
DEFAULT_SRC_PT_DENSITY = 48

BATCH_RESULT_KEYS: tuple[str, ...] = (
    "Number of valid annotated segments",
    "Mean altitude",
    "KDE reduction (%)",
    "1/3rd Octave Gain (F1)",
    "F1",
)


def resolve_acoustic_model(cli_model: AcousticModel | str | None = None) -> AcousticModel:
    if cli_model is not None:
        return AcousticModel.parse(cli_model)
    env = os.environ.get("ACOUSTIC_MODEL")
    if env is not None:
        return AcousticModel.parse(env)
    return AcousticModel.NMSIM


def resolve_aam_exe(*, aam_exe: str | None = None) -> str:
    """Return the AAM launcher from an override, then ``[project] aam``, then the Docker shim."""
    if aam_exe:
        return aam_exe
    configured = str(cfg.read("project", "aam") or "").strip()
    if configured:
        return configured
    return DEFAULT_AAM_SHIM


def build_active_space_generator(
    site_dir: str,
    dem_src: str,
    study_area: gpd.GeoDataFrame,
    ambience: pd.Series,
    mic: Microphone,
    model: AcousticModel,
    *,
    nmsim_exe: str | None = None,
    aam_shim: str | None = None,
) -> ActiveSpaceGenerator:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            try:
                from nps_active_space.active_space.aam_propagation_model import (
                    AamPropagationModel,
                )
            except ModuleNotFoundError as exc:
                if exc.name != "aam_translator":
                    raise
                raise ModuleNotFoundError(
                    "AAM support requires the optional [aam] extra "
                    "(pip install -e '.[aam]')."
                ) from exc

            propagation_model = AamPropagationModel(
                site_dir, aam_shim=resolve_aam_exe(aam_exe=aam_shim),
            )
            gen = ActiveSpaceGenerator(
                NMSIM=None,
                propagation_model=propagation_model,
                root_dir=site_dir,
                study_area=study_area,
                ambience=ambience,
                dem_src=dem_src,
            )
        case AcousticModel.NMSIM:
            exe = nmsim_exe or cfg.read("project", "nmsim")
            gen = ActiveSpaceGenerator(
                NMSIM=exe,
                root_dir=site_dir,
                study_area=study_area,
                ambience=ambience,
                dem_src=dem_src,
            )
    gen.set_dem(mic)
    return gen


def omni_stem_to_gain_db(omni_stem: str) -> float:
    sign = {"+": 1, "-": -1}
    return sign[omni_stem[-4:-3]] * int(omni_stem[-3:]) / 10


def kde_reduction_percent(points_before: int, points_after: int) -> str:
    return f"{100 * (1 - (points_after / points_before))}%"


def build_batch_run_results(
    n_annotations: int,
    altitude_m: int,
    points_before_kde: int,
    points_after_kde: int,
    *,
    best_omni: str | None = None,
    f1: float | None = None,
) -> dict[str, object]:
    return {
        "Number of valid annotated segments": n_annotations,
        "Mean altitude": altitude_m,
        "KDE reduction (%)": kde_reduction_percent(points_before_kde, points_after_kde),
        "1/3rd Octave Gain (F1)": omni_to_gain(best_omni) if best_omni is not None else None,
        "F1": f1,
    }


def write_batch_run_results(results_path: str, results: dict[str, object]) -> None:
    with open(results_path, "w") as results_file:
        json.dump(results, results_file)


def precision_recall_plot_path(
    site_dir: str,
    unit: str,
    site: str,
    year,
    altitude_m: int,
    beta: float,
    model: AcousticModel,
) -> str:
    usy = p.deployment_id(unit, site, year)
    plot_name = (
        f"PrecisionRecallPlot_{usy}_{altitude_m}m_"
        f"{str(beta).replace('.', 'p')}.png"
    )
    return os.path.join(p.precision_recall_dir(site_dir, model), plot_name)


def resolve_3d_fit_gain(
    project_dir: str,
    unit: str,
    site: str,
    year,
    model: AcousticModel = AcousticModel.NMSIM,
) -> float | None:
    csv_path = p.fits_csv(project_dir)
    if not os.path.exists(csv_path):
        return None
    df = pd.read_csv(csv_path)
    usy = p.deployment_id(unit, site, year)
    model = AcousticModel.parse(model)
    if "Model" in df.columns:
        rows = df[(df["Designator"] == usy) & (df["Model"] == model)]
    elif model is AcousticModel.NMSIM:
        rows = df[df["Designator"] == usy]
    else:
        return None
    if rows.empty:
        return None
    return float(rows.iloc[-1]["1/3rd Octave Gain (F1)"])


def upsert_site_fit(
    site_dir: str,
    designator: str,
    model: AcousticModel,
    altitude_m: int,
    density: int,
    beta: float,
    best_omni: str,
    max_fbeta: float,
    precision: float,
    recall: float,
) -> str:
    """Upsert a per-layer fit row into site ``fits.csv`` (legacy/debug; production uses project ``fits.csv``)."""
    csv_path = p.site_fits_csv(site_dir)
    row = {
        "Designator": designator,
        "Model": model,
        "Altitude_m": altitude_m,
        "Density": density,
        "1/3rd Octave Gain (F1)": omni_stem_to_gain_db(best_omni),
        f"F{beta}": max_fbeta,
        "Precision": precision,
        "Recall": recall,
        "Best_omni": best_omni,
    }
    df = pd.DataFrame([row])
    if os.path.exists(csv_path):
        existing = pd.read_csv(csv_path)
        mask = (
            (existing["Designator"] == designator)
            & (existing["Model"] == model)
            & (existing["Altitude_m"] == altitude_m)
        )
        existing = existing[~mask]
        df = pd.concat([existing, df], ignore_index=True)
    os.makedirs(site_dir, exist_ok=True)
    df.to_csv(csv_path, index=False)
    return csv_path


def upsert_project_fit(
    project_dir: str,
    designator: str,
    model: AcousticModel,
    row: pd.Series | dict,
) -> str:
    csv_path = p.fits_csv(project_dir)
    fit_row = dict(row)
    fit_row["Designator"] = designator
    fit_row["Model"] = model
    df = pd.DataFrame([fit_row])
    if os.path.exists(csv_path):
        existing = pd.read_csv(csv_path)
        if "Model" in existing.columns:
            mask = (existing["Designator"] == designator) & (existing["Model"] == model)
            existing = existing[~mask]
        else:
            existing = existing[existing["Designator"] != designator]
        df = pd.concat([existing, df], ignore_index=True)
    os.makedirs(project_dir, exist_ok=True)
    df.to_csv(csv_path, index=False)
    return csv_path


def cleanup_propagation_artifacts(site_dir: str, model: AcousticModel, max_tries: int = 5) -> None:
    """Remove NMSim control/batch/TRJ/TIS scratch files. AAM run dirs are left intact."""
    if AcousticModel.parse(model) is not AcousticModel.NMSIM:
        return
    patterns = [
        f"{site_dir}/control*",
        f"{site_dir}/batch*",
        f"{site_dir}/Input_Data/03_TRAJECTORY/*.trj",
        f"{site_dir}/{LEGACY_PREDICTIONS_SUBDIR}/*.tis",
        f"{site_dir}/{NMSIM_SCRATCH_SUBDIR}/*.tis",
    ]
    try:
        for pattern in patterns:
            for file in glob.glob(pattern):
                os.remove(file)
    except OSError:
        if max_tries > 0:
            time.sleep(1)
            cleanup_propagation_artifacts(site_dir, model, max_tries - 1)
