"""CSV cache of 1/3-octave predictions in dB, keyed by altitude, omni, and heading."""

from __future__ import annotations

import logging
import os
from collections.abc import Callable

import geopandas as gpd
import pandas as pd

from nps_active_space.utils.paths import display_path

logger = logging.getLogger(__name__)


def prediction_cache_csv_path(
    root_dir: str,
    predictions_subdir: str,
    altitude_m: int,
    omni_stem: str,
    heading: int,
) -> str:
    """Path to the merged prediction cache CSV for one alt/omni/heading combination."""
    cache_dir = os.path.join(root_dir, predictions_subdir)
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, f"{altitude_m}m_{omni_stem}_{heading}deg.csv")

_REQUIRED_CACHE_COLUMNS = {"Xpos", "Ypos", "A"}
_NO_SOUND_DB = -99.9
_CACHE_UNITS_HEADER = "# units=dB\n"


def cache_failure_reason(csv_filename: str) -> str | None:
    """Why an on-disk cache cannot be used, or None if missing or readable."""
    if not os.path.exists(csv_filename):
        return None
    if os.path.getsize(csv_filename) == 0:
        return "file is empty (0 bytes)"
    with open(csv_filename, encoding="utf-8") as cache_file:
        first_line = cache_file.readline()
    if not first_line.lstrip().startswith("# units=dB"):
        return "missing # units=dB header"
    try:
        preview = pd.read_csv(csv_filename, comment="#", nrows=1)
    except pd.errors.EmptyDataError:
        return "file has no parseable CSV content"
    if preview.empty:
        return "file has a header but no data rows"
    missing = sorted(_REQUIRED_CACHE_COLUMNS - set(preview.columns))
    if missing:
        return f"missing required columns: {missing}"
    return None


def cache_is_readable(csv_filename: str) -> bool:
    return os.path.exists(csv_filename) and cache_failure_reason(csv_filename) is None


def load_prediction_cache(
    source_pts: gpd.GeoDataFrame,
    csv_filename: str,
    altitude_m: int,
) -> tuple[pd.DataFrame, pd.DataFrame, gpd.GeoDataFrame]:
    """Split ``source_pts`` into cached spectral rows vs points still needing predict().

    Returns
    -------
    cache_df_all, cache_hits, new_pts
        All rows in the file, hits matching ``source_pts``, and unmatched source points.
        Unreadable files are deleted and treated as a miss.
    """
    if not cache_is_readable(csv_filename):
        reason = cache_failure_reason(csv_filename)
        if reason is not None:
            logger.warning(
                "Ignoring unreadable prediction cache at %s (%s). "
                "Removing file; points will be recomputed and cache rewritten after predict.",
                display_path(csv_filename),
                reason,
            )
            try:
                os.remove(csv_filename)
            except OSError as exc:
                logger.warning(
                    "Could not remove unreadable prediction cache %s: %s",
                    display_path(csv_filename),
                    exc,
                )
        return pd.DataFrame(), pd.DataFrame(), source_pts

    cache_df_all = pd.read_csv(csv_filename, comment="#")
    # Zpos is omitted on disk because it is constant within a cache file.
    cache_df_all = cache_df_all.fillna(_NO_SOUND_DB).astype("float64")
    cache_df_all["Zpos"] = altitude_m

    source_idx = pd.MultiIndex.from_frame(pd.DataFrame({
        "Xpos": source_pts.geometry.x,
        "Ypos": source_pts.geometry.y,
    }))
    prev_idx = pd.MultiIndex.from_frame(cache_df_all[["Xpos", "Ypos"]])

    cache_hits = cache_df_all[prev_idx.isin(source_idx)].drop_duplicates(["Xpos", "Ypos"])
    new_pts = source_pts[~source_idx.isin(prev_idx)].drop_duplicates("geometry")
    assert len(new_pts) + len(cache_hits) == len(source_pts)

    return cache_df_all, cache_hits, new_pts


def save_prediction_cache(cache_df_all: pd.DataFrame, csv_filename: str) -> None:
    """Write spectral predictions in dB (no Zpos; -99.9 dB as NA)."""
    if cache_df_all.empty:
        return

    cache_df_all = cache_df_all.drop_duplicates(subset=["Xpos", "Ypos"]).copy()

    db_cols = cache_df_all.loc[:, "A":"12500"].columns
    cache_df_all[db_cols] = cache_df_all[db_cols].round(1).replace(_NO_SOUND_DB, pd.NA)
    cache_df_all.drop("Zpos", axis=1, inplace=True, errors="ignore")
    with open(csv_filename, "w", encoding="utf-8", newline="") as cache_file:
        cache_file.write(_CACHE_UNITS_HEADER)
        cache_df_all.to_csv(cache_file, index=False)


def source_pts_missing_predictions(
    source_pts: gpd.GeoDataFrame,
    pred_df: pd.DataFrame,
) -> gpd.GeoDataFrame:
    """Points sent to predict() that have no row in the returned prediction DataFrame."""
    if len(source_pts) == 0:
        return source_pts.iloc[0:0]
    if len(pred_df) == 0:
        return source_pts

    source_idx = pd.MultiIndex.from_frame(pd.DataFrame({
        "Xpos": source_pts.geometry.x,
        "Ypos": source_pts.geometry.y,
    }))
    pred_idx = pd.MultiIndex.from_frame(pred_df[["Xpos", "Ypos"]].drop_duplicates())
    return source_pts[~source_idx.isin(pred_idx)].drop_duplicates("geometry")


def predict_with_cache(
    predict_fn: Callable[[gpd.GeoDataFrame], pd.DataFrame],
    source_pts: gpd.GeoDataFrame,
    csv_filename: str,
    altitude_m: int,
    job_name: str,
) -> tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """Call ``predict_fn`` only for uncached points; merge and persist the CSV.

    Returns spectral rows for ``source_pts`` that succeeded, plus points the model
    skipped (caller typically marks those inaudible).
    """
    cache_df_all, cache_hits, new_pts = load_prediction_cache(
        source_pts, csv_filename, altitude_m,
    )

    if len(new_pts) == 0:
        logger.info(
            "%s: prediction cache hit for all %d point(s)",
            job_name,
            len(source_pts),
        )
        return cache_hits, source_pts.iloc[0:0]

    new_df = predict_fn(new_pts)
    failed_pts = source_pts_missing_predictions(new_pts, new_df)
    if len(failed_pts) > 0:
        logger.warning(
            "%s: marking %d point(s) inaudible after predict failure/skip",
            job_name,
            len(failed_pts),
        )
    cache_hits = pd.concat([cache_hits, new_df], ignore_index=True)
    cache_hits = cache_hits.drop_duplicates(subset=["Xpos", "Ypos"])
    cache_df_all = pd.concat([cache_df_all, new_df], ignore_index=True)
    save_prediction_cache(cache_df_all, csv_filename)
    return cache_hits, failed_pts


__all__ = [
    "cache_failure_reason",
    "cache_is_readable",
    "load_prediction_cache",
    "prediction_cache_csv_path",
    "predict_with_cache",
    "save_prediction_cache",
    "source_pts_missing_predictions",
]
