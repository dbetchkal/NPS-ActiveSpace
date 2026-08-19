from __future__ import annotations

import geopandas as gpd
import numpy as np

from nps_active_space.viz.elevation import is_surface_track
from nps_active_space.viz.geometry import iter_plot_linestrings


def format_annotation_summary(annotations: gpd.GeoDataFrame) -> str:
    """One-line summary of a loaded annotations table for logging."""
    if annotations.empty:
        return "0 segments"
    n_surface = 0
    n_elevated = 0
    for geom in annotations.geometry:
        if geom is None or geom.is_empty:
            continue
        for line in iter_plot_linestrings(geom):
            if is_surface_track(np.array(line.coords)):
                n_surface += 1
            else:
                n_elevated += 1
    geom_types = ", ".join(
        f"{k}={v}" for k, v in annotations.geometry.geom_type.value_counts().items()
    )
    n_audible = int(annotations["audible"].sum()) if "audible" in annotations.columns else 0
    return (
        f"{len(annotations)} segments, {annotations['_id'].nunique()} tracks, "
        f"{n_audible} audible rows, {n_surface} sea-surface / {n_elevated} elevated line(s), "
        f"CRS={annotations.crs}, types: {geom_types}"
    )
