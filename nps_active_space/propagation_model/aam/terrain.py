"""AAM terrain cache, DEM preparation, and below-ground point filtering (facade)."""

from __future__ import annotations

from nps_active_space.propagation_model.aam.terrain_cache import (
    AAM_DEFAULT_FLOW_RESISTIVITY,
    AAM_INP_BASENAME,
    AAM_TERRAIN_CACHE_META,
    AOI_BOUNDS_TOLERANCE_DEG,
    AamTerrainParams,
    ensure_aam_terrain,
    log_terrain_summary,
    resolve_dem_for_aam,
    terrain_dir_for_site,
    terrain_grid_summary,
    timed_terrain_step,
)
from nps_active_space.propagation_model.aam.terrain_sampling import (
    AAM_BELOW_SURFACE_TOLERANCE_M,
    _bilinear_sample_grid,
    _elv_grid_values,
    _hop_segment_below_terrain,
    _northup_row_from_model_j,
    _split_sequential_hop_runs,
    _terrain_surface_elevation_m,
    split_below_aam_terrain,
    split_safe_aam_track_runs,
)

__all__ = [
    "AAM_BELOW_SURFACE_TOLERANCE_M",
    "AAM_DEFAULT_FLOW_RESISTIVITY",
    "AAM_INP_BASENAME",
    "AAM_TERRAIN_CACHE_META",
    "AOI_BOUNDS_TOLERANCE_DEG",
    "AamTerrainParams",
    "ensure_aam_terrain",
    "log_terrain_summary",
    "resolve_dem_for_aam",
    "split_below_aam_terrain",
    "split_safe_aam_track_runs",
    "terrain_dir_for_site",
    "terrain_grid_summary",
    "timed_terrain_step",
]
