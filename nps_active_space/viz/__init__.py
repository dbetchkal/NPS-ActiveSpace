"""3D deployment visualization (PyVista)."""

from nps_active_space.viz.annotations import format_annotation_summary
from nps_active_space.viz.cli import (
    main,
    parse_existing_file,
    parse_iso_date,
    parse_max_tracks,
    resolve_track_source_args,
    resolve_viz_plot_flags,
)
from nps_active_space.viz.markers import utm_orientation_axes_kwargs
from nps_active_space.viz.elevation import (
    DemElevationSampler,
    annotation_z_profile,
    is_airborne_track,
    is_surface_track,
    sea_surface_z_profile,
    vertex_z_from_coord,
)
from nps_active_space.viz.geometry import (
    create_polyline_3d,
    densify_linestring,
    flat_sea_surface_polyline,
    iter_plot_linestrings,
    track_points_to_linestring,
)
from nps_active_space.viz.visualizer import Visualizer

__all__ = [
    "DemElevationSampler",
    "Visualizer",
    "annotation_z_profile",
    "create_polyline_3d",
    "densify_linestring",
    "format_annotation_summary",
    "is_airborne_track",
    "is_surface_track",
    "iter_plot_linestrings",
    "main",
    "parse_existing_file",
    "parse_iso_date",
    "parse_max_tracks",
    "resolve_track_source_args",
    "resolve_viz_plot_flags",
    "sea_surface_z_profile",
    "track_points_to_linestring",
    "utm_orientation_axes_kwargs",
    "vertex_z_from_coord",
]
