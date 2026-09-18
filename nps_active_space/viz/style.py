from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class VizStyle:
    activespace_color: str = "orange"
    nmsim_activespace_color: str = "orange"
    aam_activespace_color: str = "cyan"
    mic_color: str = "white"
    audible_annotation_color: str = "deepskyblue"
    inaudible_annotation_color: str = "red"
    audible_transits_color: str = "purple"
    vessel_track_color: str = "magenta"
    flight_track_color: str = "indigo"
    z_scale_toggle_color: str = "black"
    sea_surface_offset_m: float = 5.0
    sea_surface_densify_step_m: float = 100.0
    layer_widget_dy: int = 28
    layer_checkbox_x: int = 10
    layer_checkbox_y0: int = 220
