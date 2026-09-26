from __future__ import annotations

import numpy as np
from shapely.geometry import LineString

from nps_active_space.viz.elevation import (
    annotation_z_profile,
    is_surface_track,
    sea_surface_z_profile,
)
from nps_active_space.viz.geometry import create_polyline_3d, flat_sea_surface_polyline
from nps_active_space.viz.session import VizSession


class TrackPolylineBuilder:
    def __init__(self, session: VizSession) -> None:
        self._session = session

    @staticmethod
    def flat_sea_surface_polyline(linestring: LineString, offset_m: float):
        return flat_sea_surface_polyline(linestring, offset_m)

    def annotation_polyline(self, linestring: LineString):
        """Build a 3D polyline for one annotation or flight track segment."""
        coords = np.array(linestring.coords)
        style = self._session.style
        if is_surface_track(coords):
            return flat_sea_surface_polyline(linestring, style.sea_surface_offset_m)
        return create_polyline_3d(
            linestring,
            z=annotation_z_profile(linestring, self._session.dem_sampler),
        )

    def sea_surface_polyline(self, linestring: LineString):
        """Build a 3D polyline that follows the local water/ground surface."""
        style = self._session.style
        if is_surface_track(np.array(linestring.coords)):
            return flat_sea_surface_polyline(linestring, style.sea_surface_offset_m)
        line, z_vals = sea_surface_z_profile(
            linestring,
            self._session.dem,
            self._session.crs,
            offset_m=style.sea_surface_offset_m,
            densify_step_m=style.sea_surface_densify_step_m,
        )
        return create_polyline_3d(line, z=z_vals)
