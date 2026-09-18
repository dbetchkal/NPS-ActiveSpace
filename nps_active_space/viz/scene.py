from __future__ import annotations

import numpy as np
import pyproj
import pyvista as pv
import rasterio

from nps_active_space.utils.helpers import get_deployment, load_DEM
from nps_active_space.viz.elevation import DemElevationSampler
from nps_active_space.viz.markers import utm_orientation_axes_kwargs
from nps_active_space.viz.session import VizSession


class ScenePlotter:
    def __init__(self, session: VizSession) -> None:
        self._session = session

    def plot_dem(self, show_scalar_bar: bool = False) -> None:
        session = self._session
        dem = load_DEM(session.project_dir, session.unit, session.site)
        session.dem = dem
        data = dem.read(1)
        if dem.nodata is not None:
            data[data == dem.nodata] = 0
        data[data > 9000] = 0
        session.dem_sampler = DemElevationSampler(dem, data, session.crs)

        x = np.arange(dem.shape[1])
        y = np.arange(dem.shape[0])
        x, y = np.meshgrid(x, y)
        x_coords, y_coords = rasterio.transform.xy(dem.transform, y, x, offset="center")
        x_coords = x_coords.reshape(data.shape)
        y_coords = y_coords.reshape(data.shape)
        transformer = pyproj.Transformer.from_crs(dem.crs, session.crs, always_xy=True)
        x_coords, y_coords = transformer.transform(x_coords, y_coords)

        mesh = pv.StructuredGrid()
        mesh.points = np.c_[x_coords.flatten(), y_coords.flatten(), data.flatten()]
        mesh.dimensions = (dem.shape[1], dem.shape[0], 1)
        mesh["elevation"] = data.flatten()

        session.plotter.add_mesh(
            mesh, scalars="elevation", cmap="gist_earth", show_scalar_bar=show_scalar_bar
        )

    def plot_point(self, x: float, y: float, z: float, color: str = "white") -> None:
        point = pv.PolyData(np.array([[x, y, z]]))
        self._session.plotter.add_mesh(
            point, color=color, point_size=10, render_points_as_spheres=True
        )

    def plot_mic(self) -> None:
        session = self._session
        mic = get_deployment(session.project_dir, session.unit, session.site, session.year)
        mic = mic.to_crs(session.crs)
        self.plot_point(mic.x, mic.y, mic.z, session.style.mic_color)

    def add_track_line(self, polyline: pv.PolyData, *, color: str, line_width: int = 2):
        """Add a causal track polyline over the DEM (plain lines, not tubes)."""
        return self._session.plotter.add_mesh(
            polyline,
            color=color,
            line_width=line_width,
            render_lines_as_tubes=False,
            point_size=2,
        )

    def add_annotation_lines(
        self, polylines: list[pv.PolyData], *, color: str, line_width: int = 2
    ):
        """Add many annotation segments as one mesh (much faster than per-segment tubes)."""
        polylines = [p for p in polylines if p.n_points >= 2]
        if not polylines:
            return None
        mesh = polylines[0] if len(polylines) == 1 else pv.merge(polylines)
        return self._session.plotter.add_mesh(
            mesh,
            color=color,
            line_width=line_width,
            point_size=2,
            render_lines_as_tubes=False,
        )

    def setup_orientation_widgets(self) -> None:
        """Bottom-left E/N/Z axes (+Y = north in UTM)."""
        self._session.plotter.add_axes(**utm_orientation_axes_kwargs())
