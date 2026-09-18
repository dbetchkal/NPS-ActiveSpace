from __future__ import annotations

from dataclasses import dataclass, field
from logging import Logger

import geopandas as gpd
import pyproj
import pyvista as pv

from nps_active_space.viz.elevation import DemElevationSampler
from nps_active_space.viz.style import VizStyle


@dataclass
class VizSession:
    """Shared PyVista scene state for all viz plotters."""

    unit: str
    site: str
    year: int
    deployment: str
    project_dir: str
    crs: str
    study_area: gpd.GeoDataFrame
    plotter: pv.Plotter
    logger: Logger
    style: VizStyle
    fill_layers: bool
    max_tracks: int
    to_wgs84: pyproj.Transformer
    dem: object | None = None
    dem_sampler: DemElevationSampler | None = None
    legend_models: list[tuple[str, str]] = field(default_factory=list)
    master_toggle_count: int = 0

    def status(self, message: str) -> None:
        self.logger.info(message)
