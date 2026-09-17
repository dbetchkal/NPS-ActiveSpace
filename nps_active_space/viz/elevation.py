from __future__ import annotations

import numpy as np
import pyproj
from shapely.geometry import LineString

from nps_active_space.utils.helpers import get_elevation
from nps_active_space.viz.geometry import densify_linestring


def vertex_z_from_coord(coord: tuple | np.ndarray) -> float:
    """Return vertex elevation; missing or NaN z is treated as sea level (0 m)."""
    if len(coord) < 3:
        return 0.0
    z = float(coord[2])
    return 0.0 if not np.isfinite(z) else z


def is_surface_track(coords: np.ndarray) -> bool:
    """True when every vertex is at or below sea level (vessel / surface transit)."""
    return all(vertex_z_from_coord(c) <= 0.0 for c in coords)


def is_airborne_track(coords: np.ndarray) -> bool:
    """True when every vertex has positive stored elevation (aircraft spline z)."""
    return all(vertex_z_from_coord(c) > 0.0 for c in coords)


def stored_vertex_z(linestring: LineString) -> np.ndarray:
    """Elevation from annotation geometry coordinates."""
    return np.array([vertex_z_from_coord(c) for c in linestring.coords])


class DemElevationSampler:
    """Sample elevations from an in-memory DEM band (no per-point raster I/O)."""

    def __init__(self, dem, band: np.ndarray, plot_crs: str) -> None:
        self.dem = dem
        self.band = np.asarray(band, dtype=float).copy()
        if dem.nodata is not None:
            self.band[self.band == dem.nodata] = 0.0
        self.band[self.band > 9000] = 0.0
        self._to_dem = pyproj.Transformer.from_crs(plot_crs, dem.crs, always_xy=True)

    def sample_utm_many(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        dem_x, dem_y = self._to_dem.transform(x, y)
        rows, cols = self.dem.index(dem_x, dem_y)
        rows = np.atleast_1d(rows).astype(int)
        cols = np.atleast_1d(cols).astype(int)
        in_bounds = (
            (rows >= 0)
            & (rows < self.band.shape[0])
            & (cols >= 0)
            & (cols < self.band.shape[1])
        )
        elev = np.zeros(rows.shape[0], dtype=float)
        if in_bounds.any():
            samples = self.band[rows[in_bounds], cols[in_bounds]]
            elev[in_bounds] = np.where(np.isfinite(samples), samples, 0.0)
        return elev


def annotation_z_profile(
    linestring: LineString,
    dem_sampler: DemElevationSampler | None,
    *,
    offset_m: float = 2.0,
) -> np.ndarray:
    """Per-vertex z for annotation segments.

    Aircraft splines already carry MSL elevation in geometry; only vertices at or
    below sea level need a local DEM clamp (vessels use the flat sea-surface path).
    """
    z_stored = stored_vertex_z(linestring)
    if np.all(z_stored > 0):
        return z_stored
    if dem_sampler is None:
        return z_stored
    coords = np.array(linestring.coords)
    need_dem = z_stored <= 0
    if not need_dem.any():
        return z_stored
    dem_z = dem_sampler.sample_utm_many(coords[need_dem, 0], coords[need_dem, 1])
    z_out = z_stored.copy()
    for j, i in enumerate(np.where(need_dem)[0]):
        z = z_stored[i]
        dz = float(dem_z[j])
        if z <= 0 or z < dz:
            z_out[i] = max(dz, 0.0) + offset_m
    return z_out


def safe_dem_elevation(dem, lon: float, lat: float) -> float:
    """Sample DEM elevation, returning 0 m when the point is off-raster or nodata."""
    try:
        elev = float(get_elevation(dem, lon, lat))
    except (IndexError, ValueError, TypeError):
        return 0.0
    return elev if np.isfinite(elev) else 0.0


def sea_surface_z_profile(
    linestring: LineString,
    dem,
    crs: str,
    *,
    offset_m: float = 5.0,
    densify_step_m: float = 100.0,
) -> tuple[LineString, np.ndarray]:
    """Return a densified line and z values hugging the DEM water surface (≈0 m)."""
    line = densify_linestring(linestring, densify_step_m)
    to_wgs84 = pyproj.Transformer.from_crs(crs, "epsg:4326", always_xy=True)
    z_vals = []
    for x, y in np.array(line.coords)[:, :2]:
        lon, lat = to_wgs84.transform(x, y)
        dem_z = safe_dem_elevation(dem, lon, lat)
        # NMSIM DEM uses 0 m over water; keep tracks slightly above the mesh to avoid z-fighting.
        z_vals.append(max(dem_z, 0.0) + offset_m)
    return line, np.asarray(z_vals)
