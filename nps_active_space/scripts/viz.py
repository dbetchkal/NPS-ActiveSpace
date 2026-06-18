from __future__ import annotations

import argparse
import glob
import os
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import pyvista as pv
import rasterio
from tqdm import tqdm
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon, box
from nps_active_space.utils.ais import query_ais_mxak
from nps_active_space.utils.helpers import (
    get_deployment,
    get_elevation,
    get_logger,
    load_annotations,
    load_DEM,
    load_layered_activespace,
    load_studyarea,
)
from nps_active_space.scripts.run_audible_transits import AudibleTransits
from nps_active_space.utils.models import Annotations, Tracks
from nps_active_space.utils.computation import NMSIM_bbox_utm
import nps_active_space.utils.config as cfg

# CLI argument validation =====================================

def parse_deployment(name: str) -> tuple[str, str, int]:
    """Parse deployment name like DENATRLA2024 into unit, site, and year."""
    if len(name) < 9:
        raise argparse.ArgumentTypeError(
            f"deployment must be at least 9 characters "
            f"(4-char unit + site code + 4-digit year), got {len(name)}: {name!r}"
        )
    unit, site, year_str = name[:4], name[4:-4], name[-4:]
    if not site:
        raise argparse.ArgumentTypeError(
            f"deployment must include a site code between unit and year, got {name!r}"
        )
    if len(year_str) != 4 or not year_str.isdigit():
        raise argparse.ArgumentTypeError(
            f"deployment year must be 4 digits, got {year_str!r} in {name!r}"
        )
    return unit, site, int(year_str)


def parse_iso_date(value: str, *, arg_name: str) -> str:
    """Validate YYYY-MM-DD date string."""
    try:
        datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{arg_name}: expected YYYY-MM-DD date, got {value!r}"
        ) from None
    return value


def parse_existing_file(path: str, *, arg_name: str) -> str:
    """Validate that a CLI path refers to an existing file."""
    file_path = Path(path)
    if not file_path.is_file():
        raise argparse.ArgumentTypeError(f"{arg_name}: file not found: {path}")
    return path


def parse_max_tracks(value: str) -> int:
    """Validate positive integer for --max-tracks."""
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--max-tracks must be a positive integer, got {value!r}"
        ) from None
    if n <= 0:
        raise argparse.ArgumentTypeError(
            f"--max-tracks must be a positive integer, got {n}"
        )
    return n


def resolve_viz_plot_flags(
    *,
    active_space: bool = False,
    annotations: bool = False,
    audible_transits: bool = False,
    vessels: bool = False,
    plot_all: bool = False,
    annotation_file: str | None = None,
    transits_pkl: str | None = None,
) -> tuple[bool, bool, bool, bool]:
    """Map CLI flags to Visualizer layer toggles."""
    do_active = active_space or plot_all
    do_annotations = annotations or plot_all or annotation_file is not None
    do_transits = audible_transits or plot_all or transits_pkl is not None
    do_vessels = vessels
    return do_active, do_annotations, do_transits, do_vessels


def utm_orientation_axes_kwargs() -> dict:
    """Orientation-marker kwargs for UTM scenes (+X east, +Y north).

    PyVista/VTK allow only one vtkOrientationMarkerWidget per renderer, so axes
    and north-arrow widgets cannot coexist as separate corner widgets.
    """
    return {
        "interactive": False,
        "line_width": 2,
        "xlabel": "E",
        "ylabel": "N",
        "zlabel": "Z",
        "viewport": (0, 0, 0.2, 0.2),
    }


# helper functions ============================================

def polygon_to_mesh(polygon: Polygon, z: float) -> pv.PolyData:
    exterior = np.array(polygon.exterior.coords)
    points = np.c_[exterior[:, 0], exterior[:, 1], np.full(len(exterior), z)]
    poly = pv.PolyData(points)
    poly.faces = np.hstack([[len(points)], np.arange(len(points))])
    return poly.triangulate()


def active_to_polys(active: gpd.GeoDataFrame) -> list[Polygon]:
    geometry = active.geometry.iloc[0]
    if isinstance(geometry, Polygon):
        polys = [geometry]
    elif isinstance(geometry, MultiPolygon):
        polys = geometry.geoms
    return polys


def active_to_linestrings(active: gpd.GeoDataFrame) -> list[LineString]:
    geometry = active.boundary.iloc[0]
    if isinstance(geometry, MultiLineString):
        linestrings = geometry.geoms
    elif isinstance(geometry, LineString):
        linestrings = [geometry]
    return linestrings


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


def _point_as_line(point: Point, span_m: float = 10.0) -> LineString:
    """Convert a single-point geometry into a short renderable segment."""
    x, y = point.x, point.y
    return LineString([(x, y), (x + span_m, y)])


def iter_plot_linestrings(geometry) -> list[LineString]:
    """Normalize annotation geometries to LineStrings PyVista can draw."""
    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, Point):
        return [_point_as_line(geometry)]
    if isinstance(geometry, LineString):
        if len(geometry.coords) < 2:
            return [_point_as_line(Point(geometry.coords[0]))]
        return [geometry]
    if isinstance(geometry, MultiLineString):
        lines: list[LineString] = []
        for part in geometry.geoms:
            lines.extend(iter_plot_linestrings(part))
        return lines
    return []


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


def densify_linestring(linestring: LineString, step_m: float) -> LineString:
    """Insert vertices along a line so elevation can be sampled between sparse AIS points."""
    if linestring.is_empty or linestring.length <= step_m:
        return linestring
    distances = np.arange(0, linestring.length, step_m)
    if distances[-1] < linestring.length:
        distances = np.append(distances, linestring.length)
    return LineString([linestring.interpolate(float(d)) for d in distances])


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


def create_polyline_3d(
    linestring: LineString,
    z: float | np.ndarray | None = None,
) -> pv.PolyData:
    coords = np.array(linestring.coords)
    xy = coords[:, :2]
    if z is not None:
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        if z_arr.size != coords.shape[0]:
            raise ValueError(
                f"z must be scalar or have one value per vertex ({coords.shape[0]}), got {z_arr.size}"
            )
        coords = np.column_stack((xy, z_arr))
    elif coords.shape[1] == 2:
        coords = np.column_stack((xy, np.zeros(coords.shape[0])))
    assert coords.shape[1] == 3
    coords[:, 2] = np.clip(np.nan_to_num(coords[:, 2], nan=0.0), 0, 10000)
    n_points = coords.shape[0]
    lines = np.hstack([[n_points], np.arange(n_points)])
    poly = pv.PolyData(coords)
    poly.lines = lines
    return poly


def track_points_to_linestring(track_points: gpd.GeoDataFrame) -> LineString:
    """Connect ordered track points into a 2D line in the plot CRS."""
    coords = [(geom.x, geom.y) for geom in track_points.geometry]
    if len(coords) < 2:
        return LineString(coords)
    return LineString(coords)




# main class =====================================================

class Visualizer():
    # color config
    activespace_color = "orange"
    mic_color = "white"
    audible_annotation_color = "deepskyblue"
    inaudible_annotation_color = "red"
    audible_transits_color = "purple"
    vessel_track_color = "cyan"
    z_scale_toggle_color = "black"  # button toggling z scale
    sea_surface_offset_m = 5.0
    sea_surface_densify_step_m = 100.0

    def __init__(
        self,
        unit: str,
        site: str,
        year: int,
        env: str,
        do_active: bool = False,
        gain: float | None = None,
        do_annots: bool = False,
        do_transits: bool = False,
        do_vessels: bool = False,
        annotation_file: str | None = None,
        audible_transits_pkl: str | None = None,
        vessel_start_date: str | None = None,
        vessel_end_date: str | None = None,
        terraced: bool = False,
        fill_layers: bool = False,
        max_tracks: int = 1000,
    ) -> None:
        # class metadata
        self.unit = unit
        self.site = site
        self.year = year
        self.deployment = f"{unit}{site}{year}"
        cfg.initialize(env)
        self.project_dir = cfg.read("project", "dir")
        self.fill_layers = fill_layers
        self.max_tracks = max_tracks
        self.logger = get_logger("VIZ", verbose=True)

        # study area and crs
        self.study_area = load_studyarea(self.project_dir, self.unit, self.site, self.year)
        self.crs = NMSIM_bbox_utm(self.study_area)
        self.study_area = self.study_area.to_crs(self.crs)
        self._to_wgs84 = pyproj.Transformer.from_crs(self.crs, "epsg:4326", always_xy=True)

        # plot each element
        self.plotter = pv.Plotter()
        self.plot_dem()
        self.plot_mic()
        if do_active:
            self.plot_activespace(terraced, gain)
        if do_annots:
            self.plot_annotations(annotation_file)
        if do_transits:
            self.plot_audible_transits(audible_transits_pkl)
        if do_vessels:
            self.plot_vessel_tracks(vessel_start_date, vessel_end_date)

        # configure plot parameters and display
        self.plotter.enable_terrain_style()
        self.setup_z_scale()
        self.plotter.add_title(f"{unit}{site}{year}", font_size=12)
        self.setup_orientation_widgets()
        self.plotter.reset_camera()
        self.plotter.camera.elevation = 30
        self.plotter.show()

    def _status(self, message: str) -> None:
        """Print and log a user-facing status line (visible when run via python)."""
        print(message, flush=True)
        self.logger.info(message)

    def _add_track_line(
        self, polyline: pv.PolyData, *, color: str, line_width: int = 2
    ):
        """Add a track/annotation polyline that reads clearly over the DEM."""
        return self.plotter.add_mesh(
            polyline,
            color=color,
            line_width=line_width,
            render_lines_as_tubes=True,
            point_size=2,
        )

    def _add_annotation_lines(
        self, polylines: list[pv.PolyData], *, color: str, line_width: int = 2
    ):
        """Add many annotation segments as one mesh (much faster than per-segment tubes)."""
        polylines = [p for p in polylines if p.n_points >= 2]
        if not polylines:
            return None
        mesh = polylines[0] if len(polylines) == 1 else pv.merge(polylines)
        return self.plotter.add_mesh(
            mesh,
            color=color,
            line_width=line_width,
            point_size=2,
            render_lines_as_tubes=False,
        )

    @staticmethod
    def _flat_sea_surface_polyline(
        linestring: LineString, offset_m: float
    ) -> pv.PolyData:
        """Sea-level track at a constant offset — no DEM resampling."""
        z = np.full(len(linestring.coords), offset_m)
        return create_polyline_3d(linestring, z=z)
    
    def plot_dem(self, show_scalar_bar: bool = False) -> None:
        # load DEM
        dem = load_DEM(self.project_dir, self.unit, self.site)
        self.dem = dem
        data = dem.read(1)
        if dem.nodata is not None:
            data[data == dem.nodata] = 0
        data[data > 9000] = 0  # higher than any elevation on earth, should be nodata
        self._dem_sampler = DemElevationSampler(dem, data, self.crs)
        
        # convert x and y coord of each DEM pixel to the CRS
        # we are not interpolating the DEM, just taking the points corresponding to each pixel
        # and putting them in the right place in the target CRS
        x = np.arange(dem.shape[1])
        y = np.arange(dem.shape[0])
        x, y = np.meshgrid(x, y)
        x_coords, y_coords = rasterio.transform.xy(dem.transform, y, x, offset="center")
        x_coords = x_coords.reshape(data.shape)
        y_coords = y_coords.reshape(data.shape)
        transformer = pyproj.Transformer.from_crs(dem.crs, self.crs, always_xy=True)
        x_coords, y_coords = transformer.transform(x_coords, y_coords)
        
        # create structured grid mesh from DEM
        mesh = pv.StructuredGrid()
        mesh.points = np.c_[x_coords.flatten(), y_coords.flatten(), data.flatten()]
        mesh.dimensions = (dem.shape[1], dem.shape[0], 1)
        mesh["elevation"] = data.flatten()

        # add mesh to plot
        self.plotter.add_mesh(mesh, scalars="elevation", cmap="gist_earth", show_scalar_bar=show_scalar_bar)

    def plot_point(self, x: float, y: float, z: float, color: str = "white") -> None:
        point = pv.PolyData(np.array([[x, y, z]]))
        self.plotter.add_mesh(point, color=color, point_size=10, render_points_as_spheres=True)

    def plot_mic(self) -> None:
        mic = get_deployment(self.project_dir, self.unit, self.site, self.year)
        mic = mic.to_crs(self.crs)
        self.plot_point(mic.x, mic.y, mic.z, self.mic_color)

    def plot_activespace(self, terraced: bool = False, gain: float | None = None) -> None:
        csv_3d_fits = f"{self.project_dir}/fits.csv"
        fit_results = pd.read_csv(csv_3d_fits, index_col="Designator")

        if gain is not None:
            print(f"Using gain {gain}dB")
        elif self.deployment in fit_results.index:
            gain = float(fit_results.loc[self.deployment, "1/3rd Octave Gain (F1)"])
        else:
            print(f"No fitted active space gain found in {csv_3d_fits}, skipping active space.")
            return
        
        # load activespace and plot the version specified by the user
        active_3d = load_layered_activespace(self.project_dir, self.unit, self.site, self.year,
                                             gain, self.crs)
        if terraced:
            self.plot_terraced_activespace(active_3d)
        else:
            self.plot_contoured_activespace(active_3d)

    def plot_contoured_activespace(self, active_3d) -> None:
        layer_checkboxes = []
        layer_callbacks = []
        for i, (active_z, active) in enumerate(active_3d.activespaces.items()):
            if not active.empty:
                checkbox, toggle_cb = self.plot_active_layer(active, active_z, i=i)
                layer_checkboxes.append(checkbox)
                layer_callbacks.append(toggle_cb)

        # make checkbox to toggle all active space layers
        def toggle_all_actives(flag):
            for box, toggle_cb in zip(layer_checkboxes, layer_callbacks):
                box.GetRepresentation().SetState(int(flag))
                toggle_cb(flag)
            self.plotter.render()
        
        self.plotter.add_checkbox_button_widget(
            callback=toggle_all_actives,
            value=True,
            position=(10, 5),
            size=35,
            color_on=self.activespace_color
        )

    def plot_active_layer(self, active_layer, elevation: float, i: int = 0):
        poly_actor = None
        if self.fill_layers:
            meshes = []
            for poly in active_to_polys(active_layer):
                meshes.append(polygon_to_mesh(poly, elevation))
                poly_data = pv.PolyData().merge(meshes)
            poly_actor = self.plotter.add_mesh(poly_data, color=self.activespace_color, opacity=0.5)

        line_actors = []
        for line in active_to_linestrings(active_layer):
            polyline = create_polyline_3d(line, z=elevation)
            actor = self.plotter.add_mesh(polyline, color=self.activespace_color, point_size=2, line_width=2)
            line_actors.append(actor)
        
        def toggle(flag):
            if poly_actor is not None:
                poly_actor.SetVisibility(flag)
            for actor in line_actors:
                actor.SetVisibility(flag)
        
        checkbox = self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(60 + 40*i,10),
            size=25,
            color_on=self.activespace_color
        )

        return checkbox, toggle

    def plot_terraced_activespace(
        self, active_3d, layer_thickness: int = 300, opacity: float = 1
    ) -> None:
        meshes = []
        layers = list(active_3d.activespaces.items())
        for i in range(len(layers)):
            active_z, active = layers[i]
            if active.empty:
                continue

            # add the vertical "walls" between the layers
            for poly in active_to_polys(active):
                hole_polys = [Polygon(hole) for hole in poly.interiors]
                for p in [poly] + hole_polys:
                    mesh = polygon_to_mesh(p, active_z - 0.5 * layer_thickness)
                    # Extrude the polygon upward by e.g. layer thickness
                    extruded = mesh.extrude([0, 0, layer_thickness], capping=False)
                    meshes.append(extruded)

            # add the flat part between the layers - the "floor" of this layer
            if i == 0:
                floor = active
            else:
                prev_active = layers[i-1][1].to_crs(self.crs)
                sym_diff = active.union_all().symmetric_difference(prev_active.union_all())
                floor = gpd.GeoDataFrame(geometry=[sym_diff], crs=self.crs)

            for poly in active_to_polys(floor):
                mesh = polygon_to_mesh(poly, active_z - 0.5*layer_thickness)
                meshes.append(mesh)
            
        # combine all meshes and add to plot
        stacked = pv.MultiBlock(meshes).combine()
        actor = self.plotter.add_mesh(stacked, color=self.activespace_color, opacity=opacity)

        def toggle(flag):
            if actor is not None:
                actor.SetVisibility(flag)
        self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(10, 5),
            size=35,
            color_on=self.activespace_color
        )
    
    def _annotation_polyline(self, linestring: LineString) -> pv.PolyData:
        """Build a 3D polyline for one annotation segment."""
        coords = np.array(linestring.coords)
        if is_surface_track(coords):
            return self._flat_sea_surface_polyline(linestring, self.sea_surface_offset_m)
        return create_polyline_3d(
            linestring,
            z=annotation_z_profile(linestring, self._dem_sampler),
        )

    def _sea_surface_polyline(self, linestring: LineString) -> pv.PolyData:
        """Build a 3D polyline that follows the local water/ground surface."""
        if is_surface_track(np.array(linestring.coords)):
            return self._flat_sea_surface_polyline(linestring, self.sea_surface_offset_m)
        line, z_vals = sea_surface_z_profile(
            linestring,
            self.dem,
            self.crs,
            offset_m=self.sea_surface_offset_m,
            densify_step_m=self.sea_surface_densify_step_m,
        )
        return create_polyline_3d(line, z=z_vals)

    def plot_annotations(self, annotation_file: str | None = None) -> None:
        if annotation_file is None:
            self._status(
                f"Loading annotations from project dir for {self.deployment} (valid only)"
            )
            annotations = load_annotations(
                self.project_dir, self.unit, self.site, self.year, only_valid=True
            )
        else:
            self._status(f"Loading annotations from {annotation_file} (valid only)")
            annotations = Annotations(annotation_file, only_valid=True)

        self._status(f"Parsed annotations: {format_annotation_summary(annotations)}")
        if annotations.empty:
            self._status("No annotations found, skipping.")
            return

        track_ids = annotations["_id"].drop_duplicates()
        if len(track_ids) > self.max_tracks:
            self._status(
                f"Sampling {self.max_tracks} of {len(track_ids)} tracks "
                f"({len(annotations)} segments before sample)"
            )
            selected_track_ids = track_ids.sample(self.max_tracks, replace=False, random_state=2)
            annotations = annotations[annotations["_id"].isin(selected_track_ids)]
            self._status(f"After sample: {format_annotation_summary(annotations)}")

        n_loaded = len(annotations)
        self._status(f"Reprojecting annotations to {self.crs} and clipping to study area")
        annotations = annotations.to_crs(self.crs)
        ann_bounds = annotations.total_bounds
        study_bounds = self.study_area.total_bounds
        self._status(f"Annotation bounds: {ann_bounds}")
        self._status(f"Study area bounds: {study_bounds}")

        n_empty = int(annotations.geometry.is_empty.sum())
        if n_empty:
            self._status(f"Dropping {n_empty} empty geometries before clip")
        annotations = annotations[~annotations.geometry.is_empty]
        annotations = annotations.clip(box(*study_bounds))
        annotations = annotations[~annotations.geometry.is_empty].explode(ignore_index=True)
        self._status(
            f"After clip: {len(annotations)} segments "
            f"({n_loaded - len(annotations)} removed outside study area)"
        )

        audible_actors = []
        inaudible_actors = []
        audible_polylines: list[pv.PolyData] = []
        inaudible_polylines: list[pv.PolyData] = []
        n_plotted = 0
        n_skipped = 0
        for _, annot in tqdm(
            annotations.iterrows(),
            total=len(annotations),
            desc="Building annotation lines",
            unit="segment",
        ):
            for line in iter_plot_linestrings(annot["geometry"]):
                polyline = self._annotation_polyline(line)
                if polyline.n_points < 2:
                    n_skipped += 1
                    continue
                n_plotted += 1
                if annot["audible"]:
                    audible_polylines.append(polyline)
                else:
                    inaudible_polylines.append(polyline)

        audible_actor = self._add_annotation_lines(
            audible_polylines, color=self.audible_annotation_color
        )
        if audible_actor is not None:
            audible_actors.append(audible_actor)
        inaudible_actor = self._add_annotation_lines(
            inaudible_polylines, color=self.inaudible_annotation_color
        )
        if inaudible_actor is not None:
            inaudible_actors.append(inaudible_actor)

        if n_plotted == 0:
            self._status(
                f"No annotation segments plotted ({n_loaded} loaded; "
                f"{n_skipped} degenerate lines skipped). "
                "Check CRS and study-area overlap."
            )
            return
        self._status(
            f"Plotted {n_plotted} annotation line(s) in "
            f"{len(audible_actors) + len(inaudible_actors)} mesh(es) "
            f"({len(audible_polylines)} audible, {len(inaudible_polylines)} inaudible; "
            f"{n_skipped} skipped)"
        )
        
        def toggle_audible(flag):
            for actor in audible_actors:
                actor.SetVisibility(flag)

        def toggle_inaudible(flag):
            for actor in inaudible_actors:
                actor.SetVisibility(flag)
        
        self.plotter.add_checkbox_button_widget(
            callback=toggle_audible,
            value=True,
            position=(10,100),
            size=25,
            color_on="deepskyblue"
        )
        self.plotter.add_checkbox_button_widget(
            callback=toggle_inaudible,
            value=True,
            position=(10,60),
            size=25,
            color_on="red"
        )
     
    def plot_audible_transits(self, audible_transits_pkl: str | None = None) -> None:
        if audible_transits_pkl:
            print(f"Loading audible transits from {audible_transits_pkl}")
            listener = AudibleTransits.from_pickle(audible_transits_pkl)
        else:
            print("Loading audible transits")
            matches = glob.glob(os.path.join(
                self.project_dir, self.unit+self.site, "Output_Data", "AUDIBLE_TRANSITS",
                f"3D*{self.year}-01-01*Active Space {self.year}*", "AudibleTransits_object.pkl"))
            if len(matches) == 0:
                print("No audible transits pkl file found found")
                return
            listener = AudibleTransits.from_pickle(matches[0])

        tracks = listener.tracks
        if tracks.empty:
            print("Audible transits is empty, skipping")
            return

        # downsample if too many
        print(f"{len(tracks)} audible transits")
        if len(tracks) > self.max_tracks:
            print("Too many, sampling")
            tracks = tracks.sample(self.max_tracks, random_state=4)
            print(f"Showing {len(tracks)} transits")
        
        tracks = tracks.to_crs(self.crs)

        # plot
        actors = []
        for _, track in tracks.iterrows():
            polyline = create_polyline_3d(track["interp_geometry"])
            actor = self.plotter.add_mesh(polyline, color=self.audible_transits_color, point_size=2, line_width=2)
            actors.append(actor)
        
        def toggle(flag):
            for actor in actors:
                actor.SetVisibility(flag)

        self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(10,180),
            size=25,
            color_on="purple"
        )

    def plot_vessel_tracks(
        self, start_date: str | None = None, end_date: str | None = None
    ) -> None:
        """Plot MXAK AIS vessel transits at the local sea surface."""
        start_date = start_date or f"{self.year}-01-01"
        end_date = end_date or f"{self.year}-12-31"

        try:
            ais_path = cfg.read("data", "ais")
        except KeyError:
            print("No [data] ais path in config, skipping vessel tracks.")
            return

        print(f"Querying AIS tracks from {start_date} to {end_date}")
        try:
            raw_tracks = query_ais_mxak(
                ais_path=ais_path,
                start_date=start_date,
                end_date=end_date,
                mask=self.study_area,
            )
        except AssertionError as exc:
            print(f"No AIS tracks loaded: {exc}")
            return
        tracks = Tracks(raw_tracks, id_col="event_id", datetime_col="TIME", z_col="altitude")
        tracks = tracks.to_crs(self.crs)

        track_ids = tracks["track_id"].drop_duplicates()
        print(f"{len(track_ids)} vessel transits ({len(tracks)} points)")
        if len(track_ids) > self.max_tracks:
            print(f"More than {self.max_tracks}, sampling")
            selected = track_ids.sample(self.max_tracks, replace=False, random_state=3)
            tracks = tracks[tracks["track_id"].isin(selected)]
            print(f"Showing {selected.nunique()} transits")

        actors = []
        for track_id, group in tracks.groupby("track_id", sort=False):
            group = group.sort_values("point_dt")
            line = track_points_to_linestring(group)
            if line.is_empty or line.length == 0:
                continue
            polyline = self._sea_surface_polyline(line)
            actor = self._add_track_line(polyline, color=self.vessel_track_color)
            actors.append(actor)

        if not actors:
            print("No vessel tracks to plot.")
            return

        def toggle(flag):
            for actor in actors:
                actor.SetVisibility(flag)

        self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(10, 220),
            size=25,
            color_on=self.vessel_track_color,
        )

    def setup_orientation_widgets(self) -> None:
        """Bottom-left E/N/Z axes (+Y = north in UTM)."""
        self.plotter.add_axes(**utm_orientation_axes_kwargs())

    def setup_z_scale(self) -> None:
        self.plotter.set_scale(1, 1, 2)  # easier to make out elevation differences

        # allow toggling though
        def toggle_z_scale(flag):
            self.plotter.set_scale(1, 1, 2 if flag else 1)

        self.plotter.add_checkbox_button_widget(
            callback=toggle_z_scale,
            value=True,
            position=(10,140),
            size=25,
            color_on=self.z_scale_toggle_color
        )



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # required args
    # we pass deployment as a single positional arg instead of -u DENA -s TRLA -y 2024,
    # because it saves a lot of keystrokes, and this script is often used frequently
    parser.add_argument(
        "deployment",
        type=parse_deployment,
        help="Deployment name, e.g. DENATRLA2024",
    )
    parser.add_argument("-e", "--environment", default="DENA_streamline",
                        help="Config environment name, e.g. DENA_streamline. Default is 'DENA_streamline'")

    # common args
    parser.add_argument("-s", "--active-space", action="store_true",
                        help="If included, load and plot the active space.")
    parser.add_argument("-g", "--gain", type=float, help="Active space gain, if not the default.")
    parser.add_argument("-a", "--annotations", action="store_true",
                        help="If included, load and plot annotations")
    parser.add_argument("-t", "--audible-transits", action="store_true",
                        help="If included, load and plot audible transits")
    parser.add_argument("-v", "--vessels", action="store_true",
                        help="If included, load and plot MXAK AIS vessel tracks at sea level")
    parser.add_argument("--all", action="store_true",
                        help="Load everything, shorthand for --active-space --annotations --audible-transits")

    # uncommon / special use case args
    parser.add_argument(
        "-m",
        "--max-tracks",
        type=parse_max_tracks,
        default=500,
        help="Maximum number of annotation tracks or audible transits to show.",
    )
    parser.add_argument(
        "--annotation-file",
        type=lambda p: parse_existing_file(p, arg_name="--annotation-file"),
        help="Path to .geojson annotations (implies -a).",
    )
    parser.add_argument(
        "--transits-pkl",
        type=lambda p: parse_existing_file(p, arg_name="--transits-pkl"),
        help="Path to .pkl audible transits (implies -t).",
    )
    parser.add_argument(
        "--start-date",
        type=lambda d: parse_iso_date(d, arg_name="--start-date"),
        help="AIS query start date (YYYY-MM-DD). Default: Jan 1 of deployment year.",
    )
    parser.add_argument(
        "--end-date",
        type=lambda d: parse_iso_date(d, arg_name="--end-date"),
        help="AIS query end date (YYYY-MM-DD). Default: Dec 31 of deployment year.",
    )
    parser.add_argument("--terraced", action="store_true",
                        help="If included, render the active space as the terraced surface instead of contours.")
    parser.add_argument("--fill-layers", action="store_true",
                        help="If included, fill the interior of each active space contour polygon.")
    
    args = parser.parse_args()

    unit, site, year = args.deployment
    print(unit, site, year)

    do_active, do_annotations, do_transits, do_vessels = resolve_viz_plot_flags(
        active_space=args.active_space,
        annotations=args.annotations,
        audible_transits=args.audible_transits,
        vessels=args.vessels,
        plot_all=args.all,
        annotation_file=args.annotation_file,
        transits_pkl=args.transits_pkl,
    )

    Visualizer(unit, site, year, args.environment,
               do_active, args.gain,
               do_annotations, do_transits, do_vessels,
               args.annotation_file, args.transits_pkl,
               args.start_date, args.end_date,
               args.terraced, args.fill_layers, args.max_tracks)