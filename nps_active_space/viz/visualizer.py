from __future__ import annotations

import glob
import configparser
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import pyvista as pv
import rasterio
from shapely.geometry import LineString, Polygon, box
from tqdm import tqdm

import nps_active_space.utils.config as cfg
from nps_active_space.scripts.run_audible_transits import AudibleTransits
from nps_active_space.utils.computation import NMSIM_bbox_utm
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.helpers import (
    get_deployment,
    get_logger,
    load_annotations,
    load_DEM,
    load_layered_activespace,
    load_studyarea,
)
from nps_active_space.utils.load_tracks import load_tracks
from nps_active_space.utils.models import Annotations
from nps_active_space.viz.annotations import format_annotation_summary
from nps_active_space.viz.markers import (
    WINDOW_TITLE,
    apply_window_icon,
    register_windows_taskbar_identity,
    utm_orientation_axes_kwargs,
)
from nps_active_space.viz.elevation import (
    DemElevationSampler,
    annotation_z_profile,
    is_surface_track,
    sea_surface_z_profile,
)
from nps_active_space.viz.geometry import (
    active_to_linestrings,
    active_to_polys,
    create_polyline_3d,
    flat_sea_surface_polyline,
    iter_plot_linestrings,
    polygon_to_mesh,
    track_points_to_linestring,
)


class Visualizer:
    # color config
    activespace_color = "orange"
    mic_color = "white"
    audible_annotation_color = "deepskyblue"
    inaudible_annotation_color = "red"
    audible_transits_color = "purple"
    vessel_track_color = "cyan"
    flight_track_color = "yellow"
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
        track_source: TrackSource | None = None,
        annotation_file: str | None = None,
        audible_transits_pkl: str | None = None,
        track_start_date: str | None = None,
        track_end_date: str | None = None,
        terraced: bool = False,
        fill_layers: bool = False,
        max_tracks: int = 1000,
    ) -> None:
        self.unit = unit
        self.site = site
        self.year = year
        self.deployment = f"{unit}{site}{year}"
        cfg.initialize(env)
        self.project_dir = cfg.read("project", "dir")
        self.fill_layers = fill_layers
        self.max_tracks = max_tracks
        self.logger = get_logger("VIZ", verbose=True)

        self.study_area = load_studyarea(self.project_dir, self.unit, self.site, self.year)
        self.crs = NMSIM_bbox_utm(self.study_area)
        self.study_area = self.study_area.to_crs(self.crs)
        self._to_wgs84 = pyproj.Transformer.from_crs(self.crs, "epsg:4326", always_xy=True)

        register_windows_taskbar_identity()
        self.plotter = pv.Plotter(title=WINDOW_TITLE)
        self.plot_dem()
        self.plot_mic()
        if do_active:
            self.plot_activespace(terraced, gain)
        if do_annots:
            self.plot_annotations(annotation_file)
        if do_transits:
            self.plot_audible_transits(audible_transits_pkl)
        if track_source is not None:
            self.plot_tracks(track_source, track_start_date, track_end_date)

        self.plotter.enable_terrain_style()
        self.setup_z_scale()
        self.plotter.add_title(f"{unit}{site}{year}", font_size=12)
        self.setup_orientation_widgets()
        self.plotter.reset_camera()
        self.plotter.camera.elevation = 30
        apply_window_icon(self.plotter)
        self.plotter.show()

    def _status(self, message: str) -> None:
        """Log a user-facing status line to the console (via get_logger StreamHandler)."""
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
        return flat_sea_surface_polyline(linestring, offset_m)

    def plot_dem(self, show_scalar_bar: bool = False) -> None:
        dem = load_DEM(self.project_dir, self.unit, self.site)
        self.dem = dem
        data = dem.read(1)
        if dem.nodata is not None:
            data[data == dem.nodata] = 0
        data[data > 9000] = 0
        self._dem_sampler = DemElevationSampler(dem, data, self.crs)

        x = np.arange(dem.shape[1])
        y = np.arange(dem.shape[0])
        x, y = np.meshgrid(x, y)
        x_coords, y_coords = rasterio.transform.xy(dem.transform, y, x, offset="center")
        x_coords = x_coords.reshape(data.shape)
        y_coords = y_coords.reshape(data.shape)
        transformer = pyproj.Transformer.from_crs(dem.crs, self.crs, always_xy=True)
        x_coords, y_coords = transformer.transform(x_coords, y_coords)

        mesh = pv.StructuredGrid()
        mesh.points = np.c_[x_coords.flatten(), y_coords.flatten(), data.flatten()]
        mesh.dimensions = (dem.shape[1], dem.shape[0], 1)
        mesh["elevation"] = data.flatten()

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
            self._status(f"Using gain {gain}dB")
        elif self.deployment in fit_results.index:
            gain = float(fit_results.loc[self.deployment, "1/3rd Octave Gain (F1)"])
        else:
            self._status(
                f"No fitted active space gain found in {csv_3d_fits}, skipping active space."
            )
            return

        active_3d = load_layered_activespace(
            self.project_dir, self.unit, self.site, self.year, gain, self.crs
        )
        if not active_3d.activespaces:
            self._status(
                f"No active space geometry loaded for {self.deployment} at {gain}dB; skipping."
            )
            return
        self._status(
            f"Loaded active space layers (m): {sorted(active_3d.activespaces.keys())} at {gain}dB"
        )
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
            color_on=self.activespace_color,
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
            actor = self.plotter.add_mesh(
                polyline,
                color=self.activespace_color,
                line_width=4,
                render_lines_as_tubes=True,
            )
            line_actors.append(actor)

        def toggle(flag):
            if poly_actor is not None:
                poly_actor.SetVisibility(flag)
            for actor in line_actors:
                actor.SetVisibility(flag)

        checkbox = self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(60 + 40 * i, 10),
            size=25,
            color_on=self.activespace_color,
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

            for poly in active_to_polys(active):
                hole_polys = [Polygon(hole) for hole in poly.interiors]
                for p in [poly] + hole_polys:
                    mesh = polygon_to_mesh(p, active_z - 0.5 * layer_thickness)
                    extruded = mesh.extrude([0, 0, layer_thickness], capping=False)
                    meshes.append(extruded)

            if i == 0:
                floor = active
            else:
                prev_active = layers[i - 1][1].to_crs(self.crs)
                sym_diff = active.union_all().symmetric_difference(prev_active.union_all())
                floor = gpd.GeoDataFrame(geometry=[sym_diff], crs=self.crs)

            for poly in active_to_polys(floor):
                mesh = polygon_to_mesh(poly, active_z - 0.5 * layer_thickness)
                meshes.append(mesh)

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
            color_on=self.activespace_color,
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
            position=(10, 100),
            size=25,
            color_on="deepskyblue",
        )
        self.plotter.add_checkbox_button_widget(
            callback=toggle_inaudible,
            value=True,
            position=(10, 60),
            size=25,
            color_on="red",
        )

    def plot_audible_transits(self, audible_transits_pkl: str | None = None) -> None:
        if audible_transits_pkl:
            self._status(f"Loading audible transits from {audible_transits_pkl}")
            listener = AudibleTransits.from_pickle(audible_transits_pkl)
        else:
            self._status("Loading audible transits")
            matches = glob.glob(
                os.path.join(
                    self.project_dir,
                    self.unit + self.site,
                    "Output_Data",
                    "AUDIBLE_TRANSITS",
                    f"3D*{self.year}-01-01*Active Space {self.year}*",
                    "AudibleTransits_object.pkl",
                )
            )
            if len(matches) == 0:
                self._status("No audible transits pkl file found")
                return
            listener = AudibleTransits.from_pickle(matches[0])

        tracks = listener.tracks
        if tracks.empty:
            self._status("Audible transits is empty, skipping")
            return

        self._status(f"{len(tracks)} audible transits")
        if len(tracks) > self.max_tracks:
            self._status("Too many, sampling")
            tracks = tracks.sample(self.max_tracks, random_state=4)
            self._status(f"Showing {len(tracks)} transits")

        tracks = tracks.to_crs(self.crs)

        actors = []
        for _, track in tracks.iterrows():
            polyline = create_polyline_3d(track["interp_geometry"])
            actor = self.plotter.add_mesh(
                polyline, color=self.audible_transits_color, point_size=2, line_width=2
            )
            actors.append(actor)

        def toggle(flag):
            for actor in actors:
                actor.SetVisibility(flag)

        self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(10, 180),
            size=25,
            color_on="purple",
        )

    def plot_tracks(
        self,
        source: TrackSource,
        start_date: str | None = None,
        end_date: str | None = None,
    ) -> None:
        """Plot causal tracks from GPS, ADSB, or MXAK AIS for the study window.

        Uses the same ``load_tracks`` loader as ground truthing but draws raw point
        sequences (not annotation splines) and does not apply clock-drift correction.
        Prefer explicit ``--start-date`` / ``--end-date``; default is the full year.
        """
        start_date = start_date or f"{self.year}-01-01"
        end_date = end_date or f"{self.year}-12-31"
        microphone = get_deployment(self.project_dir, self.unit, self.site, self.year)

        self._status(f"Querying {source} tracks from {start_date} to {end_date}")
        try:
            loaded = load_tracks(
                source,
                start_date=start_date,
                end_date=end_date,
                study_area=self.study_area,
                microphone=microphone,
                include_faa_paths=False,
            )
        except (configparser.NoSectionError, configparser.NoOptionError) as exc:
            self._status(f"No {source} tracks loaded: missing config option {exc!r}")
            return
        except KeyError as exc:
            self._status(f"No {source} tracks loaded: {exc}")
            return
        except (AssertionError, ValueError) as exc:
            self._status(f"No {source} tracks loaded: {exc}")
            return

        tracks = loaded.tracks.to_crs(self.crs)
        if tracks.empty:
            self._status(f"No {source} tracks loaded.")
            return

        track_ids = tracks["track_id"].drop_duplicates()
        self._status(f"{len(track_ids)} {source} tracks ({len(tracks)} points)")
        if len(track_ids) > self.max_tracks:
            self._status(f"More than {self.max_tracks}, sampling")
            selected = track_ids.sample(self.max_tracks, replace=False, random_state=3)
            tracks = tracks[tracks["track_id"].isin(selected)]
            self._status(f"Showing {selected.nunique()} tracks")

        color = (
            self.vessel_track_color
            if source is TrackSource.AIS
            else self.flight_track_color
        )
        actors = []
        for _track_id, group in tracks.groupby("track_id", sort=False):
            group = group.sort_values("point_dt")
            line = track_points_to_linestring(
                group,
                include_z=source is not TrackSource.AIS,
            )
            if line.is_empty or line.length == 0:
                continue
            match source:
                case TrackSource.AIS:
                    polyline = self._sea_surface_polyline(line)
                case TrackSource.ADSB | TrackSource.GPS:
                    polyline = self._annotation_polyline(line)
                case _:
                    raise ValueError(f"Unknown track source: {source}")
            actor = self._add_track_line(polyline, color=color)
            actors.append(actor)

        if not actors:
            self._status(f"No {source} tracks to plot.")
            return

        def toggle(flag):
            for actor in actors:
                actor.SetVisibility(flag)

        self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(10, 220),
            size=25,
            color_on=color,
        )

    def setup_orientation_widgets(self) -> None:
        """Bottom-left E/N/Z axes (+Y = north in UTM)."""
        self.plotter.add_axes(**utm_orientation_axes_kwargs())

    def setup_z_scale(self) -> None:
        self.plotter.set_scale(1, 1, 2)

        def toggle_z_scale(flag):
            self.plotter.set_scale(1, 1, 2 if flag else 1)

        self.plotter.add_checkbox_button_widget(
            callback=toggle_z_scale,
            value=True,
            position=(10, 140),
            size=25,
            color_on=self.z_scale_toggle_color,
        )
