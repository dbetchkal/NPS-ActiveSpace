from __future__ import annotations

import glob
import configparser
import os

import geopandas as gpd
import numpy as np
import pyproj
import pyvista as pv
import rasterio
from shapely.geometry import LineString, Polygon, box
from tqdm import tqdm
from vtkmodules.vtkCommonCore import vtkObject
from vtkmodules.vtkRenderingCore import vtkTextActor

import nps_active_space.utils.config as cfg
from nps_active_space.active_space.active_space_setup import resolve_3d_fit_gain
from nps_active_space.ground_truthing.load_tracks import load_tracks
from nps_active_space.scripts.run_audible_transits import AudibleTransits
from nps_active_space.utils import paths as p
from nps_active_space.utils.computation import study_area_utm_crs
from nps_active_space.utils.enums import AcousticModel, TrackSource
from nps_active_space.utils.helpers import (
    get_deployment,
    get_logger,
    load_annotations,
    load_DEM,
    load_layered_activespace,
    load_studyarea,
)
from nps_active_space.utils.models import Annotations
from nps_active_space.viz.annotations import format_annotation_summary
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
from nps_active_space.viz.markers import (
    WINDOW_TITLE,
    utm_orientation_axes_kwargs,
)


class Visualizer:
    activespace_color = "orange"
    nmsim_activespace_color = "orange"
    aam_activespace_color = "cyan"
    mic_color = "white"
    _layer_widget_dy = 28
    _layer_checkbox_x = 10
    _layer_checkbox_y0 = 220
    audible_annotation_color = "deepskyblue"
    inaudible_annotation_color = "red"
    audible_transits_color = "purple"
    vessel_track_color = "magenta"
    flight_track_color = "indigo"
    z_scale_toggle_color = "black"
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
        model: AcousticModel | None = None,
        compare_models: bool = False,
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
        self._legend_models: list[tuple[str, str]] = []
        self._master_toggle_count = 0

        self.study_area = load_studyarea(self.project_dir, self.unit, self.site, self.year)
        self.crs = study_area_utm_crs(self.study_area)
        self.study_area = self.study_area.to_crs(self.crs)
        self._to_wgs84 = pyproj.Transformer.from_crs(self.crs, "epsg:4326", always_xy=True)

        self.plotter = pv.Plotter(title=WINDOW_TITLE)
        self.plot_dem()
        self.plot_mic()
        if do_active:
            if compare_models:
                self.plot_compare_activespaces(terraced, gain)
            else:
                self.plot_activespace(terraced, gain, model=model or AcousticModel.NMSIM)
        if do_annots:
            self.plot_annotations(annotation_file)
        if do_transits:
            self.plot_audible_transits(audible_transits_pkl)
        if track_source is not None:
            self.plot_tracks(track_source, track_start_date, track_end_date)

        self.plotter.enable_terrain_style()
        self.setup_z_scale()
        self._add_color_legend(compare_models=compare_models)
        self.plotter.add_title(f"{unit}{site}{year}", font_size=12)
        self.setup_orientation_widgets()
        self.plotter.reset_camera()
        self.plotter.camera.elevation = 30
        vtkObject.GlobalWarningDisplayOff()
        try:
            self.plotter.show()
        finally:
            vtkObject.GlobalWarningDisplayOn()

    def _status(self, message: str) -> None:
        """Log a user-facing status line to the console (via get_logger StreamHandler)."""
        self.logger.info(message)

    def _add_track_line(
        self, polyline: pv.PolyData, *, color: str, line_width: int = 2
    ):
        """Add a causal track polyline over the DEM (plain lines, not tubes)."""
        return self.plotter.add_mesh(
            polyline,
            color=color,
            line_width=line_width,
            render_lines_as_tubes=False,
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

    def _site_dir(self) -> str:
        return p.site_dir(self.project_dir, self.unit, self.site)

    def _resolve_activespace_gain(self, model: AcousticModel, gain: float | None) -> float | None:
        if gain is not None:
            self._status(f"Using gain {gain} dB for {model}")
            return gain

        fitted = resolve_3d_fit_gain(
            self.project_dir, self.unit, self.site, self.year, model=model,
        )
        if fitted is not None:
            csv_3d_fits = p.fits_csv(self.project_dir)
            self._status(f"Using fitted gain {fitted} dB for {model} ({csv_3d_fits})")
            return fitted

        self._status(f"No fitted gain for {model}; pass -g or run fit_3d_active_space first.")
        return None

    def plot_activespace(
        self,
        terraced: bool = False,
        gain: float | None = None,
        model: AcousticModel = AcousticModel.NMSIM,
    ) -> None:
        gain = self._resolve_activespace_gain(model, gain)
        if gain is None:
            return

        color = (
            self.aam_activespace_color
            if AcousticModel.parse(model) is AcousticModel.AAM
            else self.nmsim_activespace_color
        )
        active_3d = load_layered_activespace(
            self.project_dir, self.unit, self.site, self.year,
            gain, self.crs, model=model,
        )
        if not active_3d.layer_dirs or active_3d.activespaces is None:
            site = p.site_dir(self.project_dir, self.unit, self.site)
            root = p.model_activespaces_dir(site, model)
            self._status(
                f"No active space layers loaded for {self._model_display_name(model)} "
                f"at gain {gain} dB (looked under {root})."
            )
            return
        prefix = self._model_display_name(model)
        self._legend_models.append((prefix, color))
        self._status(
            f"Loaded {prefix} at {gain} dB: "
            f"{', '.join(f'{z} m' for z in active_3d.activespaces)}"
        )
        if terraced:
            self.plot_terraced_activespace(active_3d, color=color, label_prefix=prefix)
        else:
            self.plot_contoured_activespace(active_3d, color=color, label_prefix=prefix)

    def plot_compare_activespaces(
        self, terraced: bool = False, gain: float | None = None
    ) -> None:
        if gain is not None:
            self._status(
                "Compare mode uses each model's fitted gain from fits.csv; ignoring -g."
            )
        widget_row = 0
        for model, color in (
            (AcousticModel.NMSIM, self.nmsim_activespace_color),
            (AcousticModel.AAM, self.aam_activespace_color),
        ):
            model_gain = self._resolve_activespace_gain(model, None)
            if model_gain is None:
                continue
            active_3d = load_layered_activespace(
                self.project_dir, self.unit, self.site, self.year,
                model_gain, self.crs, model=model,
            )
            if not active_3d.layer_dirs or active_3d.activespaces is None:
                site = p.site_dir(self.project_dir, self.unit, self.site)
                root = p.model_activespaces_dir(site, model)
                self._status(
                    f"No active space layers loaded for {self._model_display_name(model)} "
                    f"at gain {model_gain} dB (looked under {root})."
                )
                continue
            self._status(
                f"Loaded {self._model_display_name(model)} at {model_gain} dB: "
                f"{', '.join(f'{z} m' for z in active_3d.activespaces)}"
            )
            prefix = self._model_display_name(model)
            self._legend_models.append((prefix, color))
            if terraced:
                self.plot_terraced_activespace(
                    active_3d, color=color, label_prefix=prefix, widget_row=widget_row,
                )
                widget_row += 1
            else:
                widget_row = self.plot_contoured_activespace(
                    active_3d, color=color, label_prefix=prefix, widget_row=widget_row,
                )

    @staticmethod
    def _model_display_name(model: AcousticModel) -> str:
        match AcousticModel.parse(model):
            case AcousticModel.AAM:
                return "AAM"
            case AcousticModel.NMSIM:
                return "NMSim"

    def _layer_checkbox_xy(self, row: int) -> tuple[int, int]:
        return (
            self._layer_checkbox_x,
            self._layer_checkbox_y0 + self._layer_widget_dy * row,
        )

    def _add_labeled_checkbox(
        self,
        callback,
        *,
        value: bool,
        position: tuple[int, int],
        size: int,
        color_on: str,
        label: str,
    ):
        """Checkbox plus a 2D label in the same VTK display-pixel space."""
        checkbox = self.plotter.add_checkbox_button_widget(
            callback=callback,
            value=value,
            position=position,
            size=size,
            color_on=color_on,
        )
        x, y = position
        text = vtkTextActor()
        text.SetInput(label)
        text.SetTextScaleModeToNone()
        text.GetPositionCoordinate().SetCoordinateSystemToDisplay()
        text.SetPosition(x + size + 8, y + max(2, (size - 16) // 2))
        prop = text.GetTextProperty()
        prop.SetFontFamilyToArial()
        prop.SetFontSize(16)
        prop.SetColor(1.0, 1.0, 1.0)
        prop.BoldOn()
        prop.ShadowOn()
        prop.SetJustificationToLeft()
        prop.SetVerticalJustificationToBottom()
        self.plotter.renderer.AddViewProp(text)
        return checkbox

    def _add_color_legend(self, *, compare_models: bool) -> None:
        legend_entries = self._legend_models
        if not legend_entries and compare_models:
            legend_entries = [
                ("NMSim", self.nmsim_activespace_color),
                ("AAM", self.aam_activespace_color),
            ]
        if not legend_entries:
            return
        unique: list[tuple[str, str]] = []
        seen: set[str] = set()
        for name, color in legend_entries:
            if name not in seen:
                unique.append((name, color))
                seen.add(name)
        self.plotter.add_legend(
            unique,
            loc="upper right",
            bcolor=None,
            face="r",
        )

    def plot_contoured_activespace(
        self, active_3d, color=None, label_prefix: str = "", widget_row: int = 0,
    ) -> int:
        if color is None:
            color = self.activespace_color
        if active_3d is None or active_3d.activespaces is None:
            return widget_row
        layer_checkboxes = []
        layer_callbacks = []
        row = widget_row
        prefix = f"{label_prefix} " if label_prefix else ""

        for active_z, active in active_3d.activespaces.items():
            if not active.empty:
                checkbox, toggle_cb = self.plot_active_layer(
                    active,
                    active_z,
                    i=row,
                    color=color,
                    label=f"{prefix}{int(active_z)} m".strip(),
                )
                layer_checkboxes.append(checkbox)
                layer_callbacks.append(toggle_cb)
                row += 1

        def toggle_all_actives(flag):
            for box, toggle_cb in zip(layer_checkboxes, layer_callbacks):
                box.GetRepresentation().SetState(int(flag))
                toggle_cb(flag)
            self.plotter.render()

        master_x = 10 + 160 * self._master_toggle_count
        self._master_toggle_count += 1
        self._add_labeled_checkbox(
            toggle_all_actives,
            value=True,
            position=(master_x, 5),
            size=35,
            color_on=color,
            label=f"{prefix}all".strip() or "all",
        )
        return row

    def plot_active_layer(
        self, active_layer, elevation: float, i: int = 0, color=None, label: str | None = None
    ):
        if color is None:
            color = self.activespace_color
        poly_actor = None
        if self.fill_layers:
            meshes = [
                polygon_to_mesh(poly, elevation) for poly in active_to_polys(active_layer)
            ]
            if meshes:
                poly_data = pv.PolyData().merge(meshes)
                poly_actor = self.plotter.add_mesh(poly_data, color=color, opacity=0.5)

        line_actors = []
        for line in active_to_linestrings(active_layer):
            polyline = create_polyline_3d(line, z=elevation)
            actor = self.plotter.add_mesh(
                polyline,
                color=color,
                point_size=2,
                line_width=2,
            )
            line_actors.append(actor)

        def toggle(flag):
            if poly_actor is not None:
                poly_actor.SetVisibility(flag)
            for actor in line_actors:
                actor.SetVisibility(flag)

        checkbox = self._add_labeled_checkbox(
            toggle,
            value=True,
            position=self._layer_checkbox_xy(i),
            size=25,
            color_on=color,
            label=label or f"{int(elevation)} m",
        )

        return checkbox, toggle

    def plot_terraced_activespace(
        self,
        active_3d,
        layer_thickness: int = 300,
        opacity: float = 1,
        color=None,
        label_prefix: str = "",
        widget_row: int = 0,
    ) -> None:
        if color is None:
            color = self.activespace_color
        meshes = []
        layers = list(active_3d.activespaces.items())
        for i in range(len(layers)):
            active_z, active = layers[i]
            if active.empty:
                continue

            for poly in active_to_polys(active):
                hole_polys = [Polygon(hole) for hole in poly.interiors]
                for poly_part in [poly] + hole_polys:
                    mesh = polygon_to_mesh(poly_part, active_z - 0.5 * layer_thickness)
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
        actor = self.plotter.add_mesh(stacked, color=color, opacity=opacity)

        def toggle(flag):
            if actor is not None:
                actor.SetVisibility(flag)

        prefix = f"{label_prefix} " if label_prefix else ""
        terrace_label = f"{prefix}terraced".strip()
        master_x = 10 + 160 * self._master_toggle_count
        self._master_toggle_count += 1
        self._add_labeled_checkbox(
            toggle,
            value=True,
            position=(master_x, 5),
            size=35,
            color_on=color,
            label=terrace_label,
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

        self._add_labeled_checkbox(
            toggle_audible,
            value=True,
            position=(10, 90),
            size=25,
            color_on="deepskyblue",
            label="audible",
        )
        self._add_labeled_checkbox(
            toggle_inaudible,
            value=True,
            position=(10, 55),
            size=25,
            color_on="red",
            label="inaudible",
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

        self._add_labeled_checkbox(
            toggle,
            value=True,
            position=(10, 180),
            size=25,
            color_on="purple",
            label="transits",
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

        self._add_labeled_checkbox(
            toggle,
            value=True,
            position=(10, 20),
            size=25,
            color_on=color,
            label="tracks",
        )

    def setup_orientation_widgets(self) -> None:
        """Bottom-left E/N/Z axes (+Y = north in UTM)."""
        self.plotter.add_axes(**utm_orientation_axes_kwargs())

    def setup_z_scale(self) -> None:
        self.plotter.set_scale(1, 1, 2)

        def toggle_z_scale(flag):
            self.plotter.set_scale(1, 1, 2 if flag else 1)

        self._add_labeled_checkbox(
            toggle_z_scale,
            value=True,
            position=(10, 140),
            size=25,
            color_on=self.z_scale_toggle_color,
            label="2× z-scale",
        )
