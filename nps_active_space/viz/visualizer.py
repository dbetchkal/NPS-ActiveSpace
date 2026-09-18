from __future__ import annotations

import pyproj
import pyvista as pv
from shapely.geometry import LineString
from vtkmodules.vtkCommonCore import vtkObject

import nps_active_space.utils.config as cfg
from nps_active_space.utils.computation import study_area_utm_crs
from nps_active_space.utils.enums import AcousticModel, TrackSource
from nps_active_space.utils.helpers import get_logger, load_studyarea
from nps_active_space.viz.active_space_plot import ActiveSpacePlotter, model_display_name
from nps_active_space.viz.annotations_plot import AnnotationsPlotter
from nps_active_space.viz.markers import WINDOW_TITLE
from nps_active_space.viz.polyline_builders import TrackPolylineBuilder
from nps_active_space.viz.scene import ScenePlotter
from nps_active_space.viz.session import VizSession
from nps_active_space.viz.style import VizStyle
from nps_active_space.viz.tracks_plot import TracksPlotter
from nps_active_space.viz.widgets import PlotWidgets

_DEFAULT_STYLE = VizStyle()


class Visualizer:
    """PyVista deployment viewer; composes focused plotters over a shared VizSession."""

    activespace_color = _DEFAULT_STYLE.activespace_color
    nmsim_activespace_color = _DEFAULT_STYLE.nmsim_activespace_color
    aam_activespace_color = _DEFAULT_STYLE.aam_activespace_color
    mic_color = _DEFAULT_STYLE.mic_color
    _layer_widget_dy = _DEFAULT_STYLE.layer_widget_dy
    _layer_checkbox_x = _DEFAULT_STYLE.layer_checkbox_x
    _layer_checkbox_y0 = _DEFAULT_STYLE.layer_checkbox_y0
    audible_annotation_color = _DEFAULT_STYLE.audible_annotation_color
    inaudible_annotation_color = _DEFAULT_STYLE.inaudible_annotation_color
    audible_transits_color = _DEFAULT_STYLE.audible_transits_color
    vessel_track_color = _DEFAULT_STYLE.vessel_track_color
    flight_track_color = _DEFAULT_STYLE.flight_track_color
    z_scale_toggle_color = _DEFAULT_STYLE.z_scale_toggle_color
    sea_surface_offset_m = _DEFAULT_STYLE.sea_surface_offset_m
    sea_surface_densify_step_m = _DEFAULT_STYLE.sea_surface_densify_step_m

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
        project_dir = cfg.read("project", "dir")
        self.fill_layers = fill_layers
        self.max_tracks = max_tracks
        logger = get_logger("VIZ", verbose=True)

        study_area = load_studyarea(project_dir, unit, site, year)
        crs = study_area_utm_crs(study_area)
        study_area = study_area.to_crs(crs)

        self._session = VizSession(
            unit=unit,
            site=site,
            year=year,
            deployment=self.deployment,
            project_dir=project_dir,
            crs=crs,
            study_area=study_area,
            plotter=pv.Plotter(title=WINDOW_TITLE),
            logger=logger,
            style=_DEFAULT_STYLE,
            fill_layers=fill_layers,
            max_tracks=max_tracks,
            to_wgs84=pyproj.Transformer.from_crs(crs, "epsg:4326", always_xy=True),
        )
        self.project_dir = self._session.project_dir
        self.study_area = self._session.study_area
        self.crs = self._session.crs
        self.logger = self._session.logger
        self.plotter = self._session.plotter

        self._widgets = PlotWidgets(self._session)
        self._scene = ScenePlotter(self._session)
        self._polylines = TrackPolylineBuilder(self._session)
        self._active_space = ActiveSpacePlotter(self._session, self._widgets)
        self._annotations = AnnotationsPlotter(
            self._session, self._widgets, self._scene, self._polylines,
        )
        self._tracks = TracksPlotter(self._session, self._widgets, self._scene, self._polylines)

        self._scene.plot_dem()
        self._scene.plot_mic()
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
        self._widgets.setup_z_scale()
        self._widgets.add_color_legend(compare_models=compare_models)
        self.plotter.add_title(f"{unit}{site}{year}", font_size=12)
        self._scene.setup_orientation_widgets()
        self.plotter.reset_camera()
        self.plotter.camera.elevation = 30
        vtkObject.GlobalWarningDisplayOff()
        try:
            self.plotter.show()
        finally:
            vtkObject.GlobalWarningDisplayOn()

    def _status(self, message: str) -> None:
        self._session.status(message)

    def _add_track_line(self, polyline, *, color: str, line_width: int = 2):
        return self._scene.add_track_line(polyline, color=color, line_width=line_width)

    def _add_annotation_lines(self, polylines, *, color: str, line_width: int = 2):
        return self._scene.add_annotation_lines(polylines, color=color, line_width=line_width)

    @staticmethod
    def _flat_sea_surface_polyline(linestring: LineString, offset_m: float):
        return TrackPolylineBuilder.flat_sea_surface_polyline(linestring, offset_m)

    def plot_dem(self, show_scalar_bar: bool = False) -> None:
        self._scene.plot_dem(show_scalar_bar=show_scalar_bar)

    def plot_point(self, x: float, y: float, z: float, color: str = "white") -> None:
        self._scene.plot_point(x, y, z, color=color)

    def plot_mic(self) -> None:
        self._scene.plot_mic()

    def _resolve_activespace_gain(self, model: AcousticModel, gain: float | None) -> float | None:
        return self._active_space.resolve_gain(model, gain)

    def plot_activespace(
        self,
        terraced: bool = False,
        gain: float | None = None,
        model: AcousticModel = AcousticModel.NMSIM,
    ) -> None:
        self._active_space.plot(terraced, gain, model=model)

    def plot_compare_activespaces(
        self, terraced: bool = False, gain: float | None = None
    ) -> None:
        self._active_space.plot_compare(terraced, gain)

    @staticmethod
    def _model_display_name(model: AcousticModel) -> str:
        return model_display_name(model)

    def _layer_checkbox_xy(self, row: int) -> tuple[int, int]:
        return self._widgets.layer_checkbox_xy(row)

    def _add_labeled_checkbox(self, callback, **kwargs):
        return self._widgets.add_labeled_checkbox(callback, **kwargs)

    def _add_color_legend(self, *, compare_models: bool) -> None:
        self._widgets.add_color_legend(compare_models=compare_models)

    def plot_contoured_activespace(self, active_3d, color=None, label_prefix: str = "", widget_row: int = 0) -> int:
        return self._active_space.plot_contoured(active_3d, color, label_prefix, widget_row)

    def plot_active_layer(self, active_layer, elevation: float, i: int = 0, color=None, label: str | None = None):
        return self._active_space.plot_layer(active_layer, elevation, i, color, label)

    def plot_terraced_activespace(self, active_3d, layer_thickness: int = 300, opacity: float = 1, color=None, label_prefix: str = "", widget_row: int = 0) -> None:
        self._active_space.plot_terraced(active_3d, layer_thickness, opacity, color, label_prefix, widget_row)

    def _annotation_polyline(self, linestring: LineString):
        return self._polylines.annotation_polyline(linestring)

    def _sea_surface_polyline(self, linestring: LineString):
        return self._polylines.sea_surface_polyline(linestring)

    def plot_annotations(self, annotation_file: str | None = None) -> None:
        self._annotations.plot(annotation_file)

    def plot_audible_transits(self, audible_transits_pkl: str | None = None) -> None:
        self._tracks.plot_audible_transits(audible_transits_pkl)

    def plot_tracks(
        self,
        source: TrackSource,
        start_date: str | None = None,
        end_date: str | None = None,
    ) -> None:
        self._tracks.plot_tracks(
            source,
            start_date,
            end_date,
            annotation_polyline=self._annotation_polyline,
            sea_surface_polyline=self._sea_surface_polyline,
            add_track_line=self._add_track_line,
        )

    def setup_orientation_widgets(self) -> None:
        self._scene.setup_orientation_widgets()

    def setup_z_scale(self) -> None:
        self._widgets.setup_z_scale()

    @property
    def dem(self):
        return self._session.dem

    @dem.setter
    def dem(self, value) -> None:
        self._session.dem = value

    @property
    def _dem_sampler(self):
        return self._session.dem_sampler

    @_dem_sampler.setter
    def _dem_sampler(self, value) -> None:
        self._session.dem_sampler = value

    @property
    def _legend_models(self) -> list[tuple[str, str]]:
        return self._session.legend_models

    @_legend_models.setter
    def _legend_models(self, value: list[tuple[str, str]]) -> None:
        self._session.legend_models = value

    @property
    def _master_toggle_count(self) -> int:
        return self._session.master_toggle_count

    @_master_toggle_count.setter
    def _master_toggle_count(self, value: int) -> None:
        self._session.master_toggle_count = value
