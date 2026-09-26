from __future__ import annotations

import geopandas as gpd
import pyvista as pv
from shapely.geometry import Polygon

from nps_active_space.active_space.active_space_setup import resolve_3d_fit_gain
from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.utils.helpers import load_layered_activespace
from nps_active_space.viz.geometry import (
    active_to_linestrings,
    active_to_polys,
    create_polyline_3d,
    polygon_to_mesh,
)
from nps_active_space.viz.session import VizSession
from nps_active_space.viz.widgets import PlotWidgets


def model_display_name(model: AcousticModel) -> str:
    match AcousticModel.parse(model):
        case AcousticModel.AAM:
            return "AAM"
        case AcousticModel.NMSIM:
            return "NMSim"


class ActiveSpacePlotter:
    def __init__(self, session: VizSession, widgets: PlotWidgets) -> None:
        self._session = session
        self._widgets = widgets

    def resolve_gain(self, model: AcousticModel, gain: float | None) -> float | None:
        session = self._session
        if gain is not None:
            session.status(f"Using gain {gain} dB for {model}")
            return gain

        fitted = resolve_3d_fit_gain(
            session.project_dir, session.unit, session.site, session.year, model=model,
        )
        if fitted is not None:
            csv_3d_fits = p.fits_csv(session.project_dir)
            session.status(f"Using fitted gain {fitted} dB for {model} ({csv_3d_fits})")
            return fitted

        session.status(f"No fitted gain for {model}; pass -g or run fit_3d_active_space first.")
        return None

    def plot(
        self,
        terraced: bool = False,
        gain: float | None = None,
        model: AcousticModel = AcousticModel.NMSIM,
    ) -> None:
        session = self._session
        style = session.style
        gain = self.resolve_gain(model, gain)
        if gain is None:
            return

        color = (
            style.aam_activespace_color
            if AcousticModel.parse(model) is AcousticModel.AAM
            else style.nmsim_activespace_color
        )
        active_3d = load_layered_activespace(
            session.project_dir, session.unit, session.site, session.year,
            gain, session.crs, model=model,
        )
        if not active_3d.layer_dirs or active_3d.activespaces is None:
            site = p.site_dir(session.project_dir, session.unit, session.site)
            root = p.model_activespaces_dir(site, model)
            session.status(
                f"No active space layers loaded for {model_display_name(model)} "
                f"at gain {gain} dB (looked under {root})."
            )
            return
        prefix = model_display_name(model)
        session.legend_models.append((prefix, color))
        session.status(
            f"Loaded {prefix} at {gain} dB: "
            f"{', '.join(f'{z} m' for z in active_3d.activespaces)}"
        )
        if terraced:
            self.plot_terraced(active_3d, color=color, label_prefix=prefix)
        else:
            self.plot_contoured(active_3d, color=color, label_prefix=prefix)

    def plot_compare(self, terraced: bool = False, gain: float | None = None) -> None:
        session = self._session
        style = session.style
        if gain is not None:
            session.status(
                "Compare mode uses each model's fitted gain from fits.csv; ignoring -g."
            )
        widget_row = 0
        for model, color in (
            (AcousticModel.NMSIM, style.nmsim_activespace_color),
            (AcousticModel.AAM, style.aam_activespace_color),
        ):
            model_gain = self.resolve_gain(model, None)
            if model_gain is None:
                continue
            active_3d = load_layered_activespace(
                session.project_dir, session.unit, session.site, session.year,
                model_gain, session.crs, model=model,
            )
            if not active_3d.layer_dirs or active_3d.activespaces is None:
                site = p.site_dir(session.project_dir, session.unit, session.site)
                root = p.model_activespaces_dir(site, model)
                session.status(
                    f"No active space layers loaded for {model_display_name(model)} "
                    f"at gain {model_gain} dB (looked under {root})."
                )
                continue
            session.status(
                f"Loaded {model_display_name(model)} at {model_gain} dB: "
                f"{', '.join(f'{z} m' for z in active_3d.activespaces)}"
            )
            prefix = model_display_name(model)
            session.legend_models.append((prefix, color))
            if terraced:
                self.plot_terraced(
                    active_3d, color=color, label_prefix=prefix, widget_row=widget_row,
                )
                widget_row += 1
            else:
                widget_row = self.plot_contoured(
                    active_3d, color=color, label_prefix=prefix, widget_row=widget_row,
                )

    def plot_contoured(
        self, active_3d, color=None, label_prefix: str = "", widget_row: int = 0,
    ) -> int:
        session = self._session
        style = session.style
        if color is None:
            color = style.activespace_color
        if active_3d is None or active_3d.activespaces is None:
            return widget_row
        layer_checkboxes = []
        layer_callbacks = []
        row = widget_row
        prefix = f"{label_prefix} " if label_prefix else ""

        for active_z, active in active_3d.activespaces.items():
            if not active.empty:
                checkbox, toggle_cb = self.plot_layer(
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
            session.plotter.render()

        master_x = 10 + 160 * session.master_toggle_count
        session.master_toggle_count += 1
        self._widgets.add_labeled_checkbox(
            toggle_all_actives,
            value=True,
            position=(master_x, 5),
            size=35,
            color_on=color,
            label=f"{prefix}all".strip() or "all",
        )
        return row

    def plot_layer(
        self, active_layer, elevation: float, i: int = 0, color=None, label: str | None = None
    ):
        session = self._session
        style = session.style
        if color is None:
            color = style.activespace_color
        poly_actor = None
        if session.fill_layers:
            meshes = [
                polygon_to_mesh(poly, elevation) for poly in active_to_polys(active_layer)
            ]
            if meshes:
                poly_data = pv.PolyData().merge(meshes)
                poly_actor = session.plotter.add_mesh(poly_data, color=color, opacity=0.5)

        line_actors = []
        for line in active_to_linestrings(active_layer):
            polyline = create_polyline_3d(line, z=elevation)
            actor = session.plotter.add_mesh(
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

        checkbox = self._widgets.add_labeled_checkbox(
            toggle,
            value=True,
            position=self._widgets.layer_checkbox_xy(i),
            size=25,
            color_on=color,
            label=label or f"{int(elevation)} m",
        )

        return checkbox, toggle

    def plot_terraced(
        self,
        active_3d,
        layer_thickness: int = 300,
        opacity: float = 1,
        color=None,
        label_prefix: str = "",
        widget_row: int = 0,
    ) -> None:
        session = self._session
        style = session.style
        if color is None:
            color = style.activespace_color
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
                prev_active = layers[i - 1][1].to_crs(session.crs)
                sym_diff = active.union_all().symmetric_difference(prev_active.union_all())
                floor = gpd.GeoDataFrame(geometry=[sym_diff], crs=session.crs)

            for poly in active_to_polys(floor):
                mesh = polygon_to_mesh(poly, active_z - 0.5 * layer_thickness)
                meshes.append(mesh)

        stacked = pv.MultiBlock(meshes).combine()
        actor = session.plotter.add_mesh(stacked, color=color, opacity=opacity)

        def toggle(flag):
            if actor is not None:
                actor.SetVisibility(flag)

        prefix = f"{label_prefix} " if label_prefix else ""
        terrace_label = f"{prefix}terraced".strip()
        master_x = 10 + 160 * session.master_toggle_count
        session.master_toggle_count += 1
        self._widgets.add_labeled_checkbox(
            toggle,
            value=True,
            position=(master_x, 5),
            size=35,
            color_on=color,
            label=terrace_label,
        )
