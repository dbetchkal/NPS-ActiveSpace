import os
import numpy as np
import pandas as pd
import geopandas as gpd
import glob
import pyvista as pv
import rasterio
import pyproj
from shapely.geometry import box, Polygon, MultiPolygon, LineString, MultiLineString
from nps_active_space.utils.helpers import get_deployment, load_annotations, load_DEM, load_layered_activespace, load_studyarea
from nps_active_space.utils import paths as p
from nps_active_space.utils.enums import AcousticModel
from nps_active_space.active_space.active_space_setup import resolve_3d_fit_gain
from nps_active_space.scripts.run_audible_transits import AudibleTransits
from nps_active_space.utils.models import Annotations
from nps_active_space.utils.computation import study_area_utm_crs
import nps_active_space.utils.config as cfg
import argparse
from vtkmodules.vtkRenderingCore import vtkTextActor

# helper functions ============================================

def polygon_to_mesh(polygon, z):
    exterior = np.array(polygon.exterior.coords)
    points = np.c_[exterior[:, 0], exterior[:, 1], np.full(len(exterior), z)]
    poly = pv.PolyData(points)
    poly.faces = np.hstack([[len(points)], np.arange(len(points))])
    return poly.triangulate()


def active_to_polys(active):
    geometry = active.geometry.iloc[0]
    if isinstance(geometry, Polygon):
        polys = [geometry]
    elif isinstance(geometry, MultiPolygon):
        polys = geometry.geoms
    return polys


def active_to_linestrings(active):
    geometry = active.boundary.iloc[0]
    if isinstance(geometry, MultiLineString):
        linestrings = geometry.geoms
    elif isinstance(geometry, LineString):
        linestrings = [geometry]
    return linestrings


def create_polyline_3d(linestring, z=None):
    coords = np.array(linestring.coords)
    if z is not None:
        coords = np.column_stack((coords, np.full(coords.shape[0], z)))
    assert coords.shape[1] == 3
    coords[:,2] = np.clip(coords[:,2], 0, 10000)  # reduce magnitude of erroneous elevations
    n_points = coords.shape[0]
    lines = np.hstack([[n_points], np.arange(n_points)])
    poly = pv.PolyData(coords)
    poly.lines = lines
    return poly




# main class =====================================================

class Visualizer():
    # color config
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
    z_scale_toggle_color = "black"  # button toggling z scale

    def __init__(self, unit, site, year, env, do_active=False, gain=None,
                 do_annots=False, do_transits=False,
                 annotation_file=None, audible_transits_pkl=None,
                 terraced=False, fill_layers=False, max_tracks=1000,
                 model: AcousticModel | None = None, compare_models: bool = False,
                 ):
        # class metadata
        self.unit = unit
        self.site = site
        self.year = year
        self.deployment = f"{unit}{site}{year}"
        cfg.initialize(env)
        self.project_dir = cfg.read("project", "dir")
        self.fill_layers = fill_layers
        self.max_tracks = max_tracks
        self._legend_models: list[tuple[str, str]] = []
        self._master_toggle_count = 0

        # study area and crs
        self.study_area = load_studyarea(self.project_dir, self.unit, self.site, self.year)
        self.crs = study_area_utm_crs(self.study_area)
        self.study_area = self.study_area.to_crs(self.crs)

        # plot each element
        self.plotter = pv.Plotter()
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

        # configure plot parameters and display
        self.plotter.camera_position = "xz"
        self.plotter.camera.elevation = 45
        self.plotter.enable_terrain_style()
        self.setup_z_scale()
        self._add_color_legend(compare_models=compare_models)
        self.plotter.add_title(f"{unit}{site}{year}", font_size=12)
        self.plotter.show_axes()
        self.plotter.show()
    
    def plot_dem(self, show_scalar_bar=False):
        # load DEM
        dem = load_DEM(self.project_dir, self.unit, self.site)    
        data = dem.read(1)
        if dem.nodata is not None:
            data[data == dem.nodata] = 0
        data[data > 9000] = 0  # higher than any elevation on earth, should be nodata
        
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

    def plot_point(self, x, y, z, color="white"):
        point = pv.PolyData(np.array([[x, y, z]]))
        self.plotter.add_mesh(point, color=color, point_size=10, render_points_as_spheres=True)

    def plot_mic(self):
        mic = get_deployment(self.project_dir, self.unit, self.site, self.year)
        mic = mic.to_crs(self.crs)
        self.plot_point(mic.x, mic.y, mic.z, self.mic_color)

    def _site_dir(self) -> str:
        return p.site_dir(self.project_dir, self.unit, self.site)

    def _resolve_activespace_gain(self, model: AcousticModel, gain: float | None) -> float | None:
        if gain is not None:
            print(f"Using gain {gain} dB for {model}")
            return gain

        fitted = resolve_3d_fit_gain(
            self.project_dir, self.unit, self.site, self.year, model=model,
        )
        if fitted is not None:
            csv_3d_fits = p.fits_csv(self.project_dir)
            print(f"Using fitted gain {fitted} dB for {model} ({csv_3d_fits})")
            return fitted

        print(f"No fitted gain for {model}; pass -g or run fit_3d_active_space first.")
        return None

    def plot_activespace(self, terraced=False, gain=None, model: AcousticModel = AcousticModel.NMSIM):
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
            print(
                f"No active space layers loaded for {self._model_display_name(model)} "
                f"at gain {gain} dB (looked under {root})."
            )
            return
        prefix = self._model_display_name(model)
        self._legend_models.append((prefix, color))
        print(
            f"Loaded {prefix} at {gain} dB: "
            f"{', '.join(f'{z} m' for z in active_3d.activespaces)}"
        )
        if terraced:
            self.plot_terraced_activespace(active_3d, color=color, label_prefix=prefix)
        else:
            self.plot_contoured_activespace(active_3d, color=color, label_prefix=prefix)

    def plot_compare_activespaces(self, terraced=False, gain=None):
        widget_row = 0
        for model, color in (
            (AcousticModel.NMSIM, self.nmsim_activespace_color),
            (AcousticModel.AAM, self.aam_activespace_color),
        ):
            model_gain = gain if gain is not None else self._resolve_activespace_gain(model, None)
            if model_gain is None:
                continue
            active_3d = load_layered_activespace(
                self.project_dir, self.unit, self.site, self.year,
                model_gain, self.crs, model=model,
            )
            if not active_3d.layer_dirs or active_3d.activespaces is None:
                site = p.site_dir(self.project_dir, self.unit, self.site)
                root = p.model_activespaces_dir(site, model)
                print(
                    f"No active space layers loaded for {self._model_display_name(model)} "
                    f"at gain {model_gain} dB (looked under {root})."
                )
                continue
            print(
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
        # Vertically center the label on the square.
        text.SetPosition(x + size + 8, y + max(2, (size - 16) // 2))
        prop = text.GetTextProperty()
        prop.SetFontFamilyToArial()
        prop.SetFontSize(16)
        prop.SetColor(1.0, 1.0, 1.0)
        prop.BoldOn()
        prop.ShadowOn()
        prop.SetJustificationToLeft()
        prop.SetVerticalJustificationToBottom()
        self.plotter.renderer.AddActor2D(text)
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

    def plot_active_layer(self, active_layer, elevation, i=0, color=None, label: str | None = None):
        if color is None:
            color = self.activespace_color
        poly_actor = None
        if self.fill_layers:
            meshes = []
            for poly in active_to_polys(active_layer):
                meshes.append(polygon_to_mesh(poly, elevation))
                poly_data = pv.PolyData().merge(meshes)
            poly_actor = self.plotter.add_mesh(poly_data, color=color, opacity=0.5)

        line_actors = []
        for line in active_to_linestrings(active_layer):
            polyline = create_polyline_3d(line, z=elevation)
            actor = self.plotter.add_mesh(polyline, color=color, point_size=2, line_width=2)
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
        self, active_3d, layer_thickness=300, opacity=1, color=None,
        label_prefix: str = "", widget_row: int = 0,
    ):
        if color is None:
            color = self.activespace_color
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
    
    def plot_annotations(self, annotation_file=None):
        # load annotations
        if annotation_file is None:
            annotations = load_annotations(self.project_dir, self.unit, self.site, self.year, only_valid=True)
        else:
            if not os.path.isabs(annotation_file):
                site_path = os.path.join(self.project_dir, f"{self.unit}{self.site}", annotation_file)
                if os.path.isfile(site_path):
                    annotation_file = site_path
            print(f"Loading annotations from custom file: {annotation_file}")
            annotations = Annotations(annotation_file, only_valid=True)
        
        if annotations.empty:
            print("No annotations found, skipping.")
            return
        
        # downsample n tracks if too many
        track_ids = annotations["_id"].drop_duplicates()
        print(f"{len(track_ids)} tracks = {len(annotations)} annotated segments")
        if len(track_ids) > self.max_tracks:
            print(f"More than {self.max_tracks}, sampling")
            selected_track_ids = track_ids.sample(self.max_tracks, replace=False, random_state=2)
            annotations = annotations[annotations["_id"].isin(selected_track_ids)]
            print(f"Sampled tracks, showing {len(selected_track_ids)} tracks = {len(annotations)} annotated segments")

        # clip to the area - using a box() plays more nicely with z values than the study area itself
        annotations = annotations.to_crs(self.crs)
        annotations = annotations.clip(box(*self.study_area.total_bounds)).explode()

        # create and plot segments
        audible_actors = []
        inaudible_actors = []
        for _, annot in annotations.iterrows():
            color = "deepskyblue" if annot["audible"] else "red"
            polyline = create_polyline_3d(annot["geometry"])
            actor = self.plotter.add_mesh(polyline, color=color, point_size=2, line_width=2)
            if annot["audible"]:
                audible_actors.append(actor)
            else:
                inaudible_actors.append(actor)
        
        def toggle_audible(flag):
            for actor in audible_actors:
                actor.SetVisibility(flag)

        def toggle_inaudible(flag):
            for actor in inaudible_actors:
                actor.SetVisibility(flag)
        
        self._add_labeled_checkbox(
            toggle_audible,
            value=True,
            position=(10, 100),
            size=25,
            color_on="deepskyblue",
            label="audible",
        )
        self._add_labeled_checkbox(
            toggle_inaudible,
            value=True,
            position=(10, 60),
            size=25,
            color_on="red",
            label="inaudible",
        )
     
    def plot_audible_transits(self, audible_transits_pkl=None):
        if audible_transits_pkl:
            print(f"Loading audible transits from {audible_transits_pkl}")
            listener = AudibleTransits.from_pickle(audible_transits_pkl)
        else:
            print("Loading audible transits")
            matches = glob.glob(os.path.join(
                self.project_dir, self.unit+self.site, "Output_Data", "AUDIBLE_TRANSITS",
                f"3D*{year}-01-01*Active Space {self.year}*", "AudibleTransits_object.pkl"))
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

        self._add_labeled_checkbox(
            toggle,
            value=True,
            position=(10, 180),
            size=25,
            color_on="purple",
            label="transits",
        )

    def setup_z_scale(self):
        self.plotter.set_scale(1, 1, 2)  # easier to make out elevation differences

        # allow toggling though
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # required args
    # we pass deployment as a single positional arg instead of -u DENA -s TRLA -y 2024,
    # because it saves a lot of keystrokes, and this script is often used frequently
    parser.add_argument("deployment", help="Deployment name, e.g. DENATRLA2024")
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
    parser.add_argument("--all", action="store_true",
                        help="Load everything, shorthand for --active-space --annotations --audible-transits")

    # uncommon / special use case args
    parser.add_argument("-m", "--max-tracks", default=500,
                        help="Maximum number of annotation tracks or audible transits to show.")
    parser.add_argument("--annotation-file", help="Path to .geojson file to load annotations from, if not the default.")
    parser.add_argument("--transits-pkl", help="Path to .pkl file to load audible transits from, if not the default.")
    parser.add_argument("--terraced", action="store_true",
                        help="If included, render the active space as the terraced surface instead of contours.")
    parser.add_argument("--model", type=AcousticModel, choices=list(AcousticModel),
                        default=None,
                        help="Propagation model for active-space layers (default: nmsim).")
    parser.add_argument("--compare", action="store_true",
                        help="Overlay NMSim (orange) and AAM (cyan) active spaces.")
    parser.add_argument("--fill-layers", action="store_true",
                        help="If included, fill the interior of each active space contour polygon.")

    args = parser.parse_args()

    usy = args.deployment
    unit, site, year = usy[:4], usy[4:-4], usy[-4:]
    print(unit, site, year)

    do_active = args.active_space or args.all
    do_annotations = args.annotations or args.all
    do_transits = args.audible_transits or args.all

    Visualizer(unit, site, year, args.environment,
               do_active, args.gain,
               do_annotations, do_transits,
               args.annotation_file, args.transits_pkl,
               args.terraced, args.fill_layers, args.max_tracks,
               model=args.model, compare_models=args.compare)