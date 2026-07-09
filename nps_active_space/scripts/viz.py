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
from nps_active_space.scripts.run_audible_transits import AudibleTransits
from nps_active_space.utils.models import Annotations
from nps_active_space.utils.computation import NMSIM_bbox_utm
import nps_active_space.utils.config as cfg
import argparse

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
    mic_color = "white"
    audible_annotation_color = "deepskyblue"
    inaudible_annotation_color = "red"
    audible_transits_color = "purple"
    z_scale_toggle_color = "black"  # button toggling z scale

    def __init__(self, unit, site, year, env, do_active=False, gain=None,
                 do_annots=False, do_transits=False,
                 annotation_file=None, audible_transits_pkl=None,
                 terraced=False, fill_layers=False, max_tracks=1000,
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

        # study area and crs
        self.study_area = load_studyarea(self.project_dir, self.unit, self.site, self.year)
        self.crs = NMSIM_bbox_utm(self.study_area)
        self.study_area = self.study_area.to_crs(self.crs)

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

        # configure plot parameters and display
        self.plotter.camera_position = "xz"
        self.plotter.camera.elevation = 45
        self.plotter.enable_terrain_style()
        self.setup_z_scale()
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

    def plot_activespace(self, terraced=False, gain=None):
        csv_3d_fits = f"{self.project_dir}/fits.csv"
        fit_results = pd.read_csv(csv_3d_fits, index_col="Designator")

        if gain is not None:
            print(f"Using gain {gain}dB")
        elif self.deployment in fit_results.index:
            gain = float(fit_results.loc[unit+site+year, "1/3rd Octave Gain (F1)"])
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

    def plot_contoured_activespace(self, active_3d):
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

    def plot_active_layer(self, active_layer, elevation, i=0):
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

    def plot_terraced_activespace(self, active_3d, layer_thickness=300, opacity=1):
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
    
    def plot_annotations(self, annotation_file=None):
        # load annotations
        if annotation_file is None:
            annotations = load_annotations(self.project_dir, self.unit, self.site, self.year, only_valid=True)
        else:
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

        self.plotter.add_checkbox_button_widget(
            callback=toggle,
            value=True,
            position=(10,180),
            size=25,
            color_on="purple"
        )

    def setup_z_scale(self):
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
               args.terraced, args.fill_layers, args.max_tracks)