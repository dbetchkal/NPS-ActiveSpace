"""Matplotlib tricontour helpers for active-space boundary geometry."""

import geopandas as gpd
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from shapely.geometry import Polygon, box
from shapely.validation import make_valid

from nps_active_space.utils.computation import build_src_point_mesh


def audibility_contours(
    total_space: gpd.GeoDataFrame,
    levels: list[float] | np.ndarray,
) -> mpl.tri.TriContourSet:
    """Delaunay-triangulate audibility points and return matplotlib contour set.

    Uses an Agg canvas instead of pyplot so this is safe inside Windows
    multiprocessing workers (TkAgg raises ``main thread is not in main loop``).
    """
    fig = Figure()
    FigureCanvasAgg(fig)
    ax = fig.add_subplot()
    tri = mpl.tri.Triangulation(
        total_space.geometry.x.tolist(),
        total_space.geometry.y.tolist(),
    )
    return ax.tricontour(tri, total_space.audible.tolist(), levels=levels)


def contour_active_space(
    total_space: gpd.GeoDataFrame,
    altitude_m: int,
) -> gpd.GeoDataFrame:
    """
    Use triangulation to select points along the audible/inaudible line of the active space to more precisely
    define the boundaries.

    Parameters
    ----------
    total_space : gpd.GeoDataFrame
        All points that have been tested for audibility so far -- both audible and inaudible.
    altitude_m : int
        The altitude (in meters) we are calculating the active space at.

    Returns
    -------
    points: gpd.GeoDataFrame
        GeoDataFrame of points to pass into the propagation model along the audible/inaudible line.
    """
    # Uses Delaunay triangulation
    cs = audibility_contours(total_space, levels=[0.5])  # contour with arbitrary point cloud

    contour_path = cs.get_paths()[0]
    x = contour_path.vertices[:, 0]
    y = contour_path.vertices[:, 1]
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x, y, altitude_m), crs=total_space.crs)

    return pts


def build_active_space(total_space: gpd.GeoDataFrame, crs: str) -> gpd.GeoDataFrame:
    """
    Build the final active space polygon given the audibility of all tested points.

    Parameters
    ----------
    total_space : gpd.GeoDataFrame
        A GeoDataFrame of all tested points and their audibility to be used to build the final active space
        polygon from.
    crs : str
        The crs of the points. Of the format 'epsg:XXXX...'

    Returns
    -------
    A GeoDataFrame of the active space.
    """
    # augment total space with inaudible points around the boundary, to avoid contour artifacts
    minx, miny, maxx, maxy = total_space.total_bounds
    region = gpd.GeoDataFrame({"geometry": [box(minx, miny, maxx, maxy)]}, crs=total_space.crs)
    region.geometry = region.buffer(100)  # small buffer so that the active space boundary occurs just outside the study area
    new_points = build_src_point_mesh(region)
    new_points = new_points[~new_points.geometry.within(region.union_all())]
    new_points["audible"] = 0
    total_space = pd.concat([total_space, new_points], ignore_index=True)

    # create the triangulated irregular network and contour lines
    levels = np.linspace(0, 1, 10, endpoint=False)
    cs = audibility_contours(total_space, levels=levels)

    # pick some contour level to choose as the active space boundary... 0.5 is somewhat arbitrary
    level_ind = np.where(cs.levels == 0.5)[0][0]

    # iterate through all contour paths in the `TriContourSet` at level_ind
    # https://matplotlib.org/stable/api/tri_api.html#matplotlib.tri.TriContourSet
    contour_path = cs.get_paths()[level_ind]  # in recent versions of `matplotlib` there is 1:1 correspondence `cs.levels` : `Path`
    polygons = [Polygon(P) for P in contour_path.to_polygons()]  # convert to `shapely.Polygon`

    # mark inner polygons as holes - otherwise when we dissolve or union_all(), they will disappear
    outer_polys = []
    for i, outer in enumerate(polygons):
        # skip if this is contained in another polygon
        if any(outer.within(other) for j, other in enumerate(polygons) if i != j):
            continue
        # find polygons inside this one
        inner_polys = [inner for j, inner in enumerate(polygons) if j != i and inner.within(outer)]
        holes = [inner.exterior.coords for inner in inner_polys]
        # construct polygon with holes
        poly = Polygon(shell=outer.exterior.coords, holes=holes)
        outer_polys.append(poly)

    # ensure valid geometries
    outer_polys = [make_valid(poly) if not poly.is_valid else poly for poly in outer_polys]

    return gpd.GeoDataFrame(data={"geometry": outer_polys}, geometry="geometry", crs=crs)
