import datetime as dt
from typing import Iterable, List, Optional, TYPE_CHECKING

import geopandas as gpd
import numpy as np
from scipy import interpolate
from shapely.geometry import Point

if TYPE_CHECKING:
    from nps_active_space.utils.models import Tracks
    from shapely.geometry import Polygon, Point


__all__ = [
    'audible_time_delay',
    'build_src_point_mesh',
    'climb_angle',
    'coords_to_utm',
    'interpolate_spline',
    'NMSIM_bbox_utm'
]


def NMSIM_bbox_utm(study_area: gpd.GeoDataFrame) -> str:
    """
    NMSIM references an entire project to the westernmost extent of the elevation (or landcover) file.
    Given that, return the UTM Zone the project will eventually use. NMSIM uses NAD83 as its geographic
    coordinate system, so the study area will be projected into NAD83 before calculating the UTM zone.

    Parameters
    ----------
    study_area : gpd.GeoDataFrame
        A study area (Polygon) to find the UTM zone of the westernmost extent for.

    Returns
    -------
    UTM zone projection name (e.g.  'epsg:26905' for UTM 5N) that aligns with the westernmost extent of a study area.
    """
    if study_area.crs.to_epsg() != 4269:
        study_area = study_area.to_crs(epsg='4269')
    study_area_bbox = study_area.geometry.iloc[0].bounds  # (minx, miny, maxx, maxy)
    lat = study_area_bbox[3]  # maxy
    lon = study_area_bbox[0]  # minx

    return coords_to_utm(lat, lon)


def coords_to_utm(lat: float, lon: float) -> str:
    """
    Takes the latitude and longitude of a point and outputs the EPSG code corresponding to the UTM zone of the point.

    Parameters
    ----------
    lat : float
        Latitude of a point in decimal degrees in a WGS84 projection (EPSG:4326)
    lon : float
        Longitude of a point in decimal degrees in a WGS84 projection (EPSG:4326)

    Returns
    -------
    utm_proj : str
        UTM zone projection name (e.g.  'epsg:26905' for UTM 5N)

    Notes
    -----
    Remember: x=longitude, y=latitude
    """
    # 6 degrees per zone; add 180 because zone 1 starts at 180 W.
    utm_zone = int((lon + 180) // 6 + 1)

    # 269 = northern hemisphere, 327 = southern hemisphere
    utm_proj = 'epsg:269{:02d}'.format(utm_zone) if lat > 0 else 'epsg:327{:02d}'.format(utm_zone)
    return utm_proj


def climb_angle(v: Iterable) -> np.ndarray:
    """
    Compute the 'climb angle' of a vector.
    A = 𝑛•𝑏=|𝑛||𝑏|𝑠𝑖𝑛(𝜃)

    Parameters
    ----------
    v : array-like
        Vector to compute the climb angle for.

    Returns
    -------
    degrees : ndarray of floats
        Corresponding climb angle value in degrees.
    """
    n = np.array([0, 0, 1])  # A unit normal vector perpendicular to the xy plane
    degrees = np.degrees(np.arcsin(np.dot(n, v) / np.linalg.norm(n) * np.linalg.norm(v)))
    return degrees


def interpolate_spline(points: 'Tracks', ds: int = 1) -> gpd.GeoDataFrame:
    """
    Interpolate points with a cubic spline between flight points, if possible.
    See https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html#spline-interpolation for docs

    Parameters
    ----------
    points : Tracks
        A Track gpd.GeoDataframe object containing known track points in a path. A minimum of 2 points is required.
    ds : int, default 1
        The second interval in which to calculate the spline for.
        E.g. ds = 1 is "calculate a spline point at every 1 second delta"

    Returns
    -------
    gpd.GeoDataFrame of all points in the interpolated spline.
    Columns: point_dt, geometry

    Raises
    ------
    AssertionError if there is fewer than 1 Track point.
    """
    # Calculate the order of polynomial to fit to the spline. The maximum is a cubic spline. If there are fewer than
    #  3 points, a cubic spline cannot be fit and lower order must be chosen.
    assert points.shape[0] > 1, "A minimum of 2 points is required for calculate a spline."
    k = min(points.shape[0] - 1, 3)

    points.sort_values(by='point_dt', ascending=True, inplace=True)
    starttime = points.point_dt.iat[0]
    endtime = points.point_dt.iat[-1]
    flight_times = (points.point_dt - starttime).dt.total_seconds().values  # Seconds after initial point

    coords = [points.geometry.x, points.geometry.y, points.z] if 'z' in points else [points.geometry.x, points.geometry.y]
    tck, u = interpolate.splprep(x=coords, u=flight_times, k=k)

    # Parametric interpolation on the time interval provided.
    duration = (endtime - starttime).total_seconds()
    tnew = np.arange(0, duration + ds, ds)
    spl_out = interpolate.splev(tnew, tck)
    track_spline = gpd.GeoDataFrame({'point_dt': [starttime + dt.timedelta(seconds=offset) for offset in tnew]},
                                    geometry=[Point(xyz) for xyz in zip(spl_out[0], spl_out[1], spl_out[2])],
                                    crs=points.crs)
    return track_spline


def audible_time_delay(points: gpd.GeoDataFrame, time_col: str, target: Point,
                       m1: float = 343., drop_cols: bool = False) -> gpd.GeoDataFrame:
    """
    Given a set of points and a target location, calculate when a sound made at each point could be heard at
    the target.

    **IMPORTANT**: The points GeoDataFrame and the target Point should be in the same crs for accurate calculations.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        A gpd.GeoDataFrame of sound location points.
    time_col : str
        Name of the column in the points gpd.GeoDataFrame with time of sound occurrence at each point.
    target : Point
        The target point.
    m1 : float, default 343 m/s
        The speed of sound to use for calculations. Make sure this value uses the same units as the crs of
        the points GeoDataFrame and the target Point.
    drop_cols : bool, default False
        If True, drop the intermediate columns used to determine time of audibility.

    Returns
    -------
    The points GeoDataFrame with added columns:
    Standard: time_audible
    Optional: distance_to_target, audible_delay_sec
    """
    points['distance_to_target'] = points.geometry.apply(lambda geom: target.distance(geom))
    points['audible_delay_sec'] = points['distance_to_target'] / m1
    points['time_audible'] = points.apply(lambda row: row[time_col] + dt.timedelta(seconds=row.audible_delay_sec), axis=1)

    if drop_cols:
        points.drop(['distance_to_target', 'audible_delay_sec'], inplace=True)

    return points


def build_src_point_mesh(area: 'Polygon', density: int = 48, altitude: Optional[int] = None) -> List['Point']:
    """
    Given a polygon and a density, create a square mesh of evenly spaced points throughout the polygon.

    Parameters
    ----------
    area : Polygon
        A shapely Polygon of the area to create the square point mesh over.
    density : int
        The number of points along each mesh axis. The mesh will contain density x density points.
    altitude : int, default None
        A standard altitude to apply to every point in the mesh.

    Returns
    -------
    mesh points : List[Point]
        A list of shapely Points in the mesh.
    """
    # Start out with a grid of N = density x density points. Polygon bounds:  (minx, miny, maxx, maxy)
    x = np.linspace(area.bounds[0], area.bounds[2], density)
    y = np.linspace(area.bounds[1], area.bounds[3], density)
    x_ind, y_ind = np.meshgrid(x, y)

    # Create an array of mesh points. np.ravel linearly indexes an array into a row.
    mesh_points = np.array([np.ravel(x_ind), np.ravel(y_ind)]).T

    # Convert coordinate tuples into shapely points.
    mesh_points = [Point(point[0], point[1]) if not altitude
                   else Point(point[0], point[1], altitude) for point in mesh_points]

    return mesh_points
