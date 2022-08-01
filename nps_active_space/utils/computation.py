import datetime as dt
from typing import Iterable, Optional, Tuple, Union

import geopandas as gpd
import numpy as np
from scipy import interpolate
from shapely.geometry import Point

from nps_active_space.utils.models import Tracks


__all__ = [
    'audible_time_delay',
    'climb_angle',
    'coords_to_utm',
    'interpolate_spline'
]


def coords_to_utm(lat: float, lon: float) -> str:
    """
    Takes the latitude and longitude of a point and outputs the EPSG code corresponding to the UTM zone of the point.

    Parameters
    ----------
    lat : float
        Latitude of a point in decimal degrees in a WGS84 projection (EPSG:4326)
    lon: float
        Longitude of a point in decimal degrees in a WGS84 projection (EPSG:4326)

    Returns
    -------
    utm_proj : str
        UTM zone projection name (e.g. 'epsg:26905' for UTM 5N)

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


def interpolate_spline(points: Tracks, ds: float = 1):
    """
    # TODO

    Parameters
    ----------
    points: Tracks
    ds : float
        delta seconds
    """
    # Calculate the order of polynomial to fit to the spline. The maximum is a cubic spline. If there are fewer than
    # 3 points, a cubic spline cannot be fit and lower order must be chosen.
    assert points.shape[0] > 1, "A minimum of 2 points is required for calculate a spline."
    k = min(points.shape[0] - 1, 3)

    points.sort_values(by='point_dt', ascending=True, inplace=True)
    starttime = points.point_dt.iat[0]
    endtime = points.point_dt.iat[-1]
    flight_times = (points.point_dt - starttime).dt.total_seconds().values  # Seconds after initial point

    coords = [points.geometry.x, points.geometry.y, points.z] if 'z' in points else [points.geometry.x, points.geometry.y]
    tck, u = interpolate.splprep(coords, u=flight_times, k=k)

    # parametric interpolation on the time interval provided
    duration = (endtime - starttime).total_seconds()
    tnew = np.arange(0, duration + ds, ds)  # new spline fit with dS
    spl_out = interpolate.splev(tnew, tck)

    track_spline = gpd.GeoDataFrame({'point_dt': [starttime + dt.timedelta(seconds=offset) for offset in tnew]},
                                    geometry=[Point(xyz) for xyz in zip(spl_out[0], spl_out[1], spl_out[2])],
                                    crs=points.crs)

    return track_spline


def audible_time_delay(points: gpd.GeoDataFrame, time_col: str, target: Point, m1: float = 343.,
                       drop_cols: Optional[bool] = False):
    # TODO: same crs??
    # TODO: does this give the right value?

    points['distance_to_target'] = points.geometry.apply(lambda geom: target.distance(geom))  # TODO What are the metrics here.
    points['audible_delay_sec'] = points['distance_to_target'] / m1
    points['time_audible'] = points.apply(lambda row: row[time_col] + dt.timedelta(seconds=row.audible_delay_sec), axis=1)

    if drop_cols:
        points.drop(['distance_to_target', 'audible_delay_sec'], inplace=True)

    return points