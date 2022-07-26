import datetime as dt
from typing import Iterable, Optional, Tuple, Union

import geopandas as gpd
import numpy as np
from scipy import interpolate
from shapely.geometry import Point


__all__ = [
    'calculate_spline',
    'climb_angle',
    'coords_to_utm'
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


# TODO
def calculate_spline(idx, points, crs, dT=1):
    pass
