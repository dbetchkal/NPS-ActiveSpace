import logging
from typing import List, Optional, TYPE_CHECKING

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from nps_active_space import _ACTIVE_SPACE_DIR
from nps_active_space.utils import coords_to_utm, Microphone

if TYPE_CHECKING:
    from sqlalchemy.engine import Engine


__all__ = [
    'get_deployment',
    'get_logger',
    'get_omni_sources',
    'query_tracks'
]


def get_deployment(unit: str, site: str, year: int, filename: str) -> Microphone:
    """
    Obtain all metadata for a specific microphone deployment from a metadata file.

    Parameters
    ----------
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA'
    year : int
        Deployment year. YYYY
    filename : str
        Absolute path to microphone deployment metadata text file. '/path/to/metadata.txt'

    Returns
    -------
    mic : Microphone
        A Microphone object containing the mic deployment site metadata from the specific unit/site/year combination.
    """
    metadata = pd.read_csv(filename, delimiter='\t', encoding='ISO-8859-1')
    site_meta = metadata.loc[(metadata['unit'] == unit) & (metadata['code'] == site) & (metadata['year'] == year)]

    # Microphone coordinates are stored in WGS84, epsg:4326
    mic = Microphone(
        lat=site_meta.lat.iat[0],
        lon=site_meta.long.iat[0],
        z=site_meta.elevation.iat[0],
        name=f"{unit}{site}{year}"
    )

    return mic


def query_tracks(engine: 'Engine', start_date: str, end_date: str,
                 mask: Optional[gpd.GeoDataFrame] = None) -> gpd.GeoDataFrame:
    """
    Query flight tracks from the FlightsDB for a specific date range and optional within a specific area.

    Parameters
    ----------
    engine : sqlalchemy Engine
        SQLAlchemy Engine instance for connecting to the overflights DB.
    start_date : str
        ISO date string (YYYY-mm-dd) indicating the beginning of the date range to query within
    end_date : str
        ISO date string (YYYY-mm-dd) indicating the end of the date range to query within
    mask : gpd.GeoDataFrame, default None
        Geopandas.GeoDataframe instance to spatially filter query results.

    Returns
    -------
    data : gpd.GeoDataFrame
        A GeoDataFrame of flight track points.
    """
    wheres = [f"fp.ak_datetime::date BETWEEN '{start_date}' AND '{end_date}'"]

    if mask is not None:
        if mask.crs.to_epsg() != 4326:  # If mask is not already in WGS84, project it.
            mask = mask.to_crs(epsg='4326')
        mask['dissolve_field'] = 1
        mask_wkt = mask.dissolve(by='dissolve_field').squeeze()['geometry'].wkt
        wheres.append(f"ST_Intersects(geom, ST_GeomFromText('{mask_wkt}', 4326))")

    query = f"""
        SELECT
            f.id as flight_id,
            fp.altitude_ft * 0.3048 as altitude_m,
            fp.ak_datetime,
            fp.geom, 
            date_trunc('hour', fp.ak_datetime) as ak_hourtime
        FROM flight_points as fp
        JOIN flights f ON f.id = fp.flight_id
        WHERE {' AND '.join(wheres)}
        ORDER BY fp.ak_datetime asc
        """

    flight_tracks = gpd.GeoDataFrame.from_postgis(query, engine, geom_col='geom', crs='epsg:4326')

    data = flight_tracks.loc[~(flight_tracks.geometry.is_empty)]
    return data


class _TqdmStreamHandler(logging.StreamHandler):
    """
    A Logger Stream Hanlder so Tqdm loading bars work with python loggers.
    https://stackoverflow.com/questions/14897756/python-progress-bar-through-logging-module
    """
    @classmethod
    def write(cls, msg):
        tqdm.write(msg, end='')


def get_logger(name: str, level: str = 'INFO') -> logging.Logger:
    """
    General purpose function for creating a console logger.

    Parameters
    ----------
    name : str
        Logger name
    level : str, default INFO
        Logger message severity

    Returns
    -------
    logger : logging.Logger
        A python logger object
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = _TqdmStreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def get_omni_sources(lower: int, upper: int) -> List[str]:
    """
    Get a list of omni source files for tuning NMSim within a specific gain range.
    Source files are provided in the data directory for gains between -30 and +50.

    Parameters
    ----------
    lower: int
        The lowest gain omni source file to pull.
    upper : int
        The high gain omni source file to pull

    Returns
    -------
    A list of omni source files within the specified gain range.

    Raises
    ------
    AssertionError if the lower or upper gain bound is out of range or of the upper gain bound is lower than
    the lower gain bound.
    """
    assert -30 <= upper <= 50 and -30 <= lower <= 50 and upper > lower, "Bounds must be between [-30, 50]"

    omni_source_dir = f"{_ACTIVE_SPACE_DIR}/data/tuning"
    omni_sources = []

    for i in range(lower*10, upper*10+5, 5):
        if i < 0:
            omni_sources.append(f"{omni_source_dir}/O_{i:04}.src")
        elif i >= 0:
            omni_sources.append(f"{omni_source_dir}/O_+{i:03}.src")

    return omni_sources
