import glob
import logging
import os
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING, Union
import pandas as pd
import matplotlib.pyplot as plt

import geopandas as gpd
import numpy as np
import sqlalchemy
from sqlalchemy.engine import URL
from tqdm import tqdm
import re
import rasterio
import rasterio.plot
from pyproj import Transformer
from shapely.geometry import box

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.active_space import LayeredActiveSpace
from nps_active_space.utils.models import Adsb, EarlyAdsb, Microphone, Annotations
from nps_active_space.setup.site_decoder import decode_sit_geographic_coords, parse_sit_coords_line
from nps_active_space.setup.site_writer import (
    NMSIM_SITES_DIR,
    deployment_sit_name,
    sit_file_path,
)
from nps_active_space.utils.computation import NMSIM_bbox_utm

if TYPE_CHECKING:
    from sqlalchemy.engine import Engine


__all__ = [
    'omni_to_gain',
    'load_activespace',
    'load_layered_activespace',
    'load_DEM',
    'get_elevation',
    'load_study_area',
    'get_deployment',
    'create_overflights_engine',
    'query_tracks',
    'query_adsb',
    'load_annotations',
    'get_logger',
    'get_omni_sources',
    'estimate_line_count'
]


def _glob_path_or_raise(pattern: str, *, description: str) -> str:
    """Return the path matching ``pattern``, or raise ``FileNotFoundError``."""
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(
            f"No {description} found matching:\n  {pattern}"
        )
    return matches[0]


def omni_to_gain(omni_source: str) -> float:
    """
    Converts an omni source name to the corresponding gain.
    Pure omni strings ("O_+125") or paths "directory/O_+125.src" can be passed, since regex is used.
    """
    match = re.search(r"O_([+-]\d\d\d)", omni_source)
    return int(match.group(1)) / 10


def load_layered_activespace(project_dir, unit, site, year, gain=None, crs="epsg:4326"):
    prefix = os.path.join(project_dir, f"{unit}{site}", "Output_Data", "ACTIVESPACES")
    layer_dirs = {}
    output_dirs = glob.glob(os.path.join(prefix, f"{unit}{site}{year}_*m"))
    for dir in output_dirs:
        altitude = int(os.path.basename(dir).split("_")[1].split("m")[0])
        layer_dirs[altitude] = dir
    study_area = load_studyarea(project_dir, unit, site, year)
    return LayeredActiveSpace(unit+site+year, layer_dirs, study_area, gain, crs)


def load_activespace(project_dir, unit, site, year, gain, altitude_m=None, crs=None):
    """
    Load in the active space for a given unit, site, year, and gain

    Parameters
    ----------
    project_dir: str
        Path to the project directory (contains site subfolders named e.g. DENATRLA/)
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA', '009'
    year : int
        Deployment year. YYYY
    gain : float
        The optimal gain, or scaling factor, of the active space, determined during ground truthing
    altitude_m : int, default None
        The altitude of the active space, in meters. If not provided, the middle altitude of existing active
        space layers is used.
    crs : string
        Optional argument to provide a coordinate reference system to convert the active space to (e.g.  'epsg:26905'). Defaults to None

    Returns
    -------
    active_space : `gpd.GeoDataFrame`
        A dataframe containing the geometry of the active space. Can be a single polygon or a multipolygon.
    """

    prefix = os.path.join(project_dir, unit + site, "Output_Data", "ACTIVESPACES")

    # pick middle altitude if no altitude provided
    if altitude_m is None:
        altitude_dirs = glob.glob(os.path.join(prefix, f"{unit}{site}{year}_*m"))
        altitudes = []
        for dir in altitude_dirs:
            altitudes.append(int(os.path.basename(dir).split("_")[1].split("m")[0]))
        altitudes.sort()
        altitude_m = altitudes[len(altitudes) // 2]
        print(f"No altitude specified, using {altitude_m}m")

    # read activespace
    sign = "-" if gain < 0 else "+"
    gain_string = str(np.abs(int(10*gain))).zfill(3)
    usy = f"{unit}{site}{year}"
    path = os.path.join(prefix, f"{usy}_{altitude_m}m", f"{usy}_O_{sign}{gain_string}.geojson")
    active_space = gpd.read_file(path)

    if crs is not None:
        active_space = active_space.to_crs(crs)

    return active_space


def load_DEM(project_dir: str, unit: str, site: str):
    """Loads the `NMSIM` digital elevation model using the project info

    Parameters
    ----------
    project_dir: str
        Path to the project directory (contains site subfolders named e.g. DENATRLA/)
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA', '009'
    year : int
        Deployment year. YYYY

    Returns
    -------
    DEM
        A rasterio Dataset object for reading the DEM data
    """

    dem_glob = os.path.join(
        project_dir, unit + site, "Input_Data", "01_ELEVATION", "elevation_m_nad83_utm*.tif"
    )
    raster_path = _glob_path_or_raise(
        dem_glob,
        description=(
            f"elevation DEM for {unit}{site} "
            f"(expected GeoTIFF in Input_Data/01_ELEVATION/elevation_m_nad83_utm*.tif)"
        ),
    )
    return rasterio.open(raster_path)


def get_elevation(DEM, lon: float, lat: float) -> float:
    """Read elevation at a certain lat/lon from a DEM.
    """
    proj = Transformer.from_crs("epsg:4326", DEM.crs, always_xy=True)
    x, y = proj.transform(lon, lat)
    row, col = DEM.index(x, y)
    elevation = DEM.read(1)[row, col]
    return elevation


def load_studyarea(project_dir: str, unit: str, site: str, year: int, crs: str = None):
    """
    Load in the study area that contains the active space for a given unit, site, and year. 
    At present this function is overridden by `query_tracks()`, which is set to provide a study area = active space buffered by 25km.

    Parameters
    ----------
    project_dir: str
        Path to the project directory (contains site subfolders named e.g. DENATRLA/)
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA', '009'
    year : int
        Deployment year. YYYY
    crs : string
        Optional argument to provide a coordinate reference system to convert the study area to (e.g.  'epsg:26905'). Defaults to None

    Returns
    -------
    study_area : `gpd.GeoDataFrame`
        A dataframe containing the geometry of the study area. A single polygon.
    """

    study_glob = os.path.join(
        project_dir, unit + site, f"{unit}{site}*study*area*.shp"
    )
    study_area_path = _glob_path_or_raise(
        study_glob,
        description=f"study area shapefile for {unit}{site}{year}",
    )
    study_area = gpd.read_file(study_area_path)

    if crs is not None:
        study_area = study_area.to_crs(crs)

    return study_area


def get_deployment(project_dir: str, unit: str, site: str, year: int, elevation: bool = True) -> Microphone:
    """
    Obtain all metadata for a specific microphone deployment from a metadata file.

    Parameters
    ----------
    project_dir: str
        Path to the project directory (contains site subfolders named e.g. DENATRLA/)
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA', '009'
    year : int
        Deployment year. YYYY
    elevation : bool, default True
        If True, the microphone z value will be set to its elevation. If False, the microphone z value will be
        set to the microphone's height from the ground.

    Returns
    -------
    mic : Microphone
        A Microphone object containing the mic deployment site metadata from the specific unit/site/year combination.
    """

    # Read the .SIT file containing x, y, and z AGL in NMSIM's CRS
    mic_location_path = sit_file_path(
        Path(project_dir) / (unit + site),
        deployment_sit_name(unit, site, year),
    )
    if not mic_location_path.is_file():
        raise FileNotFoundError(
            f"Microphone site file not found:\n  {mic_location_path}\n"
            f"Expected NMSIM .sit file at {NMSIM_SITES_DIR}/{deployment_sit_name(unit, site, year)}.sit"
        )
    raw_text = Path(mic_location_path).read_text()
    sit_lines = raw_text.splitlines()
    if len(sit_lines) < 3:
        raise ValueError(
            f"Microphone site file has invalid format (expected coordinates on line 3):\n"
            f"  {mic_location_path}"
        )
    coords_line = sit_lines[2]
    x, y, z_agl = parse_sit_coords_line(coords_line)

    # get lat/lon so we can initialize a Microphone object, need to get crs of the .SIT file first
    study_area = load_studyarea(project_dir, unit, site, year)
    mic_crs = NMSIM_bbox_utm(study_area)
    lon, lat = decode_sit_geographic_coords(
        x,
        y,
        study_area,
        sit_path=mic_location_path,
        unit=unit,
        site=site,
        year=year,
    )

    # calculate elevation from the DEM if necessary
    if elevation:
        DEM = load_DEM(project_dir, unit, site)
        z = z_agl + get_elevation(DEM, lon, lat)
    else:
        z = z_agl

    mic = Microphone(
        name=f"{unit}{site}{year}",
        lat=lat,
        lon=lon,
        z=z,
        crs=mic_crs
    )

    return mic


def create_overflights_engine(db: dict[str, str]) -> "Engine":
    """SQLAlchemy engine for the overflights PostgreSQL database."""
    return sqlalchemy.create_engine(
        URL.create(
            drivername="postgresql",
            username=db["username"],
            password=db["password"],
            host=db["host"],
            port=db["port"],
            database=db["name"],
        )
    )


def query_tracks(engine: 'Engine', start_date: str, end_date: str,
                 mask: Optional[gpd.GeoDataFrame] = None,
                 mask_buffer_distance: Optional[int] = None) -> gpd.GeoDataFrame:
    """
    Query flight tracks from the FlightsDB for a specific date range and optional within a specific area.

    Parameters
    ----------
    engine : sqlalchemy Engine
        SQLAlchemy Engine instance for connecting to the overflights DB.
    start_date : str
        ISO date string (YYYY-mm-dd) indicating the beginning of the date range to query within
    end_date : str
        ISO date string (YYYY-mm-dd) indicating the end of the date range to query within, inclusive
    mask : gpd.GeoDataFrame, default None
        Geopandas.GeoDataframe instance to spatially filter query results.

    Returns
    -------
    data : gpd.GeoDataFrame
        A GeoDataFrame of flight track points.
    """
    wheres = [f"fp.ak_datetime::date BETWEEN '{start_date}' AND '{end_date}'"]  # start and end date are inclusive

    if mask is not None:
        if mask.crs.to_epsg() != 4326:  # If mask is not already in WGS84, project it.
            mask = mask.to_crs(epsg='4326')
        if mask_buffer_distance:
            ak_albers_mask = mask.to_crs(epsg=3338)
            mask.geometry = ak_albers_mask.buffer(
                mask_buffer_distance).to_crs(epsg=4326)
        wheres.append(
            f"ST_Intersects(geom, ST_GeomFromText('{mask.union_all().wkt}', 4326))")

    query = f"""
        SELECT
            f.flight_id as flight_id,
            fp.altitude_ft * 0.3048 as altitude_m,
            fp.ak_datetime,
            fp.geom, 
            date_trunc('hour', fp.ak_datetime) as ak_hourtime
        FROM flight_points as fp
        JOIN flights f ON f.id = fp.flight_id
        WHERE {' AND '.join(wheres)}
        ORDER BY fp.ak_datetime asc
        """

    flight_tracks = gpd.GeoDataFrame.from_postgis(
        query, engine, geom_col='geom', crs='epsg:4326')

    data = flight_tracks.loc[~(flight_tracks.geometry.is_empty)]
    return data


def query_adsb(adsb_path: str,  start_date: str, end_date: str,
               mask: Optional[gpd.GeoDataFrame] = None,
               mask_buffer_distance: Optional[int] = None,
               exclude_early_ADSB: Optional[bool] = False) -> Union[Adsb, EarlyAdsb]:
    """
    Query flight tracks from ADSB files for a specific date range and optional within a specific area.

    Parameters
    ----------
    adsb_path : str
        Absolute path to a directory with adsb data files to read in.
    start_date : str
        ISO date string (YYYY-mm-dd) indicating the beginning of the date range to query within
    end_date : str
        ISO date string (YYYY-mm-dd) indicating the end of the date range to query within, inclusive
    mask : gpd.GeoDataFrame, default None
        Geopandas.GeoDataframe instance to spatially filter query results.

    Returns
    -------
    adsb : ADSB or EarlyADSB
        An ADSB or EarlyADSB object of flight track points.
    """

    if mask is not None:
        # if not mask.crs.to_epsg() == 4326:  # If mask is not already in WGS84, project it.
        #     mask = mask.to_crs(epsg='4326')
        if mask_buffer_distance:
            ak_albers_mask = mask.to_crs(epsg=3338)
            mask.geometry = ak_albers_mask.buffer(
                mask_buffer_distance).to_crs(epsg=4326)
        mask = mask.to_crs("epsg:4326")

    # ADSB file formats changed after 2019.
    if (int(start_date[:4]) <= 2019) & (exclude_early_ADSB == False):
        adsb_files = glob.glob(os.path.join(adsb_path, "*.txt"))
        assert len(adsb_files) > 0, f"No ADSB files found in {adsb_path}"
        adsb = EarlyAdsb(adsb_files)
        adsb.set_crs(epsg='4326', inplace=True)
        if mask is not None:
            adsb = gpd.clip(adsb, mask)
    else:
        adsb_files = glob.glob(os.path.join(adsb_path, "*.TSV"))
        assert len(adsb_files) > 0, f"No ADSB files found in {adsb_path}"
        start_dt = pd.Timestamp(start_date)
        # Adsb() takes an exclusive end date (midnight) so it can do comparisons like: time < end_dt
        # so convert inclusive end date to exclusive end date by adding one day
        end_dt = pd.Timestamp(end_date) + pd.Timedelta(days=1)
        adsb = Adsb(adsb_files, mask, start_dt, end_dt)
    
    assert not adsb.empty, f"No ADSB data loaded for {start_date} to {end_date}, please check the time period and/or region mask"

    end_ts = pd.Timestamp(end_date) + pd.Timedelta(days=1)
    start_ts = pd.Timestamp(start_date)
    adsb = adsb.loc[(adsb["TIME"] >= start_ts) & (adsb["TIME"] < end_ts)]
    adsb = adsb.loc[~(adsb.geometry.is_empty)]
    return adsb


def load_annotations(project_dir: str, unit: str, site: str, year: str, only_valid: bool = True):
    """Utility for locating and loading ground-truthing annotation files in a directory.
    If multiple files exist, they are combined into one GeoDataFrame.

    Parameters
    ----------
    project_dir: str
        Path to the project directory (contains site subfolders named e.g. DENATRLA/)
    unit : str
        Four letter park service unit code E.g. 'DENA'
    site : str
        Deployment site character code. E.g. 'TRLA', '009'
    year : int
        Deployment year. YYYY
    only_valid : bool, default True
        Whether to only load valid annotations.
    
    Returns
    -------
    annotations: Annotations
        An Annotations object containing the loaded annotations.
    """
    # Verify that annotation files exist for the unit/site location. If they do exist, load them into memory.
    annotation_files = glob.glob(f"{project_dir}/{unit}{site}/{unit}{site}{year}*saved_annotations*.geojson")
    print("Found these annotation files:", list(map(lambda f: os.path.basename(f), annotation_files)))
    if len(annotation_files) == 0:
        return gpd.GeoDataFrame()
    annotations = []
    for file in tqdm(annotation_files, desc='Loading annotation files', unit='files', colour='white'):
        annotations.append(Annotations(file, only_valid=only_valid))
    return pd.concat(annotations, ignore_index=True)


class _TqdmStream:
    """
    A Logger Stream so Tqdm loading bars work with python loggers.
    https://github.com/tqdm/tqdm/issues/313#issuecomment-346819396
    """
    def write(cls, msg: str):
        tqdm.write(msg, end='', )
    write = classmethod(write)


class _BufferedHandler(logging.Handler):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.logs = []

    def emit(self, record):
        log_entry = self.format(record)
        self.logs.append(log_entry)

    def save(self, filename):
        with open(filename, 'w') as f:
            f.write('\n'.join(self.logs))

    def clear(self):
        self.logs.clear()


def get_logger(name: str, verbose: bool = False, logfile: str = None, make_log_buffer=False) -> logging.Logger:
    """
    General purpose function for creating a logger.

    Parameters
    ----------
    name : str
        Logger name
    verbose : bool, default False
        If verbose is True, all messages will be printed to the console and saved in a log file.
        If verbose is False, only high-level messages will be printed to the console; everything is still saved to a log file.
    logfile: str, default None
        Path to a file to save the log in. If None, the log will not be saved to a file.
    make_log_buffer: bool, default False
        If True, save a full copy of the log in memory, that can later be saved to a file

    Returns
    -------
    logger : logging.Logger
        A python logger object
    
    If make_log_buffer is True:
    Tuple of (logger , log_buffer). You can call log_buffer.save(filename) to save the buffered log to a file.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)  # let the handlers do the filtering

    # Clear existing handlers from previous loggers with the same name
    if logger.hasHandlers():
        logger.handlers.clear()

    # Handler to print to the console
    console_handler = logging.StreamHandler(stream=_TqdmStream)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(console_handler)

    # Handler to save to a file during logging
    if logfile is not None:
        file_handler = logging.FileHandler(logfile)
        # always print everything to the log file
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(
            fmt='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
        logger.addHandler(file_handler)

    # Handler to just remember logs, for future saving
    if make_log_buffer:
        buffer_handler = _BufferedHandler()
        buffer_handler.setLevel(logging.DEBUG)
        buffer_handler.setFormatter(logging.Formatter(
            fmt='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
        logger.addHandler(buffer_handler)

    if make_log_buffer:
        return logger, buffer_handler
    return logger


def get_omni_sources(lower: float, upper: float) -> List[str]:
    """
    Get a list of omni source files for tuning NMSim within a specific gain range.
    Source files are provided in the data directory for gains between -30 and +50.

    O_+005 = .5
    O_+050 = 5
    O_+500 = 50

    Parameters
    ----------
    lower: float
        The lowest gain omni source file to pull.
    upper : float
        The high gain omni source file to pull

    Returns
    -------
    A list of omni source files within the specified gain range.

    Raises
    ------
    AssertionError if the lower or upper gain bound is out of range or of the upper gain bound is lower than
    the lower gain bound.
    """
    assert -30 <= upper <= 50 and - \
        30 <= lower <= 50 and upper >= lower, "Bounds must be ordered and between [-30, 50]."
    assert upper % .5 == 0, "Invalid upper limit. Value must be divisible by 0.5."
    assert lower % .5 == 0, "Invalid lower limit. Value must be divisible by 0.5."

    omni_source_dir = os.path.join(ACTIVE_SPACE_DIR, "data", "tuning")
    omni_sources = []

    for i in range(int(lower*10), int(upper*10+5), 5):
        if i < 0:
            omni_sources.append(os.path.join(omni_source_dir, f"O_{i:04}.src"))
        elif i >= 0:
            omni_sources.append(os.path.join(omni_source_dir, f"O_+{i:03}.src"))

    return omni_sources


def estimate_line_count(filename, sample_size=1024 * 1024):
    """Use a 1MB sample to estimate the number of lines in a large file"""
    file_size = os.path.getsize(filename)
    with open(filename, 'rb') as f:
        sample = f.read(sample_size)
    newlines = sample.count(b'\n')
    if not newlines:
        return 0
    return int((file_size / sample_size) * newlines)



def plot_activespace_fit(project_dir, unit, site, year, gain, altitude_m=None,
                         ax=None, dem=None, mic=None, active=None, annotations=None):
    if ax is None:
        fig, ax = plt.subplots()

    if dem is None:
        dem = load_DEM(project_dir, unit, site)
    if mic is None:
        mic = get_deployment(project_dir, unit, site, year)
    if active is None:
        active = load_activespace(project_dir, unit, site, year, gain, altitude_m)
    if annotations is None:
        annotations = load_annotations(project_dir, unit, site, year)
    
    mic = mic.to_crs(dem.crs)
    

    ax.set_title(f"{unit}{site}{year} Gain {gain}dB")
    rasterio.plot.show(dem, ax=ax, alpha=.4, cmap='Blues_r')

    if not active.empty:
        active = active.to_crs(dem.crs)
        active.boundary.plot(ax=ax, color="black", label=f"{gain}", zorder=5)

    dem_extent = box(dem.bounds.left, dem.bounds.bottom, dem.bounds.right, dem.bounds.top)
    annotations = annotations.clip(dem_extent)
    if not annotations.empty and annotations["valid"].sum() > 0:
        valid_segments = annotations[annotations.valid]
        color = valid_segments.apply(lambda x: "deepskyblue" if x["audible"] else "red", axis=1)
        valid_segments.plot(
            ax=ax,
            color=color,
            alpha=0.5,
            markersize=3,
            zorder=2
        )
    
    mic.plot(ax=ax, color="black", markersize=15, marker='o', zorder=10)

