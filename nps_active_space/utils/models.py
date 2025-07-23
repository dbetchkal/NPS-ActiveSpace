import time
from types import GeneratorType
import concurrent.futures
from tzwhere import tzwhere
from tqdm import tqdm
from pyproj import Transformer
import datetime as dt
import glob
import pytz
import re
import os
import csv
import json
import matplotlib.pyplot as plt
from warnings import warn
from dataclasses import dataclass, field
from typing import List, Optional, Union
from decimal import Decimal, ROUND_FLOOR

import geopandas as gpd
from shapely.geometry import Point, box
import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True
pd.set_option('future.no_silent_downcasting', True)

__all__ = [
    'Adsb',
    'Ais',
    'Annotations',
    'EarlyAdsb',
    'Microphone',
    'Nvspl',
    'Tracks'
]


@dataclass
class Microphone:
    """
    An object representing a microphone deployment location.

    Parameters
    ----------
    name : str
        A name for the Microphone instance.
    lat : float
        The latitude of the microphone deployment location in WGS84 (epsg:4326)
    lon : float
        The longitude of the microphone deployment location in WGS84 (epsg:4326)
    z : float
        The elevation of the microphone deployment location in meters.
    crs : str, default None
        Epsg projected coordinated system to calculate the x, y values in. E.g. 'epsg:4326'
        Latitude and Longitude will not be projected if no crs is provided.

    Instance Variables
    ------------------
    x : float
        The longitude value projected into the current crs.
    y : float
        The latitude value projected into the current crs.
    """
    name: str
    lat: float
    lon: float
    z: float
    crs: str = None
    x: float = field(init=False)
    y: float = field(init=False)

    def __repr__(self):
        return f"Microphone(name={self.name})"

    def __post_init__(self):
        """Set x,y coordinates and instance name."""
        if self.crs:
            self.to_crs(self.crs)

    def to_crs(self, crs: str, inplace: bool = False) -> Optional['Microphone']:
        """
        Project instance x,y values to a new coordinate system.

        Parameters
        ----------
        crs : str
            The coordinate system to project the instance to.
                Format: epsg:XXXX. E.g. epsg:26906
        inplace : bool, default False
            If True, crs will be updated and no instance will be returned.
            If False, crs will be updated an the updated instance will be returned.
        """
        projection = Transformer.from_crs('epsg:4326', crs, always_xy=True)
        self.x, self.y = projection.transform(self.lon, self.lat)
        self.crs = crs
        if not inplace:
            return self

    def plot(self, **kwargs):
        """Plot this microphone's location using geopandas

        Parameters
        ----------
        **kwargs
            Any parameters that GeoDataFrame.plot() accepts
        """
        gdf = gpd.GeoDataFrame(geometry=[Point(self.x, self.y)], crs=self.crs)
        gdf.plot(**kwargs)


class Nvspl(pd.DataFrame):
    """
    A pandas DataFrame wrapper class to ensure consistent NVSPL data.

    Parameters
    ----------
    filepaths_or_data : List, str, or pd.DataFrame
        A directory containing NVSPL files, a list of NVSPL files, or an existing pd.DataFrame of NVSPL data.
    """

    standard_fields = {
        'SiteID', 'dbA', 'dbC', 'dbF',
        'Voltage', 'WindSpeed', 'WindDir', 'TempIns',
        'TempOut', 'Humidity', 'INVID', 'INSID',
        'GChar1', 'GChar2', 'GChar3', 'AdjustmentsApplied',
        'CalibrationAdjustment', 'GPSTimeAdjustment',
        'GainAdjustment', 'Status'
    }

    octave_regex = re.compile(r"^H[0-9]+$|^H[0-9]+p[0-9]$")

    def __init__(self, filepaths_or_data: Union[List[str], str, pd.DataFrame]):
        data = self._read(filepaths_or_data)
        super().__init__(data=data)

    def parseNvspl(self, nvsplFileEntry, state=(None, None, 1)):

        timestamps, columns, index_index = state

        df = pd.read_csv(str(nvsplFileEntry),
                         engine='c',
                         parse_dates=True,
                         index_col=index_index,
                         usecols=columns
                         )

        # Make column names slightly nicer
        df.index.name = "date"
        renamedColumns = {column: column.replace('H', '').replace(
            'p', '.') for column in df.columns if re.match(r"H\d+p?\d*", column) is not None}
        df.rename(columns=renamedColumns, inplace=True)

        # Coerce numeric columns to floats, in case of "-Infinity" values
        try:
            numericCols = [
                '12.5', '15.8', '20', '25', '31.5', '40', '50', '63', '80', '100',
                '125', '160', '200', '250', '315', '400', '500', '630', '800', '1000',
                '1250', '1600', '2000', '2500', '3150', '4000', '5000', '6300', '8000',
                '10000', '12500', '16000', '20000', 'dbA', 'dbC', 'dbF',
                'Voltage', 'WindSpeed', 'WindDir', 'TempIns', 'TempOut', 'Humidity'
            ]
            presentNumericCols = df.columns.intersection(numericCols)
            if len(presentNumericCols) > 0:
                df[presentNumericCols].astype(
                    'float32', copy=False, errors='ignore')

        except KeyError:
            pass

        self._validate(df, False)
        return df

    def _read(self, filepaths_or_data: Union[List[str], str, pd.DataFrame, GeneratorType]):
        """
        Read in and validate the NVSPL data.

        # TODO: for speed and memory improvements, use usecols, define datatypes, and drop empty columns.

        Parameters
        ----------
        filepaths_or_data : List, str, or pd.DataFrame
            A directory containing NVSPL files, a list of NVSPL files, or an existing pd.DataFrame of NVSPL data.

        Raises
        ------
        AssertionError if directory path or file path does not exists or is of the wrong format.
        """
        if isinstance(filepaths_or_data, pd.DataFrame):
            self._validate(filepaths_or_data.columns)
            data = filepaths_or_data

        else:

            if str(type(filepaths_or_data)) == "<class 'iyore.Subset'>":
                filepaths_or_data = [str(entry)
                                     for entry in list(iter(filepaths_or_data))]

            if isinstance(filepaths_or_data, str):
                assert os.path.isdir(
                    filepaths_or_data), f"{filepaths_or_data} does not exist."
                filepaths_or_data = glob.glob(f"{filepaths_or_data}/*.txt")

            else:
                for file in filepaths_or_data:
                    assert os.path.isfile(file), f"{file} does not exist."
                    assert file.endswith(
                        '.txt'), f"Only .txt NVSPL files accepted."

            with concurrent.futures.ThreadPoolExecutor() as pool:
                parts = list(tqdm(pool.map(self.parseNvspl, filepaths_or_data), total=len(
                    filepaths_or_data), unit="NVSPL files"))

            data = pd.concat(parts)

        octave_columns = {c: c.replace('H', '').replace(
            'p', '.') for c in filter(self.octave_regex.match, data.columns)}
        data.rename(columns=octave_columns, inplace=True)

        # we deliberately sort the DatetimeIndex to ensure it is monotonic
        # this avoids a `KeyError` when selecting using position report timestamps later on
        data.sort_index(inplace=True)

        return data

    def _validate(self, columns: List[str], verifyNonStandardOctave):
        """
        Ensure that the provided data has only the standard

        Parameters
        ----------
        columns : List of strs
            List of NVSPL DataFrame columns.

        Raises
        ------
        AssertionError if any standard column is missing or if any non-standard and non-octave column is present.
        """
        # Verify that all NVSPL standard columns exist.
        missing_standard_cols = self.standard_fields - set(columns)
        assert missing_standard_cols == set(
        ), f"Missing the following standard NVSPL columns: {missing_standard_cols}"

        # Verify all non-standard columns are octave columns. Use verifyNonStandardOctave=False to allow extra columns
        if verifyNonStandardOctave:
            only_standard_cols = all(re.match(self.octave_regex, col)
                                     for col in (set(columns) - self.standard_fields))
            assert only_standard_cols is True, "NVSPL data contains unexpected NVSPL columns."


class Ais(gpd.GeoDataFrame):
    """
    A geopandas GeoDataFrame wrapper class to ensure consistent AIS data.

    Parameters
    ----------
    filepaths_or_data : List, str, or gpd.GeoDataFrame
        A directory containing AIS 

    """

    def __init__(self, filepaths_or_data: Union[List[str], str, gpd.GeoDataFrame]):
        data = self._read(filepaths_or_data)
        super().__init__(data=data)

    def parseAis(self, aisFileEntry, state=(None, None, 1)):

        timestamps, columns, index_index = state

        df = pd.read_csv(str(aisFileEntry),
                         engine='c',
                         usecols=columns,
                         low_memory=False
                         )  # read the .csv

        if 'mmsi' in df.columns:  # Then this is the more modern version of AIS file...

            # We must rename the headers to match
            df.columns = ['Base station time stamp', 'MMSI', 'callsign', 'IMO', 'Ship name',
                          'Navigational status (text)', 'Latitude', 'Longitude', 'Course over ground',
                          'Speed over ground', 'Destination', 'eta', 'Type of ship (text)', 'Draught',
                          'length', 'width']

            # We drop the fields that don't exist in the legacy files
            df.drop(['callsign', 'eta', 'length', 'width'],
                    axis=1, inplace=True, errors="ignore")

        # if there are 1090 MHz jet ADS-B points mixed into this dataset
        # this is a convenient way to make sure they are removed
        # we strictly require a 9-digit MMSI code
        df = df.loc[df["MMSI"] >= 100000000, :].copy()

        # tidy up all the header field names
        mask = df.iloc[:, 0].isin(['Base station time stamp'])
        df = df[~mask]
        header_list = ['Base station time stamp']
        import_header = df.axes[1]
        result = any(elem in import_header for elem in header_list)
        if result:
            pass
        else:
            raise KeyError

        # Standardize key field names and remove extra columns collected by the AIS logger
        if 'Base station time stamp' in df.columns:
            df = df.rename(columns={'Base station time stamp': "TIME"})

        df.drop(['IMO', 'Ship name', 'Type of ship (text)', 'Size A',
                'Size B', 'Size C', 'Size D', 'Draught', 'Destination', 'Heading',
                 'Navigational status (text)', 'Country (AIS)', 'Target class (text)',
                 'Data source type (text)', 'Data source region'], axis=1, inplace=True, errors="ignore")

        if 'Course over ground' in df.columns:
            df = df.rename(columns={'Course over ground': "heading"})

        if 'Latitude' in df.columns:
            df = df.rename(columns={'Latitude': "lat"})

        if 'Longitude' in df.columns:
            df = df.rename(columns={'Longitude': "lon"})

        if 'Speed over ground' in df.columns:
            df = df.rename(columns={'Speed over ground': "velocity"})

        # Delete duplicate records
        df.drop_duplicates(inplace=True)
        df.dropna(how="any", axis=0, inplace=True)

        # it is possible that upon removing ADS-B no points remain,
        # in which case we're done... we return an empty `pd.DataFrame` with formatted field names
        if (len(df) == 0):
            return df

        else:
            # For now, we assume the vessel's z-position is "at sea level"
            # a slower, but more accurate z-position would be derived
            # using the NOAA CO-OPS Data Retrieval API
            # https://api.tidesandcurrents.noaa.gov/api/prod/
            df["altitude"] = 0.0  # meters

            # The MXAK has released data with many different time formats
            # we can safely assume that a single fine has a consistent time format over all rows
            # our best approach is to check the first row against the known patterns:
            test = df["TIME"].iloc[0]
            time_pattern = test.split(" ")

            # Files exist with the following patterns:
            MXAK_timestamp_patterns = {1: "%d %b %Y %H:%M:%S UTC",
                                       2: "%Y-%m-%d %H:%M:%S UTC",
                                       3: "%m-%d-%Y %H:%M:%S UTC",
                                       4: "%Y-%m-%d %H:%M:%S AKST",
                                       5: "%Y-%m-%d %H:%M:%S AKDT",
                                       6: "%Y-%m-%d %H:%M:%S"}

            # Begin conditional timestamp formatting...
            if ((len(time_pattern) == 5) & (time_pattern[-1] == "UTC")):
                df["TIME"] = pd.to_datetime(
                    df["TIME"], format=MXAK_timestamp_patterns[1])

            elif ((len(time_pattern) == 3) & (time_pattern[-1] == "UTC")):
                try:
                    df["TIME"] = pd.to_datetime(
                        df["TIME"], format=MXAK_timestamp_patterns[2])

                except ValueError:
                    df["TIME"] = pd.to_datetime(
                        df["TIME"], format=MXAK_timestamp_patterns[3])

            elif ((len(time_pattern) == 3) & (time_pattern[-1] == "AKST")):
                try:
                    df["TIME"] = pd.to_datetime(
                        df["TIME"], format=MXAK_timestamp_patterns[4]) + dt.timedelta(hours=9)
                except ValueError:
                    # we must handle the day in November where we change from AKST to AKDT
                    df.loc[df["TIME"].str[-4:] == "AKST", "TIME"] = pd.to_datetime(df["TIME"].loc[df["TIME"].str[-4:] == "AKST"],
                                                                                   format=MXAK_timestamp_patterns[4]) + dt.timedelta(hours=9)
                    df.loc[df["TIME"].str[-4:] == "AKDT", "TIME"] = pd.to_datetime(df["TIME"].loc[df["TIME"].str[-4:] == "AKDT"],
                                                                                   format=MXAK_timestamp_patterns[5]) + dt.timedelta(hours=8)
                    # apparently we still need to nudge these into datetime format?
                    df["TIME"] = pd.to_datetime(df["TIME"])

            elif ((len(time_pattern) == 3) & (time_pattern[-1] == "AKDT")):
                try:
                    df["TIME"] = pd.to_datetime(
                        df["TIME"], format=MXAK_timestamp_patterns[5]) + dt.timedelta(hours=8)
                except ValueError:
                    # we must handle the day in March where we change from AKDT to AKST
                    df.loc[df["TIME"].str[-4:] == "AKST", "TIME"] = pd.to_datetime(df["TIME"].loc[df["TIME"].str[-4:] == "AKST"],
                                                                                   format=MXAK_timestamp_patterns[4]) + dt.timedelta(hours=9)
                    df.loc[df["TIME"].str[-4:] == "AKDT", "TIME"] = pd.to_datetime(df["TIME"].loc[df["TIME"].str[-4:] == "AKDT"],
                                                                                   format=MXAK_timestamp_patterns[5]) + dt.timedelta(hours=8)
                    df["TIME"] = pd.to_datetime(df["TIME"])

            elif ((len(time_pattern) == 2)):
                df["TIME"] = pd.to_datetime(
                    df["TIME"], format=MXAK_timestamp_patterns[6])

            else:
                raise ValueError(
                    "The file's timestamp format is not recognized!")

            df["TIME"].dt.tz_localize(tz="UTC")
            df["DATE"] = df["TIME"].dt.strftime("%Y%m%d")

            # Sort records by MMSI and TIME then reset dataframe index
            df.sort_values(["MMSI", "TIME"], inplace=True, ignore_index=True)

            # Calculate time difference between sequential waypoints for each watercraft
            df["TIME"] = pd.to_datetime(arg=df["TIME"], errors="coerce")
            df["dur_secs"] = df.groupby(
                "MMSI")["TIME"].diff().dt.total_seconds()
            df["dur_secs"] = df["dur_secs"].fillna(0)

            # Drop any identical waypoints in a single input file based on MMSI, time, lat, and lon
            df.drop_duplicates(
                subset=['MMSI', 'TIME', 'lat', 'lon'], keep='last')

            # Use threshold waypoint duration value to identify separate flights by a vessel then sum the number of "true" conditions to assign unique ID's
            df['diff_event'] = df['dur_secs'] >= 1200  # ( = 20 minutes)
            df['cumsum'] = df.groupby('MMSI')['diff_event'].cumsum()
            df['event_id'] = df['MMSI'].astype(
                'str') + "_" + df['cumsum'].astype(str) + "_" + df['DATE'].astype(str)

            # Let us only consider events with more than 2 AIS points
            df = df[df.groupby("event_id").event_id.transform(len) > 2]

            return df

    def _read(self, filepaths_or_data: Union[List[str], str, gpd.GeoDataFrame]):
        """
        Read in AIS points as formatted by the Alaska Marine Exchange (www.mxak.org).

        Parameters
        ----------
        filepaths_or_data : List, str, or gpd.GeoDataFrame
            A directory containing AIS files, a list of AIS files, or an existing gpd.GeoDataFrame of AIS data.

        Raises
        ------
        AssertionError if directory path or file path does not exists or is of the wrong format.
        """
        if isinstance(filepaths_or_data, gpd.GeoDataFrame):
            data = filepaths_or_data.to_crs("epsg:4326")

        else:
            if isinstance(filepaths_or_data, str):
                assert os.path.isdir(
                    filepaths_or_data), f"{filepaths_or_data} does not exist."
                filepaths_or_data = glob.glob(f"{filepaths_or_data}/*.csv")

            else:
                for file in filepaths_or_data:
                    assert os.path.isfile(file), f"{file} does not exist."
                    assert file.endswith(
                        '.csv'), f"Only .csv AIS files accepted."

            with concurrent.futures.ThreadPoolExecutor() as pool:
                parts = list(tqdm(pool.map(self.parseAis, filepaths_or_data), total=len(
                    filepaths_or_data), unit="AIS files"))

            data = pd.concat(parts)
            data = gpd.GeoDataFrame(
                data,
                geometry=gpd.points_from_xy(data["lon"], data["lat"]),
                crs="epsg:4326"
            )

        return data


class Adsb(gpd.GeoDataFrame):
    """
    A geopandas GeoDataFrame wrapper class to ensure consistent ADS-B data.

    Parameters
    ----------
    filepaths_or_data : List, str, or gpd.GeoDataFrame
        A directory containing ADS-B TSV files, a list of ADS-B TSV files, or an existing gpd.GeoDataFrame of ADS-B data.
    region : gpd.GeoDataFrame, default None
        A geodataframe containing the spatial region of interest. The associated geometry should be a polygon or multipolygon.
        ADS-B points inside this region will be loaded, and points outside will not. If None, all ADS-B data will be loaded.

    Raises
    ------
    AssertionError if directory path or file path does not exists or is of the wrong format, or if region is of the wrong type.
    """

    lat_lon_grid_resolution = 0.01
    index_file_name = "index.txt"

    def __init__(self, filepaths_or_data: Union[List[str], str, gpd.GeoDataFrame], region: gpd.GeoDataFrame = None):
        if region is not None:
            assert isinstance(region, gpd.GeoDataFrame), "Region is not a GeoDataFrame"
            assert region.geometry.geom_type.isin(["Polygon", "MultiPolygon"]).all(), "Region geometry must be Polygon or MultiPolygon"
            region = region.to_crs("epsg:4326")

        if isinstance(filepaths_or_data, gpd.GeoDataFrame):
            data = filepaths_or_data.to_crs("epsg:4326")

        else:
            # convert filepath list or directory to just a filepath list
            filepaths = []
            if isinstance(filepaths_or_data, str):
                assert os.path.isdir(
                    filepaths_or_data), f"{filepaths_or_data} does not exist."
                filepaths = glob.glob(f"{filepaths_or_data}/*.TSV")
            else:
                for file in filepaths_or_data:
                    assert os.path.isfile(file), f"{file} does not exist."
                    assert (file.endswith('.txt') | file.endswith(
                        '.TSV')), f"Only .TSV ADS-B files accepted."
                    filepaths.append(file)

            # read files
            data = self._read(filepaths, region)

        super().__init__(data=data)

    def _read(self, filepaths: List[str], region: gpd.GeoDataFrame = None):
        """
        Read in ADS-B points as formatted by NPS data loggers.

        Parameters
        ----------
        filepaths_or_data : List[str]
            A list of ADS-B files.
        region : gpd.GeoDataFrame
            The spatial region of interest

        Returns
        -------
        data: gpd.GeoDataFrame
            A GeoDataFrame containing the ADSB data.
        """

        t_file_load = 0
        t_process = 0
        t_final = 0

        # organize files by directory, since we maintain a separate index file for each directory
        dirs = {}
        for path in filepaths:
            dir = os.path.dirname(path)
            if dir not in dirs:
                dirs[dir] = []
            dirs[dir].append(path)

        dataframes = []
        pbar = tqdm(total=len(filepaths), desc='Loading ADS-B files',
                    unit='files', colour='green')
        for dir in dirs:
            # attempt to load that directory's index file (which may or may not exist)
            index, ranges = self._load_index(dir, region)

            index_updated = False
            for filepath in dirs[dir]:
                basename = os.path.basename(filepath)
                t0 = time.perf_counter()
                if basename in ranges:
                    # use the index to speed up file reading
                    df = self._read_tsv_ranges(filepath, ranges[basename])
                else:
                    # up-to-date index data doesn't exist, so read entire file and update the index
                    df = self._read_tsv_and_update_index(filepath, index)
                    index_updated = True
                
                t_file_load += (time.perf_counter() - t0)

                t0 = time.perf_counter()

                if len(df) > 0:
                    # convert to geodataframe and clip to spatial region before doing data processing
                    # this way the data processing (which is slow) has a lot less to do
                    df["lat"] = df["lat"].astype(int) / 1e7
                    df["lon"] = df["lon"].astype(int) / 1e7
                    gdf = gpd.GeoDataFrame(
                        data=df,
                        geometry=gpd.points_from_xy(df["lon"], df["lat"]),
                        crs="epsg:4326"
                    )
                    if region is not None:
                        gdf = gpd.sjoin(gdf, region, predicate="within", how="inner")
                    # Resetting the index is crucial for small spatial regions so that excess data doesn't get dropped.
                    # No idea why though
                    gdf.reset_index(drop=True, inplace=True)

                    gdf = self._process_raw_dataframe(gdf)
                    if gdf is not None:
                        dataframes.append(gdf)

                t_process += (time.perf_counter() - t0)

                pbar.update(1)

            if index_updated:
                self._save_index(dir, index)

        pbar.close()

        t0 = time.perf_counter()

        if len(dataframes) == 0:
            data = gpd.GeoDataFrame(geometry=[])
        else:
            data = pd.concat(dataframes, ignore_index=True)
            data.drop_duplicates(subset=['TIME', 'ICAO_address'], inplace=True, keep='last')

        t_final = (time.perf_counter() - t0)

        print(f"time loading files: {t_file_load:.3f} s")
        print(f"time processing: {t_process:.3f} s")
        print(f"time finishing up: {t_final:.3f} s")

        return data

    def _load_index(self, directory: str, region: gpd.GeoDataFrame = None):
        """Given a directory containing ADSB TSV files and their associated index file,
        attempts to load the parts of the index that correspond to the spatial region of interest.

        Parameters
        ----------
        directory: str
            The directory containing .TSV ADSB files and their associated index file, if it exists.
        region: gpd.GeoDataFrame, default None
            A polygon representing the spatial region of interest. All ADSB points inside this region will be loaded,
            and some points outside of the region may also be loaded.

        Returns
        -------
        index: dict
            A dictionary describing where to find the ADSB entries that occur in each spatial grid cell.
            If no index file exists or the index is invalid, this will be an empty dictionary.
        ranges: dict
            A dictionary describing for each file, which byte ranges are relevant to the spatial query.
            Only files that have been previously indexed and haven't been changed since then will be included.
            If no index file exists or the index is invalid, this will be an empty dictionary.
        """

        index_path = os.path.join(directory, self.index_file_name)
        if not os.path.exists(index_path):
            return {}, {}

        with open(index_path, "r", encoding="utf-8") as f:
            # read index headers
            f.readline()
            grid_res = float(f.readline())  # grid resolution of the index
            f.readline()  # for human use only
            file_mtimes = json.loads(f.readline())  # what files are present and when they were last modified (mtime)
            f.readline()  # for human use only
            grid_cell_byte_offsets = json.loads(f.readline())  # what grid cells are present and their byte offset in the index file
            f.readline()  # for human use only
            index_start = f.tell()

            # make sure we are using the same grid resolution, if not, index is invalid and should be replaced
            if grid_res != self.lat_lon_grid_resolution:
                return {}, {}

            # check if any files were modified since they were last indexed
            # if so, will need to remove them from the index later
            out_of_date = []
            for file in file_mtimes:
                full_path = os.path.join(directory, file)
                up_to_date = os.path.exists(full_path) and os.path.getmtime(
                    full_path) == file_mtimes[file]
                if not up_to_date:
                    out_of_date.append(file)

            # figure out which grid cells contain a part of the spatial region
            if region is not None:
                region_cells = self._grid_cells_intersecting_region(region)
            else:
                region_cells = grid_cell_byte_offsets.keys() # all cells

            # determine which grid cells to read from the index file
            # if there are out of date files, we need to read everything so we can remove the out of date stuff
            if len(out_of_date) > 0:
                cells_to_read = grid_cell_byte_offsets.keys()  # all cells
            else:
                cells_to_read = region_cells

            # read the necessary grid cells
            index = {}
            for grid_cell in cells_to_read:
                if grid_cell not in grid_cell_byte_offsets:
                    continue
                offset = grid_cell_byte_offsets[grid_cell]
                f.seek(index_start + offset)
                index = index | json.loads(f.readline())

            # delete out of date index records and re-save the index
            if len(out_of_date) > 0:
                for file in out_of_date:
                    for grid_cell in index:
                        if file in index[grid_cell]:
                            del index[grid_cell][file]
                self._save_index(directory, index)
            
            # for each file, determine which ranges should be read
            # we may have read more cells into the index than needed, so iterate over region cells specifically for determining this
            ranges = {}
            # make sure each file gets a range, even if it has no records in the region of interest
            # to communicate that it was indexed before
            for file in file_mtimes:
                ranges[file] = []
            for grid_cell in region_cells:
                if grid_cell not in index:
                    continue
                for file in index[grid_cell]:
                    ranges[file] += index[grid_cell][file]
            # sort ranges by start offset, probably helps a bit with file-reading speed
            for file in ranges:
                ranges[file].sort(key=lambda x: x[0])

        return index, ranges

    def _save_index(self, directory: str, index: dict):
        """Updates/creates an ADSB index file in a directory containing ADSB TSV files.

        Parameters
        ----------
        directory: str
            A directory containing ADSB files
        index: dict
            An index dict containing information about the TSV files in the directory.
            Note that this can be partial information, not all TSV files must be represented.
            The index dict will be merged with existing index data in the index file if it exists.
        """
        index_path = os.path.join(directory, self.index_file_name)
        with open(index_path, "w", encoding="utf-8", newline="\n") as f:
            # note that we need newline="\n" so that windows carriage returns aren't added,
            # which would otherwise mess with computing byte offsets

            # write header recording the grid resolution
            f.write("grid resolution (degrees)\n")
            f.write(f"{self.lat_lon_grid_resolution}\n")

            # write header that describes which files are present in the index and when they were last modified
            files_present = set()
            for grid_cell in index:
                for file in index[grid_cell]:
                    files_present.add(file)
            mtime_header = {}
            for file in files_present:
                full_path = os.path.join(directory, file)
                mtime = os.path.getmtime(full_path)
                mtime_header[file] = mtime
            f.write("last modified\n")
            json.dump(mtime_header, f)
            f.write("\n")

            # prepare lines of the index for writing
            index_lines = []
            byte_offsets = {}
            current_byte_offset = 0
            for grid_cell in index:
                line = json.dumps({grid_cell: index[grid_cell]}) + "\n"
                index_lines.append(line)
                byte_offsets[grid_cell] = current_byte_offset
                current_byte_offset += len(line.encode("utf-8"))

            # write header describing which grid cells are present and their byte offset in the index file
            f.write("grid cell byte offsets in this file\n")
            json.dump(byte_offsets, f)
            f.write("\n")

            # write the index
            f.write("index\n")
            f.writelines(index_lines)

    def _add_range_to_index(self, index, grid_cell, filepath, start, length):
        """Utility function for inserting items into an index."""
        file = os.path.basename(filepath)
        if not grid_cell in index:
            index[grid_cell] = {}
        if not file in index[grid_cell]:
            index[grid_cell][file] = []
        index[grid_cell][file].append((start, length))

    def _read_tsv_and_update_index(self, filepath, index):
        """Reads a TSV file and updates the index."""

        tqdm.write(f"Reading full file {filepath}")
        with open(filepath, "r", encoding="utf-8-sig") as f:
            header_line = f.readline()
            fieldnames = next(csv.reader([header_line], delimiter="\t"))
            rows = []
            prev_grid_cell = None
            start_offset = f.tell()
            while True:
                offset = f.tell()
                line = f.readline()
                if not line:
                    break
                row = next(csv.DictReader([line], fieldnames, delimiter="\t"))
                # check for a logging blip causing two rows to get collapsed, resulting in extra values in a row
                # these extra values get collected into a list referenced by the key None
                # also check for extra header inserted by the logger
                if (None in row) or (row["timestamp"] == "timestamp"):
                    # don't include the broken line in the index, end the last range and start a new one
                    if prev_grid_cell is not None:
                        self._add_range_to_index(
                            index, prev_grid_cell, filepath, start_offset, offset-start_offset)
                    prev_grid_cell = None
                    start_offset = f.tell()
                    continue
                rows.append(row)

                # if grid cell changed, write the previous grid cell's byte range to the index
                lat, lon = int(row["lat"]) / 1e7, int(row["lon"]) / 1e7
                grid_cell = self._get_grid_cell(lat, lon)
                if prev_grid_cell is not None and prev_grid_cell != grid_cell:
                    self._add_range_to_index(
                        index, prev_grid_cell, filepath, start_offset, offset-start_offset)
                    start_offset = offset
                prev_grid_cell = grid_cell

            # add the final range
            self._add_range_to_index(
                index, grid_cell, filepath, start_offset, offset-start_offset)

        return pd.DataFrame(rows).convert_dtypes()

    def _read_tsv_ranges(self, filepath, ranges):
        """Reads sections of a TSV file specified by the `ranges` parameter."""

        with open(filepath, "r", encoding="utf-8-sig") as f:
            header_line = f.readline()
            fieldnames = next(csv.reader([header_line], delimiter="\t"))
            lines = []
            for (start, length) in ranges:
                f.seek(start)
                while f.tell() < start + length:
                    line = f.readline()
                    if not line:
                        break
                    lines.append(line)
        
        rows = csv.DictReader(lines, fieldnames, delimiter="\t")
        return pd.DataFrame(rows).convert_dtypes()

    def _get_grid_cell(self, lat, lon):
        """Determine grid cell for a certain coordinate using Decimal library to avoid annoying floating point precision problems"""
        res = Decimal(str(self.lat_lon_grid_resolution))
        lat = (Decimal(str(lat)) / res).quantize(0, ROUND_FLOOR) * res
        lon = (Decimal(str(lon)) / res).quantize(0, ROUND_FLOOR) * res
        return f"{lat},{lon}"
    
    def _grid_cells_intersecting_region(self, region: gpd.GeoDataFrame, visualize=False):
        """Determine which grid cells intersect a spatial region."""

        intersecting_cells = []
        polys = []
        res = self.lat_lon_grid_resolution

        # simplify geometry to make computation tractable
        # also combine geometry in case multiple polygons were provided
        geom = region.simplify(0.1 * res).union_all()

        # Get the bounding box of the geometry
        minx, miny, maxx, maxy = geom.bounds
        left_cell_x = np.floor(minx / res) * res
        top_cell_y = np.floor(miny / res) * res
        
        # Generate all possible grid cells and test intersection
        for x in np.arange(left_cell_x, maxx, res):
            for y in np.arange(top_cell_y, maxy, res):
                cell = box(x, y, x + res, y + res)
                if cell.intersects(geom):
                    # use center point for getting cell name to avoid floating point rounding issues
                    cell_name = self._get_grid_cell(y + 0.5*res, x + 0.5*res)
                    intersecting_cells.append(cell_name)
                    polys.append(cell)
        
        if visualize:
            fig, ax = plt.subplots(figsize=(8,8))
            grid_gdf = gpd.GeoDataFrame(geometry=polys, crs=region.crs)
            region.plot(ax=ax, color='lightblue', edgecolor='blue', label='Original Geometry')
            grid_gdf.plot(ax=ax, facecolor='none', edgecolor='red', linewidth=1, label='Grid Cells')
            plt.show()
        
        return intersecting_cells

    def _process_raw_dataframe(self, df):
        """Processes a raw dataframe read from the TSV file, and does some data cleaning.

        Parameters
        ----------
        df: pd.DataFrame
            The dataframe read from an ADSB TSV file.
        
        Returns
        -------
        df: pd.DataFrame or None
            The cleaned dataframe, or None if the cleaning removed all dataframe rows
        """

        # remove extra header rows inserted by the ADSB logger
        mask = df.iloc[:, 0].isin(["TIME", "timestamp"])
        df = df[~mask]
        if len(df) == 0:
            return None

        # verify a time header exists
        header_list = ["TIME", "timestamp"]
        import_header = df.axes[1]
        result = any(elem in import_header for elem in header_list)
        if not result:
            raise KeyError

        # Standardize key field names and remove extra columns collected by the ADS-B df logger
        if "timestamp" in df.columns:
            df = df.rename(columns={"timestamp": "TIME"})
        if "valid_flags" in df.columns:
            df = df.rename(columns={"valid_flags": "validFlags"})
        df.drop(["squawk", "altitude_type", "alt_type", "altType", "callsign",
                 "emitter_type", "emitterType"], axis=1, inplace=True, errors="ignore")

        # Delete duplicate and NA records
        df.drop_duplicates(inplace=True)
        df.dropna(how="any", axis=0, inplace=True)
        if len(df) == 0:
            return None

        # Unpack validFLags and convert the 2-byte flag field into a list of Boolean values
        flags_names = ["valid_BARO", "valid_VERTICAL_VELOCITY", "SIMULATED_REPORT", "valid_IDENT",
                       "valid_CALLSIGN", "valid_VELOCITY", "valid_HEADING", "valid_ALTITUDE", "valid_LATLON"]
        flags = df["validFlags"].apply(
            lambda t: list(bin(int(t, 16))[2:].zfill(9)[-9:]))
        flags_df = pd.DataFrame(list(flags), columns=flags_names).replace(
            {'0': False, '1': True}).infer_objects(copy=False)
        df = pd.concat(
            [df.drop("validFlags", axis=1), flags_df], axis=1)

        # Keep only those records with valid latlon and altitude values based on validFlags
        df.dropna(how="any", axis=0, inplace=True)
        # if df["valid_LATLON"].sum() == len(df.index):
        #     invalidLatLon = 0
        # else:
        #     invalidLatLon = round(
        #         100 - df["valid_LATLON"].sum() / len(df.index) * 100, 2)
        # if df["valid_ALTITUDE"].sum() == len(df.index):
        #     invalidAltitude = 0
        # else:
        #     invalidAltitude = round(
        #         100 - df["valid_ALTITUDE"].sum() / len(df.index) * 100, 2)
        df.drop(df[df["valid_LATLON"] == "False"].index, inplace=True)
        df.drop(df[df["valid_ALTITUDE"] ==
                "False"].index, inplace=True)
        if len(df) == 0:
            return None
        
        # Ensure remaining field values except TIME are in proper numeric format
        df.replace('-', np.nan, inplace=True)
        df["ICAO_address"] = df["ICAO_address"].astype(str)
        # df["lat"] = df["lat"].astype(int)
        # df["lon"] = df["lon"].astype(int)
        df["altitude"] = df["altitude"].astype(int)
        df["heading"] = df["heading"].astype(int)
        df["hor_velocity"] = df["hor_velocity"].astype(int)
        df["ver_velocity"] = df["ver_velocity"].astype(int)
        df["tslc"] = df["tslc"].astype(int)

        # Convert Unix timestamp to datetime objects in UTC and re-scale selected variable values
        df["TIME"] = pd.to_datetime(df["TIME"].astype(int), unit="s")
        df["DATE"] = df["TIME"].dt.strftime("%Y%m%d")
        # df["lat"] = df["lat"] / 1e7
        # df["lon"] = df["lon"] / 1e7
        df["altitude"] = df["altitude"] / 1e3
        df["heading"] = df["heading"] / 1e2
        df["hor_velocity"] = df["hor_velocity"] / 1e2
        df["ver_velocity"] = df["ver_velocity"] / 1e2

        # Keep only those records with TSLC values of 1 or 2 seconds
        # invalidTslc = len(
        #     df.query("tslc >= 3 or tslc == 0")) / df.shape[0] * 100
        df.drop(df[df["tslc"] >= 3].index, inplace=True)
        df.drop(df[df["tslc"] == 0].index, inplace=True)
        if len(df) == 0:
            return None

        # Keep only those records with realistic altitudes
        # 10000 meters = 32808 feet; this should encompass most flights
        # NOTE: some jet aircraft may be eliminated by this process
        df = df.loc[(df["altitude"] > 0) & (
            df["altitude"] <= 10000), :]
        if len(df) == 0:
            return None

        # Sort records by ICAO Address and TIME then reset dfframe index
        df.sort_values(["ICAO_address", "TIME"],
                       inplace=True, ignore_index=True)

        # Calculate time difference between sequential waypoints for each aircraft
        df["dur_secs"] = df.groupby("ICAO_address")[
            "TIME"].diff().dt.total_seconds()
        df["dur_secs"] = df["dur_secs"].fillna(0)

        # Count then delete any identical waypoints in a single input file based on ICAO_address, time, lat, and lon
        # duplicateWaypoints = 100 - \
        #     (len(df.drop_duplicates(
        #         subset=['ICAO_address', 'TIME', 'lat', 'lon'])) / len(df) * 100)
        df.drop_duplicates(
            subset=['ICAO_address', 'TIME', 'lat', 'lon'], keep='last')

        # Use threshold waypoint duration value to identify separate flights by an aircraft then sum the number of "true" conditions to assign unique ID's
        df['diff_flight'] = df['dur_secs'] >= 900
        df['cumsum'] = df.groupby('ICAO_address')[
            'diff_flight'].cumsum()
        df['flight_id'] = df['ICAO_address'] + "_" + \
            df['cumsum'].astype(str) + "_" + df['DATE']

        # Remove records where there is only one recorded waypoint for an aircraft and fields that are no longer needed
        df = df[df.groupby("flight_id").flight_id.transform(len) > 1]
        df = df.drop(columns=['tslc', 'dur_secs', 'diff_flight', 'cumsum', 'valid_BARO', 'valid_VERTICAL_VELOCITY', 'SIMULATED_REPORT',
                              'valid_IDENT', 'valid_CALLSIGN', 'valid_VELOCITY', 'valid_HEADING', 'valid_ALTITUDE', 'valid_LATLON', 'DATE'])

        return df


class EarlyAdsb(gpd.GeoDataFrame):
    """
    A geopandas GeoDataFrame wrapper class to ensure consistent ADS-B data.

    Parameters
    ----------
    filepaths_or_data : List, str, or gpd.GeoDataFrame
        A directory containing ADS-B TSV files, a list of ADS-B TSV files, or an existing gpd.GeoDataFrame of ADS-B data.
    """

    def __init__(self, filepaths_or_data: Union[List[str], str, gpd.GeoDataFrame]):
        data = self._read(filepaths_or_data)
        data.drop_duplicates(subset=['TIME'], inplace=True, keep='last')
        super().__init__(data=data)

    def _read(self, filepaths_or_data: Union[List[str], str, gpd.GeoDataFrame]):
        """
        Read in ADS-B points as formatted by early-development NPS data loggers (circa 2019).

        Parameters
        ----------
        filepaths_or_data : List, str, or gpd.GeoDataFrame
            A directory containing ADS-B files, a list of ADS-B files, or an existing gpd.GeoDataFrame of ADS-B data.

        Raises
        ------
        AssertionError if directory path or file path does not exists or is of the wrong format.
        """
        if isinstance(filepaths_or_data, gpd.GeoDataFrame):
            data = filepaths_or_data.to_crs("epsg:4326")

        else:
            if isinstance(filepaths_or_data, str):
                assert os.path.isdir(
                    filepaths_or_data), f"{filepaths_or_data} does not exist."
                filepaths_or_data = glob.glob(f"{filepaths_or_data}/*.TSV")

            else:
                for file in filepaths_or_data:
                    assert os.path.isfile(file), f"{file} does not exist."
                    assert (file.endswith('.txt')
                            ), f"Only .txt ADS-B files accepted."

            data = pd.DataFrame()
            for file in tqdm(filepaths_or_data, desc='Loading ADS-B files', unit='files', colour='green'):
                df = pd.read_csv(file, sep="\t")

                df.columns = ["ICAO_address", "TIME", "lat", "lon", "altitude"]
                df["TIME"] = df["TIME"].apply(lambda t: dt.datetime.strptime(
                    t, "%Y/%m/%d %H:%M:%S.%f").replace(microsecond=0))
                df["DATE"] = df["TIME"].dt.strftime("%Y%m%d")

                # unlike later loggers, EarlyAdsb was collected in feet MSL
                # we need to convert altitude from feet to meters!
                df["altitude"] = 0.3048*df["altitude"]

                # Keep only those records with realistic altitudes
                # 10000 meters = 32808 feet; this should encompass most flights
                # NOTE: some jet aircraft may be eliminated by this process
                df = df.loc[(df["altitude"] > 0) & (
                    df["altitude"] <= 10000), :]

                # Sort records by ICAO Address and TIME then reset dataframe index
                df.sort_values(["ICAO_address", "TIME"],
                               inplace=True, ignore_index=True)

                # Calculate time difference between sequential waypoints for each aircraft
                df["dur_secs"] = df.groupby("ICAO_address")[
                    "TIME"].diff().dt.total_seconds()
                df["dur_secs"] = df["dur_secs"].fillna(0)

                # Use threshold waypoint duration value to identify separate flights by an aircraft
                # then sum the number of "true" conditions to assign unique ID's
                df['diff_flight'] = df['dur_secs'] >= 900
                df['cumsum'] = df.groupby('ICAO_address')[
                    'diff_flight'].cumsum()
                df['flight_id'] = df['ICAO_address'] + "_" + \
                    df['cumsum'].astype(str) + "_" + df['DATE']

                # Remove records where there is only one recorded waypoint for an aircraft
                df = df[df.groupby("flight_id").flight_id.transform(len) > 1]

                data = pd.concat([data, df], ignore_index=True)

            data = gpd.GeoDataFrame(
                data,
                geometry=gpd.points_from_xy(data["lon"], data["lat"]),
                crs="epsg:4326"
            )

        return data


class FAAReleasable():
    """Use a pre-downloaded copy of the U.S. Federal Aviation Administration's releasable aircraft database
    (https://www.faa.gov/licenses_certificates/aircraft_certification/aircraft_registry/releasable_aircraft_download)
    to glean various properties associated with a set of aircraft tracks.
    """

    def __init__(self, FAA_path, aircraft_corrections_path=None, n_numbers=None, icao_addresses=None, warnings=True):
        """Load FAA data into a DataFrame for a subset of aircraft (or all aircraft).
        The first time this is run, creates an index file that lives in the same directory as the FAA path to speed up future database reading.

        After initializing, use the .data attribute to access the loaded dataframe.

        Parameters
        ----------
        FAA_path: str
            Path to the MASTER.txt file downloaded from the FAA
        aircraft_corrections_path: str, default None
            Path to a corrections text file for fixing errors with aircraft type in the FAA database. Formatted as a JSON object where keys are N-numbers and values are the correct aircraft type (e.g. "Fixed-wing")
        n_numbers: array_like of str, default None
            If provided, only load data for these N-numbers.
        icao_addresses: array_like of str, default None
            If provided, only load data for these ICAO 24-bit addresses.
        warnings: bool, default True
            If True, warns the user if aircraft are not found in the FAA database.

        Notes
        -----
        [1] If both `n_numbers` and `icao_addresses` are passed, this class will load aircraft with a N-Number or a matching ICAO address. If a provided N-number and ICAO address refer to the same aircraft, the loaded dataframe will have duplicate rows.
        """

        self.FAA_path = FAA_path
        self.aircraft_corrections_path = aircraft_corrections_path
        self.n_numbers = n_numbers
        self.icao_addresses = icao_addresses
        self.warnings = warnings

        self.index_path = os.path.join(
            os.path.dirname(FAA_path), "MASTER_index.json")

        # If we only want some of the data and we have an index, use it. Otherwise load all data
        index_loaded = self._load_index()
        if index_loaded and (n_numbers is not None or icao_addresses is not None):
            self._read_using_index()
        else:
            self._read_and_build_index()

        # Convert aircraft type to human-readable format
        Type_Map = {"4": "Fixed-wing", "5": "Jet", "6": "Helicopter"}
        self.data['TYPE AIRCRAFT'] = self.data['TYPE AIRCRAFT'].apply(
            lambda l: Type_Map[str(l)] if str(l) in Type_Map else l)

        self._apply_corrections()

    def _load_index(self):
        """Attempt to load the index from disk.

        Returns
        -------
        success: bool
            True if an index exists, is up to date, and has been loaded successfully into self.index. False otherwise.
        """
        if not os.path.exists(self.index_path):
            return False
        with open(self.index_path, "r") as f:
            index = json.load(f)

        # check that database hasn't changed
        if "database_last_modified" not in index:
            return False
        if index["database_last_modified"] != os.path.getmtime(self.FAA_path):
            return False
        # no need to keep metadata in the index, might be confusing?
        del index["database_last_modified"]

        self.index = index
        return True

    def _estimate_line_count(self, filename, sample_size=1024 * 1024):
        """Use a 1MB sample to estimate the number of lines in a large file"""
        file_size = os.path.getsize(filename)
        with open(filename, 'rb') as f:
            sample = f.read(sample_size)
        newlines = sample.count(b'\n')
        if not newlines:
            return 0
        if file_size < sample_size:
            return newlines
        return int((file_size / sample_size) * newlines)

    def _read_and_build_index(self):
        """Read data, building an index along the way. Assigns data to self.data, and saves the index to a file."""
        self.index = {
            "n_number": {},
            "icao": {}
        }
        n_lines = self._estimate_line_count(self.FAA_path)
        pbar = tqdm(total=n_lines-1, desc="Loading all FAA Data")

        with open(self.FAA_path, 'r', encoding='utf-8-sig') as f:
            header_line = f.readline()
            rows = []
            while True:
                offset = f.tell()
                line = f.readline()
                if not line:
                    break
                row = next(csv.DictReader([header_line, line]))
                row = {k: v.strip() if v is not None else v for k,
                       v in row.items()}  # remove extra whitespace
                rows.append(row)
                self.index["n_number"][row["N-NUMBER"]] = offset
                self.index["icao"][row["MODE S CODE HEX"]] = offset
                pbar.update(1)
        pbar.close()

        self.data = pd.DataFrame(rows).convert_dtypes()

        # save index to file
        to_save = self.index.copy()
        to_save["database_last_modified"] = os.path.getmtime(self.FAA_path)
        with open(self.index_path, "w") as f:
            json.dump(to_save, f)
        print(
            f"Saved FAA database index to {os.path.abspath(self.index_path)}")

        # filter data
        self._check_aircraft_filter()
        if self.n_numbers is not None or self.icao_addresses is not None:
            selection = np.zeros(len(self.data), dtype=bool)
            if self.n_numbers is not None:
                selection = selection | self.data['N-NUMBER'].isin(
                    self.n_numbers)
            if self.icao_addresses is not None:
                selection = selection | self.data['MODE S CODE HEX'].isin(
                    self.icao_addresses)
            self.data = self.data[selection]

    def _read_using_index(self):
        assert self.index is not None and "n_number" in self.index and "icao" in self.index, \
            f"Something is wrong with the FAA index, please delete {os.path.abspath(self.index_path)} and try again"

        self._check_aircraft_filter()

        # load file offsets from the index
        offsets = []
        if self.n_numbers is not None:
            for n in self.n_numbers:
                offsets.append(self.index["n_number"][str(n)])
        if self.icao_addresses is not None:
            for code in self.icao_addresses:
                offsets.append(self.index["icao"][str(code)])

        offsets.sort()  # probably improves disk access speed

        with open(self.FAA_path, mode="r", encoding="utf-8-sig") as f:
            header_line = f.readline()
            rows = []
            for offset in offsets:
                f.seek(offset)
                row = next(csv.DictReader([header_line, f.readline()]))
                row = {k: v.strip() if v is not None else v for k,
                       v in row.items()}  # remove extra whitespace
                rows.append(row)

        self.data = pd.DataFrame(rows).convert_dtypes()

    def _check_aircraft_filter(self):
        """Checks that self.n_numbers and self.icao_addresses are actually in the FAA database. If not, warns the user and removes that value."""

        if self.n_numbers is not None:
            found_n_numbers = []
            for n in self.n_numbers:
                if str(n) in self.index["n_number"]:
                    found_n_numbers.append(n)
                elif self.warnings:
                    warn(
                        f"N-number {n} not found in the FAA database, skipping")
            self.n_numbers = found_n_numbers

        if self.icao_addresses is not None:
            found_codes = []
            for code in self.icao_addresses:
                if str(code) in self.index["icao"]:
                    found_codes.append(code)
                elif self.warnings:
                    warn(
                        f"ICAO Address {code} not found in the FAA database, skipping")
            self.icao_addresses = found_codes

    def _apply_corrections(self):
        if self.aircraft_corrections_path is None:
            return

        # Open aircraft corrections from specified path
        with open(self.aircraft_corrections_path) as f:
            aircraft_corrections = json.load(f)

        # Correct the aircraft lookup table.
        for n_number in aircraft_corrections:
            affected_rows = self.data['N-NUMBER'] == n_number
            correct_type = aircraft_corrections[n_number]
            self.data.loc[affected_rows, 'TYPE AIRCRAFT'] = correct_type


class Tracks(gpd.GeoDataFrame):
    """
    A geopandas GeoDataFrame wrapper class to standardize track points.

    Parameters
    ----------
    data : gpd.GeoDataFrame
        A GeoDataFrame of track points.
    id_col : str
        The name of the column containing aa unique identifier to group track points by.
        This column will be given the standardized name of track_id and converted to a string.
            E.g. flight id, license plate
    datetime_col : str
        A tracks GeoDataFrame is required to have a column with the datetime of each track point.
        This column will be given the standardized name of "point_dt".
    z_col : str, default None
        A tracks GeoDataFrame can have a column with the altitude of the points.
        This column will be given the standardized name of "z".

    Notes
    -----
    Currently, there is a bug with GeoPandas where running to_crs() will delete the z values of Points as mentioned
    in this post https://stackoverflow.com/questions/72987452/geopands-to-crs-dropping-z-values. Therefore, z values must
    be kept in a separate standard column until this bug has been resolved.
    """

    def __init__(self, data: gpd.GeoDataFrame, id_col: str, datetime_col: str, z_col: Optional[str] = None):
        col_renames = {id_col: 'track_id', datetime_col: 'point_dt'}
        if z_col:
            col_renames[z_col] = 'z'
        data.rename(columns=col_renames, inplace=True)
        if 'geometry' not in data:
            data.rename_geometry('geometry', inplace=True)
        data['track_id'] = data.track_id.astype(str)
        data.sort_values(by=['track_id', 'point_dt'],
                         ascending=True, inplace=True)
        super().__init__(data=data)


class Annotations(gpd.GeoDataFrame):
    """
    A geopandas GeoDataFrame wrapper class to standardize track annotations.

    Parameters
    ----------
   filename : str, default None
       Filename to read annotation data from. If no filename is passed, an empty Annotations GeoDataFrame
       will be created.
    only_valid : bool, default False
        If True and an annotation filename was passed, only valid records will be loaded.
    """

    def __init__(self, filename: Optional[str] = None, only_valid: bool = False):

        if filename:
            data = gpd.read_file(filename).astype(
                {'start_dt': 'datetime64[ns]', 'end_dt': 'datetime64[ns]'})

            # Sometimes the annotation file is read in with the valid and audible columns as booleans and other times
            #  as objects depending on what values are stored.
            try:
                data.valid.replace({'1': True, '0': False}, inplace=True)
                data.audible.replace({'1': True, '0': False}, inplace=True)
            except TypeError:
                pass

            if only_valid:
                data = data[data.valid == True]

        else:

            data = gpd.GeoDataFrame(columns=['_id', 'start_dt', 'end_dt', 'valid', 'audible', 'geometry', 'note'],
                                    geometry='geometry', crs='epsg:4326')

        super().__init__(data=data, crs=data.crs)
