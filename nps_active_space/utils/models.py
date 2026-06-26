import time
from datetime import datetime
from types import GeneratorType
import concurrent.futures
from tqdm import tqdm
from pyproj import Transformer
import datetime as dt
import glob
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
from .computation import contiguous_regions
pd.options.mode.copy_on_write = True
pd.set_option('future.no_silent_downcasting', True)

__all__ = [
    'Microphone',
    'Nvspl',
    'Adsb',
    'EarlyAdsb',
    'FAAReleasable',
    'Tracks',
    'Annotations',
    'Srcid',
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
        'SiteID', 'dbA', 'dbC', 
        'Voltage', 'WindSpeed', 'WindDir', 'TempIns',
        'TempOut', 'Humidity', 'INVID', 'INSID',
        'GChar1', 'AdjustmentsApplied',
        'CalibrationAdjustment', 'GPSTimeAdjustment',
        'GainAdjustment', 'Status'
    } # to ensure compatability with PamGuide-based conversions the set of fields {'GChar2', 'GChar3', 'dbF'} are removed from consideration

    octave_regex = re.compile(r"^H[0-9]+$|^H[0-9]+p[0-9]$")

    def __init__(self, filepaths_or_data: Union[List[str], str, pd.DataFrame]):
        if isinstance(filepaths_or_data, list) and len(filepaths_or_data) == 0:
            raise ValueError("No NVSPL files found")
        data = self._read(filepaths_or_data)
        super().__init__(data=data)

    def parseNvspl(self, nvsplFileEntry, state=(None, None, 1)):
        timestamps, columns, index_index = state

        # occasionally NVSPL files can get corrupted when copying, etc.
        # when this happens, it's nice to know which one had the issue
        try:
            df = pd.read_csv(str(nvsplFileEntry),
                            engine='c',
                            parse_dates=True,
                            index_col=index_index,
                            usecols=columns
                            )
        except Exception as e:
            print(f" (!!!) Error occurred reading {nvsplFileEntry}")
            raise e

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
    start_dt : datetime.datetime, default None
        A naive datetime object representing the earliest time that should be loaded.
    end_dt : datetime.datetime, default None
        A naive datetime object representing the latest time that should be loaded.

    Notes
    -----
    [1] If both `start_dt` and `end_dt` are None, all data will be loaded (in the region of interest).

    [2] If `start_dt` is passed but `end_dt` is None, all points after `start_dt` will be returned, and the same for if only `end_dt` is passed.

    Raises
    ------
    AssertionError if directory path or file path does not exists or is of the wrong format, or if region is of the wrong type.
    """
    
    lat_lon_grid_resolution = 0.1
    index_file_name = "index.txt"

    def __init__(self, filepaths_or_data: Union[List[str], str, gpd.GeoDataFrame],
                 region: gpd.GeoDataFrame = None,
                 start_dt: datetime = None,
                 end_dt: datetime = None):
        
        if region is not None:
            assert isinstance(region, gpd.GeoDataFrame), "Region is not a GeoDataFrame"
            assert region.geometry.geom_type.isin(["Polygon", "MultiPolygon"]).all(), "Region geometry must be Polygon or MultiPolygon"
            region = region.to_crs("epsg:4326")
        
        if start_dt is not None and end_dt is not None:
            assert start_dt < end_dt, "Start time must be earlier than end time"
        
        if isinstance(filepaths_or_data, gpd.GeoDataFrame):
            data = filepaths_or_data.to_crs("epsg:4326")
        elif isinstance(filepaths_or_data, pd.DataFrame):
            # this is a weird quirk of subclassing a geodataframe, gets triggered sometimes when you try to print this object
            data = filepaths_or_data
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
            assert len(filepaths) > 0, "No .TSV files found to read, check the file path(s)"
            data = self._read(filepaths, region, start_dt, end_dt)

        super().__init__(data=data)

    def _read(self, filepaths: List[str], region: gpd.GeoDataFrame = None,
              start_dt: datetime = None, end_dt: datetime = None):
        """
        Read in ADS-B points as formatted by NPS data loggers.

        Parameters
        ----------
        filepaths_or_data : List[str]
            A list of ADS-B files.
        region : gpd.GeoDataFrame, default None
            The spatial region of interest
        start_dt : datetime.datetime, default None
        end_dt : datetime.datetime, default None

        Returns
        -------
        data: gpd.GeoDataFrame
            A GeoDataFrame containing the ADSB data.
        """

        # organize files by directory, since we maintain a separate index file for each directory
        dirs = {}
        for path in filepaths:
            dir = os.path.dirname(path)
            if dir not in dirs:
                dirs[dir] = []
            dirs[dir].append(path)

        dataframes = []
        for dir in dirs:
            print(f"Loading ADSB files in {dir}")
            dir_files = dirs[dir]

            # attempt to load that directory's index file (which may or may not exist)
            index, timestamp_ranges, ranges = self._load_index(dir, region, start_dt, end_dt)

            # for `ranges` and `region`, will pass the same object to all threads
            # this is okay because the threads are only reading these objects, not writing
            ranges_copy = [ranges for _ in range(len(dir_files))]
            region_copy = [region for _ in range(len(dir_files))]
            start_dt_copy = [start_dt for _ in range(len(dir_files))]
            end_dt_copy = [end_dt for _ in range(len(dir_files))]
            
            with concurrent.futures.ThreadPoolExecutor() as pool:
                results = list(tqdm(
                    pool.map(self._read_file, dir_files, ranges_copy, region_copy, start_dt_copy, end_dt_copy),
                    total=len(dir_files), unit="files"))

            # process results from all the threads
            index_updated = False
            for (df, index_update, timestamp_update) in results:

                if df is not None:
                    dataframes.append(df)

                if index_update is not None:
                    assert timestamp_update is not None, "Timestamp update should be part of the indexing process, what happened??"

                    # merge this index update with the rest of the index
                    for grid_cell in index_update:
                        if grid_cell not in index:
                            index[grid_cell] = {}
                        files = list(index_update[grid_cell].keys())
                        assert len(files) == 1, "Somehow the index update from one file contains multiple files??"
                        file = files[0]
                        # can overwrite index[grid_cell][file] if it exists b/c we should know everything about that file
                        index[grid_cell][file] = index_update[grid_cell][file]

                    # merge timestamp update
                    timestamp_ranges = timestamp_ranges | timestamp_update

                    index_updated = True

            if index_updated:
                self._save_index(dir, index, timestamp_ranges)

        if len(dataframes) == 0:
            data = gpd.GeoDataFrame(geometry=[])
        else:
            data = pd.concat(dataframes, ignore_index=True)
            data.drop_duplicates(subset=['TIME', 'ICAO_address'], inplace=True, keep='last')
            data = gpd.GeoDataFrame(
                data=data,
                geometry=gpd.points_from_xy(data["lon"], data["lat"]),
                crs="epsg:4326")

        return data

    def _read_file(self, filepath: str, ranges: dict, region: gpd.GeoDataFrame = None,
                   start_dt: datetime = None, end_dt: datetime = None):
        """Read a single ADSB .TSV file. Uses file byte ranges from the index to go faster if provided.
        If no ranges are provided, indexes the file while loading it and returns that part of the index.
        Clips to the region and datetime range, if provided.
        
        Returns
        -------
        tuple of (df, index_update, datetime_range)
            df: pd.DataFrame or None
                The read contents of the ADSB file. If nothing was read of adequate quality that matches the region, is None.
            index_update: dict or None
                If needed, a spatial index for this file. Will be combined with the index from other files to build an overall spatial index.
                If this file was previously indexed, is None.
            timestamp_range_update: tuple of (float, float) or None
                Range of unix timestamps present in the file. Is None if the file wasn't being indexed.
        """

        index_update = None
        timestamp_range_update = None

        basename = os.path.basename(filepath)
        if basename in ranges:
            # use the index to speed up file reading
            df = self._read_tsv_ranges(filepath, ranges[basename])
        else:
            # up-to-date index data doesn't exist, so read entire file and update the index
            index_update = {}
            df = self._read_tsv_and_update_index(filepath, index_update)
        
        if len(df) == 0:
            return None, index_update, timestamp_range_update
        
        timestamps = pd.to_numeric(df["TIME" if "TIME" in df else "timestamp"], errors="coerce")
        if index_update is not None:
            # should update the timestamp range too
            timestamp_range_update = {
                basename: (float(timestamps.min()), float(timestamps.max()))}  # json doesn't like numpy types somehow
        
        # clip to datetime range before doing processing
        if start_dt is not None:
            df = df[pd.Series(timestamps >= start_dt.timestamp(), index=df.index)]
        if end_dt is not None:
            df = df[pd.Series(timestamps <= end_dt.timestamp(), index=df.index)]
        
        # clip to spatial region before doing data processing
        # this way the data processing (which is slow) has a lot less to do
        if region is not None:
            # make a temporary geodataframe to figure out which rows are within the region
            lat = df["lat"].astype(int) / 1e7
            lon = df["lon"].astype(int) / 1e7
            gdf = gpd.GeoDataFrame(
                geometry=gpd.points_from_xy(lon, lat),
                crs="epsg:4326"
            )
            mask = gdf.within(region.union_all())
            mask.index = df.index
            df = df[mask]
            # Resetting the index is crucial for small spatial regions so that excess data doesn't get dropped.
            # No idea why though
            df.reset_index(drop=True, inplace=True)

        df = self._process_raw_dataframe(df)
            
        return df, index_update, timestamp_range_update

    def _load_index(self, directory: str, region: gpd.GeoDataFrame = None,
                    start_dt: datetime = None, end_dt: datetime = None):
        """Given a directory containing ADSB TSV files and their associated index file,
        attempts to load the parts of the index that correspond to the spatial region of interest.

        Parameters
        ----------
        directory: str
            The directory containing .TSV ADSB files and their associated index file, if it exists.
        region: gpd.GeoDataFrame, default None
            A polygon representing the spatial region of interest. All ADSB points inside this region will be loaded,
            and some points outside of the region may also be loaded.
        start_dt : datetime.datetime, default None
        end_dt : datetime.datetime, default None

        Returns
        -------
        index: dict
            A dictionary describing where to find the ADSB entries that occur in each spatial grid cell.
            If no index file exists or the index is invalid, this will be an empty dictionary.
        timestamp_ranges: dict
            A dictionary that lists for each file, a tuple storing the range of unix timestamps that are present.
        ranges: dict
            A dictionary describing for each file, which byte ranges are relevant to the spatial query.
            Only files that have been previously indexed and haven't been changed since then will be included.
            If no index file exists or the index is invalid, this will be an empty dictionary.
        """

        index_path = os.path.join(directory, self.index_file_name)
        if not os.path.exists(index_path):
            return {}, {}, {}

        with open(index_path, "r", encoding="utf-8") as f:
            # read index headers
            f.readline()
            grid_res = float(f.readline())  # grid resolution of the index
            f.readline()  # for human use only
            file_mtimes = json.loads(f.readline())  # what files are present and when they were last modified (mtime)
            f.readline()  # for human use only
            timestamp_ranges = json.loads(f.readline())  # what timestamp range each file contains
            f.readline()  # for human use only
            grid_cell_byte_offsets = json.loads(f.readline())  # what grid cells are present and their byte offset in the index file
            f.readline()  # for human use only
            index_start = f.tell()

            # make sure we are using the same grid resolution, if not, index is invalid and should be replaced
            if grid_res != self.lat_lon_grid_resolution:
                return {}, {}, {}

            # check if any files were modified since they were last indexed
            # if so, will need to remove them from the index later
            out_of_date = []
            for file in file_mtimes:
                full_path = os.path.join(directory, file)
                up_to_date = os.path.exists(full_path) and \
                    os.path.getmtime(full_path) == file_mtimes[file] and \
                    file in timestamp_ranges
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
            
            # determine which files are relevant to the temporal query
            # again, if something is out of date, we need to read the whole index
            files_to_read = []
            if len(out_of_date) > 0:
                files_to_read = list(file_mtimes.keys())
            else:
                for file in timestamp_ranges:
                    file_dts = pd.to_datetime(timestamp_ranges[file], unit="s")
                    no_overlap = (start_dt is not None and start_dt > max(file_dts)) or \
                        (end_dt is not None and end_dt < min(file_dts))
                    if not no_overlap:
                        files_to_read.append(file)

            # read the necessary grid cells
            index = {}
            for grid_cell in cells_to_read:
                if grid_cell not in grid_cell_byte_offsets:
                    continue
                # go to the position of this grid cell in the file
                cell_offset = grid_cell_byte_offsets[grid_cell]
                f.seek(index_start + cell_offset)

                # determine the byte offsets for each file in this grid cell
                f.readline()  # human-readable grid cell label
                file_offsets = json.loads(f.readline())
                file_ranges_start = f.tell()

                # for each file that we want to read, go to its location in this grid cell and read it
                grid_cell_json = {}
                for file in file_offsets:
                    if file not in files_to_read:
                        continue
                    f.seek(file_ranges_start + file_offsets[file])
                    grid_cell_json = grid_cell_json | json.loads(f.readline())

                index = index | {grid_cell: grid_cell_json}
            
            # for each file, determine which ranges should be read
            ranges = {}
            # make sure each file gets a range, even if it has no records in the region of interest
            # to communicate that it was indexed before
            for file in file_mtimes:
                ranges[file] = []
            for grid_cell in region_cells:  # we may have read more cells into the index than that match the region query, so iterate over just region cells
                if grid_cell not in index:
                    continue
                for file in index[grid_cell]:
                    ranges[file] += index[grid_cell][file]

            for file in ranges:
                # clear the ranges from files not matching the datetime query
                if file in timestamp_ranges:
                    file_dts = pd.to_datetime(timestamp_ranges[file], unit="s")
                    if (start_dt is not None and start_dt > max(file_dts)) or \
                        (end_dt is not None and end_dt < min(file_dts)):
                        ranges[file] = []

                # sort ranges by start offset, probably helps a bit with file-reading speed
                ranges[file].sort(key=lambda x: x[0])
            
            # delete out of date index records and re-save the index
            if len(out_of_date) > 0:
                for file in out_of_date:
                    for grid_cell in index:
                        if file in index[grid_cell]:
                            del index[grid_cell][file]
                    # also remove it from timestamp ranges
                    if file in timestamp_ranges:
                        del timestamp_ranges[file]
                    # also remove it from ranges, its range isn't valid anymore
                    del ranges[file]
                self._save_index(directory, index, timestamp_ranges)
        
        return index, timestamp_ranges, ranges

    def _save_index(self, directory: str, index: dict, timestamp_ranges: dict):
        """Updates/creates an ADSB index file in a directory containing ADSB TSV files.

        Parameters
        ----------
        directory: str
            A directory containing ADSB files
        index: dict
            An index dict containing information about the TSV files in the directory.
            Note that this can be partial information, not all TSV files must be represented.
            The index dict will be merged with existing index data in the index file if it exists.
        timestamp_ranges: dict
            A dict listing the timestamp range for each file. Keys are file basenames, values are
            tuples of integers representing the start/end unix timestamp.
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

            # write file timestamp range information
            f.write("file timestamp ranges\n")
            json.dump(timestamp_ranges, f)
            f.write("\n")

            # prepare chunks of the index for writing
            index_chunks = []
            grid_byte_offsets = {}
            cur_grid_offset = 0

            for grid_cell in index:
                lines = []
                file_byte_offsets = {}
                cur_file_offset = 0

                for file in index[grid_cell]:
                    line = json.dumps({file: index[grid_cell][file]}) + "\n"
                    lines.append(line)
                    file_byte_offsets[file] = cur_file_offset
                    cur_file_offset += len(line.encode("utf-8"))

                chunk = f"{grid_cell}\n{json.dumps(file_byte_offsets)}\n{''.join(lines)}"
                index_chunks.append(chunk)
                grid_byte_offsets[grid_cell] = cur_grid_offset
                cur_grid_offset += len(chunk.encode("utf-8"))

            # write header describing which grid cells are present and their byte offset in the index file
            f.write("grid cell byte offsets in this file\n")
            json.dump(grid_byte_offsets, f)
            f.write("\n")

            # write the index
            f.write("index\n")
            f.writelines(index_chunks)

    def _add_range_to_index(self, index: dict, grid_cell: str, filepath: str, start: int, length: int):
        """Utility function for inserting items into an index."""
        file = os.path.basename(filepath)
        if not grid_cell in index:
            index[grid_cell] = {}
        if not file in index[grid_cell]:
            index[grid_cell][file] = []
        index[grid_cell][file].append((start, length))

    def _read_tsv_and_update_index(self, filepath: str, index: dict):
        """Reads a TSV file and updates the index."""

        # tqdm.write(f"Reading full file {filepath}")

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
                # Detect if row is trash
                # - check for a logging blip causing two rows to get collapsed, resulting in extra values in a row
                #   these extra values get collected into a list referenced by the key None
                # - check for extra header inserted by the logger
                # - check for nonsensical timestamp
                extra_values = None in row
                extra_header = row[fieldnames[0]] == fieldnames[0]
                impossible_time = False
                if not extra_values and not extra_header:
                    timestamp = pd.to_numeric(row["TIME" if "TIME" in row else "timestamp"], errors="coerce")
                    dt = pd.Timestamp(timestamp, unit="s")
                    impossible_time = pd.isna(dt) or (dt < pd.Timestamp("2019-01-01") or (dt > pd.Timestamp.now()))
                if extra_values or extra_header or impossible_time:
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

    def _read_tsv_ranges(self, filepath: str, ranges: list):
        """Reads sections of a TSV file specified by the `ranges` parameter."""
        if len(ranges) == 0:
            return pd.DataFrame()

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

    def _get_grid_cell(self, lat: float, lon: float):
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

    def _process_raw_dataframe(self, df: pd.DataFrame):
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
        
        # Keep only those records with TSLC values of 1 or 2 seconds
        df["tslc"] = df["tslc"].astype(int)
        # invalidTslc = len(
        #     df.query("tslc >= 3 or tslc == 0")) / df.shape[0] * 100
        invalid_tslc_mask = (df["tslc"] >= 3) | (df["tslc"] == 0)
        df = df[~invalid_tslc_mask]
        if len(df) == 0:
            return None
                
        # Keep only those records with realistic altitudes
        # 10000 meters = 32808 feet; this should encompass most flights
        # NOTE: some jet aircraft may be eliminated by this process
        df["altitude"] = df["altitude"].astype(int) / 1e3
        df = df[(df["altitude"] > 0) & (df["altitude"] <= 10000)]
        if len(df) == 0:
            return None
                
        # Unpack validFlags and convert the 2-byte flag field into a list of Boolean values
        flags_names = ["valid_BARO", "valid_VERTICAL_VELOCITY", "SIMULATED_REPORT", "valid_IDENT",
                "valid_CALLSIGN", "valid_VELOCITY", "valid_HEADING", "valid_ALTITUDE", "valid_LATLON"]
        int_column = df["validFlags"].apply(int, base=16).values[:, None]
        shifts = np.arange(8, -1, -1)
        bits = (int_column >> shifts) & 1
        flags_df = pd.DataFrame(bits.astype(bool), index=df.index, columns=flags_names)
        df = pd.concat(
            [df.drop("validFlags", axis=1), flags_df], axis=1)
                
        # Keep only those records with valid latlon and altitude values based on validFlags
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
        df = df[(df["valid_LATLON"]) & (df["valid_ALTITUDE"])]
        if len(df) == 0:
            return None
                
        df.replace('-', np.nan, inplace=True)
        df.dropna(how="any", axis=0, inplace=True)
                
        # Ensure remaining field values except TIME are in proper numeric format
        df["ICAO_address"] = df["ICAO_address"].astype(str)
        df["lat"] = df["lat"].astype(int) / 1e7
        df["lon"] = df["lon"].astype(int) / 1e7
        df["heading"] = df["heading"].astype(int) / 1e2
        df["hor_velocity"] = df["hor_velocity"].astype(int) / 1e2
        df["ver_velocity"] = df["ver_velocity"].astype(int) / 1e2
        
        # Convert naive Unix timestamp to naive datetime objects and re-scale selected variable values
        df["TIME"] = pd.to_datetime(df["TIME"].astype(int), unit="s")
        df["DATE"] = df["TIME"].dt.year.astype(str) + \
            df["TIME"].dt.month.astype(str).str.zfill(2) + \
            df["TIME"].dt.day.astype(str).str.zfill(2)
        
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
            subset=['ICAO_address', 'TIME', 'lat', 'lon'], inplace=True, keep='last')
                
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
        
        if not self.data.empty: # could be empty if n_numbers=[]
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

    def __init__(self, data: gpd.GeoDataFrame, id_col: str = None, datetime_col: str = None, z_col: Optional[str] = None):
        # FYI this function needs to be able to handle the case when __init__ gets called with subsets of data,
        # which occurs as a side effect whenever you boolean index into a Tracks object or try to print it (I think).
        # This is caused by subclassing GeoDataFrame

        # check if we're using it with original intended init behavior
        # if not (if only data was passed as a side effect as described above, just call super().__init__())
        if id_col is not None and datetime_col is not None:
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

        # handle init being called as a side-effect by various pandas methods
        # in that case it will pass one argument - the data to use
        if isinstance(filename, gpd.GeoDataFrame):
            data = filename
            super().__init__(data=data, crs=data.crs)
            return
        elif isinstance(filename, pd.DataFrame):
            data = filename
            super().__init__(data=data)
            return
        
        def convert_to_bool(value):
            """
            A rigorous function to ensure proper typing of incoming fields. 
            Converts <'str'> and <'int'> to <'bool'> for better import of various file versions.
            """
            if isinstance(value, bool):
                return value
            elif isinstance(value, int):
                return bool(value)
            elif isinstance(value, str):
                if value.lower() in ['true', '1']:
                    return True
                elif value.lower() in ['false', '0']:
                    return False
                else:
                    raise ValueError(f"Invalid string value for conversion: {value}")
            else:
                raise TypeError(f"Unsupported type: {type(value)}")

        if filename:
            data = gpd.read_file(filename).astype(
                {'start_dt': 'datetime64[ns]', 'end_dt': 'datetime64[ns]'})

            # Sometimes the annotation file is read in with the valid and audible columns as booleans and other times
            #  as objects depending on what values are stored.
            data["valid"] = data["valid"].apply(convert_to_bool)
            data["audible"] = data["audible"].apply(convert_to_bool)

            if only_valid:
                data = data[data.valid == True]

        else:

            data = gpd.GeoDataFrame(columns=['_id', 'start_dt', 'end_dt', 'valid', 'audible', 'geometry', 'note'],
                                    geometry='geometry', crs='epsg:4326')

        super().__init__(data=data, crs=data.crs)


class Srcid():
    """
    A pandas DataFrame wrapper class to ensure consistent SRCID data.

    The ``nvsplDate``, ``hr``, and ``secs`` columns are combined into a single DatetimeIndex for the DataFrame and dropped.
    The ``len`` column (length of the noise event) is converted to a pandas Timedelta.

    Example Resulting DataFrame
    ---------------------------

    <class 'pandas.core.frame.DataFrame'>
    DatetimeIndex: 971 entries, 2015-05-12 05:40:03 to 2015-07-14 23:55:11
    Data columns (total 10 columns):
    len         971 non-null timedelta64[ns]
    srcID       971 non-null float64
    Hz_L        971 non-null float64
    Hz_U        971 non-null int64
    MaxSPL      971 non-null float64
    SEL         971 non-null float64
    MaxSPLt     971 non-null float64
    SELt        971 non-null float64
    userName    971 non-null object
    tagDate     971 non-null datetime64[ns]
    dtypes: datetime64[ns](1), float64(6), int64(1), object(1), timedelta64[ns](1)
    memory usage: 83.4+ KB

    Parameters
    ----------
    srcid_path: str
        Path to a single SRCID file one wishes to load.

    merge: bool, default True
        
    Notes
    -----
    [1] The SRCID format, short for "Source Identification File" was initially developed by the U.S. National Park Service
          circa 2008 as the primary output of `SPLAT.exe` (Sound Pressure Level Annotation Tool).
    """

    def __init__(self, srcid_path, merge=True):

        self.srcid_path = srcid_path
        self.merge = merge
        self.data = None

        src_raw = self._read()
        self.data = src_raw # if `merge = False` then we stop here...
        
        if(merge): # ...however, by default we'd like to work with events spanning hours (or even days) as if they were annotated in their entirety
            src = self._merge_SRCID(self.data)
            src = src.astype(src_raw.dtypes.to_dict())  # convert dtypes, they got lost somehow during merging
            self.data = src
        
    def _read(self):
        """
        Read data, format and tidy fields, handle various edge cases.
        """

        with open(self.srcid_path) as f:

            # Determine version; older versions immediately start with header, newer has version comment
            firstline = f.readline()
            if not firstline.startswith(r"%%"):
                f.seek(0)   # Rewind, as we just consumed the header row.
                            # Otherwise, we consumed the comment row, which is fine

            data = pd.read_csv(f,
                                engine= "c",
                                sep= "\t",
                                # skiprows= 1,
                                parse_dates= False)

        # Combine nvsplDate, hr, secs columns into one DatetimeIndex
        dates = pd.to_datetime(data.nvsplDate)
        hrs   = pd.to_timedelta(data.hr, unit= "h")
        secs  = pd.to_timedelta(data.secs, unit= "s")

        data.drop(["nvsplDate", "hr", "secs"], axis= 1, inplace= True)
        data.index = dates + hrs + secs

        # Turn len into timedelta
        data.len = pd.to_timedelta(data.len, unit= "s")

        if 'sID' in data.columns:
            # Some bizzare old files have a different name for srcID (really only DENAUSLC2008 so far)
            data.rename(columns= {'sID': 'srcID'}, inplace= True)

        # Days with no noise events are entered with the `MaxSPLt` and `SELt` columns skipped,
        # so they get filled with userName and tagDate, giving them a mixed type instead of float
        # Resolve this by finding zero-noise rows and setting all their values to NaN (more appropriate than 0)
        if 'MaxSPLt' in data.columns:
            if data.MaxSPLt.dtype == np.object_:
                # There are noise-free days in this dataset
                converted = pd.to_numeric(data[['MaxSPLt', 'SELt']], errors= "coerce")
                noisefree = converted.isnull().all(axis= 1)
                data.loc[noisefree, ["userName", "tagDate"]] = data.loc[noisefree, ['MaxSPLt', 'SELt']].values
                data[['MaxSPLt', 'SELt']] = converted

                nanCols = data.columns.difference(("userName", "tagDate"))
                data.loc[ noisefree, nanCols ] = np.nan

        # Parse tagDate to datetime (though old versions don't have tagDate)
        if 'tagDate' in data.columns:
            data.tagDate = pd.to_datetime(data.tagDate)

        return data

    def _join_srcID_rows(self, df):

        '''
        Join a sequence of SRCID-based spectrogram annotations into a single annotation.
        '''

        new_begin = df.head(1).index[0]
        new_length = df['len'].sum()
        new_srcID = df["srcID"].values[0]
        new_L = df['Hz_L'].min()
        new_U = df['Hz_U'].max()
        new_MaxA = "{0:.01f}".format(df["MaxSPL"].max())

        new_user = df.head(1)["userName"].values[0]
        new_tagdate = df.tail(1)["tagDate"].values[0]

        # because SEL values are already normalized you just logarithmically add them
        # this gives the total energy dose
        new_SELA = "{0:.01f}".format(10*np.log10(np.power(10, df["SEL"]/10).sum()))

        # repeat for truncated values
        # older SRCID files do not have these values, resulting in a KeyError, so we exept those cases
        try:
            new_MaxT = "{0:.01f}".format(df["MaxSPLt"].max())
            new_SELT = "{0:.01f}".format(10*np.log10(np.power(10, df["SELt"]/10).sum()))

            joined = pd.DataFrame([new_length, new_srcID, new_L, new_U, new_MaxA, new_SELA, new_MaxT, new_SELT, new_user, new_tagdate], 
                        columns=[new_begin], index=df.columns[:-1]).T

        except KeyError:
            joined = pd.DataFrame([new_length, new_srcID, new_L, new_U, new_MaxA, new_SELA, new_user, new_tagdate], 
                columns=[new_begin], index=df.columns[:-1]).T

        return joined
    
    def _merge_SRCID(self, src):

        '''
        Find SRCID-based spectrogram annotations that break across hours, and join them to create a
        final, merged SRCID for more accurate calculations.
        '''

        # these are only the events that end during the last second of the hour
        end_at_hour = src.loc[((src.index + src["len"]).dt.minute==59)&
                            ((src.index + src["len"]).dt.second==59)]

        start_at_hour = src.loc[(src.index.minute==0)&(src.index.second==0)]

        if((len(end_at_hour) > 0)&(len(start_at_hour) > 0)):

            # these are only the events that start during the first second of the hour
            starts_at_hour_break = src.index[(src.index.minute==0)&(src.index.second==0)].to_series()

            # IMPORTANT: these are the starts where there is definitely an ending one second before
            # conveniently the next two lines handle both (1) matching srcid values, 
            # and (2) lots of back-to-back-to-back hour long annotations!
            matches_end = start_at_hour.copy()
            matches_end.index = (start_at_hour.index + start_at_hour["len"] + dt.timedelta(seconds=1))

            # the .isin method is extremely helpful for working backwards to get the matched starts
            matches_start = end_at_hour[(end_at_hour.index + end_at_hour["len"] + dt.timedelta(seconds=1)).isin(start_at_hour.index)].copy()

            # these are the real break-point annotations
            cons = pd.concat([matches_start, matches_end]).sort_index().dropna()
            cons = cons.drop_duplicates(keep="first") # frankly, it doesn't matter

            # assign groups to each set of consecutive annotations 
            # it's a huge benefit to group by source type first!
            group = 1
            for srcID, source_group in cons.groupby("srcID"):

                for ts, annotation in source_group.sort_index().iterrows():

                    ends = ts + annotation["len"]

                    if((ts.minute!=0)&(ts.second!=0)&
                    (ends.minute==59)&(ends.second==59)):

                        group = group + 1
                        cons.loc[ts, "group"] = group

                    elif((ts.minute==0)&(ts.second==0)&
                    (ends.minute!=59)&(ends.second!=59)):

                        cons.loc[ts, "group"] = group
                        group = group + 1

                    else:
                        cons.loc[ts, "group"] = group

                # at the end of the current source type, 
                # "flush" the current group
                group = group + 1

            frames = []

            # now we can actually perform the joining operations
            for group_number, pieces in cons.groupby('group'):

                if(len(pieces) < 2):
                    # this removes single values - they don't actually need to be merged
                    cons.drop(pieces.index, inplace=True)

                else:
                    frames.append(self._join_srcID_rows(pieces))

            # here are all the joined data
            merged_breaks = pd.concat(frames)

            # now that everything is neat and tidy, we can get the lines
            # not representing true breaks
            no_breaks = src.loc[~src.index.isin(cons.index)]

            # final SRCID file with events across hour breaks merged
            merged_src = pd.concat([merged_breaks, no_breaks])
            merged_src = merged_src.sort_index()

        elif((len(end_at_hour) == 0)|(len(start_at_hour) == 0)):
            merged_src = src.copy()


        return merged_src
    
    def get_observation_periods(self):
        """
        Get periods of time with an acoustic record (there may be gaps). Periods are separated by
        the passing of at least a full day (midnight to midnight) without any SPLAT annotations.

        Returns
        -------
        periods: np.ndarray
            Numpy array of shape (# periods, 2) containing pairs of strings "yyyy-mm-dd" bounding each period.
            The start and end dates are inclusive, meaning there are annotations on both of these dates.
        """
        annotated_dates = np.sort(np.unique(self.data.index.date))

        # make np array of all dates
        # add one extra day to end date for np.arange exclusive indexing
        all_dates = np.arange(annotated_dates[0], annotated_dates[-1] + dt.timedelta(days=1))

        # Get contiguous regions. Returns a numpy array of shape (len(all_dates), 2),
        # containing indices of the beginning and end of each contiguous region (exclusive end index)
        period_indices = contiguous_regions(np.isin(all_dates, annotated_dates))
        # convert end indices to inclusive
        period_indices[:,1] -= 1

        # index into all dates and convert to strings
        obs_periods = all_dates[period_indices]
        return np.datetime_as_string(obs_periods, unit="D")
