import glob
import re
import os
from dataclasses import dataclass, field
from typing import List, Optional, Union

import geopandas as gpd
import pandas as pd
from pyproj import Transformer
from tqdm import tqdm


__all__ = [
    'Microphone',
    'Nvspl',
    'Tracks'
]


@dataclass
class Microphone:
    """An object representing a microphone deployment location and time."""

    unit: str
    site: str
    year: int   # TODO: Should this be a date instead of year?
    lat: float  # WGS84, epsg=4326
    lon: float  # WGS84, epsg=4326
    z: float    # In meters
    name: str = None
    crs: str = None
    x: float = field(init=False)
    y: float = field(init=False)

    def __repr__(self):
        return f"Microphone(name={self.name})"

    def __post_init__(self):
        """Set x,y coordinates and instance name."""
        if self.crs:
            self.to_crs(self.crs)
        if not self.name:
            self.name = f"{self.unit}{self.site}{self.year}"

    def to_crs(self, crs: str):
        """Project instance x,y values to a new coordinate system.

        Parameters
        ----------
        crs : str
            The coordinate system to project the instance to. Format: epsg:XXXX. E.g. epsg:26906
        """
        projection = Transformer.from_crs('epsg:4326', crs, always_xy=True)
        self.x, self.y = projection.transform(self.lon, self.lat)
        self.crs = crs


class Nvspl(pd.DataFrame):
    """A pandas DataFrame wrapper class to ensure consistent NVSPL data."""

    standard_fields = {
        'SiteID', 'STime', 'dbA', 'dbC', 'dbF',
        'Voltage', 'WindSpeed', 'WindDir', 'TempIns',
        'TempOut', 'Humidity', 'INVID', 'INSID',
        'GChar1', 'GChar2', 'GChar3', 'AdjustmentsApplied',
        'CalibrationAdjustment', 'GPSTimeAdjustment',
        'GainAdjustment', 'Status'
    }

    octave_regex = re.compile(r"^H[0-9]+$|^H[0-9]+p[0-9]$")

    def __init__(self, filepaths_or_data: Union[List[str], str, pd.DataFrame]):
        """

        Parameters
        ----------
        filepaths_or_data : List, str, or pd.DataFrame
            A directory containing NVSPL files, a list of NVSPL files, or an existing pd.DataFrame of
            NVSPL data.
        """
        data = self._read(filepaths_or_data)
        data.set_index('STime', inplace=True)
        super().__init__(data=data)

    def _read(self, filepaths_or_data):
        # TODO: for speed, usecols and define datatype, drop empty columns

        if isinstance(filepaths_or_data, pd.DataFrame):
            self._validate(filepaths_or_data.columns)
            return filepaths_or_data

        elif isinstance(filepaths_or_data, str):
            assert os.path.isdir(filepaths_or_data), f"{filepaths_or_data} does not exist."
            filepaths_or_data = glob.glob(f"{filepaths_or_data}/*.txt")

        else:
            for file in filepaths_or_data:
                assert os.path.isfile(file), f"{file} does not exist."
                assert file.endswith('.txt'), f"Only .txt NVSPL files accepted."

        data = pd.DataFrame()
        for file in tqdm(filepaths_or_data, desc='Loading NVSPL files', unit='files', colour='green'):
            df = pd.read_csv(file)
            self._validate(df.columns)
            data = data.append(df)

        return data

    def _validate(self, columns):
        """
        Parameters
        ----------
        columns
        """
        # Verify that all NVSPL standard columns exist.
        missing_standard_cols = self.standard_fields - set(columns)
        assert missing_standard_cols == set(), f"Missing the following standard NVSPL columns: {missing_standard_cols}"

        # Verify all non-standard columns are octave columns.
        only_standard_cols = all(re.match(self.octave_regex, col) for col in (set(columns) - self.standard_fields))
        assert only_standard_cols is True, "NVSPL data contains unexpected NVSPL columns."


class Tracks(gpd.GeoDataFrame):
    """A geopandas GeoDataFrame wrapper class to standardize track points.

    Parameters
    ----------
    data : gpd.GeoDataFrame
    id_col : str, default None
        name of the columns
    """

    def __init__(self, data: gpd.GeoDataFrame, id_col: Optional[str] = None):
        # TODO: time column

        data.set_index(id_col, inplace=True)

        super().__init__(data=self._read(data))
