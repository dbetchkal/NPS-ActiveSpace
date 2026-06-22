import numpy as np
import numpy.typing as npt
import pyproj
import rasterio
from rasterio.io import DatasetReader


def process_dem(
    dem: DatasetReader,
    crs: str,
) -> tuple[npt.NDArray[np.floating], npt.NDArray[np.floating], np.ma.MaskedArray | npt.NDArray[np.floating]]:
    data = dem.read(1)

    if dem.nodata is not None:
        data = np.ma.masked_equal(data, dem.nodata)
    
    # convert pixel coords to spatial coords
    x = np.arange(0, data.shape[1])
    y = np.arange(0, data.shape[0])
    x_coords, y_coords = np.meshgrid(x, y)  # pixel coords right now
    x_coords, y_coords = rasterio.transform.xy(dem.transform, y_coords, x_coords, offset="center") # now spatial coords
    x_coords = x_coords.reshape(data.shape)
    y_coords = y_coords.reshape(data.shape)

    # convert x,y to correct CRS
    transformer = pyproj.Transformer.from_crs(dem.crs, crs, always_xy=True)
    x_coords, y_coords = transformer.transform(x_coords, y_coords)

    return x_coords, y_coords, data
