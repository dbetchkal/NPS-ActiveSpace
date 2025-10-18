#geospatial libraries
import gdal
from gdalconst import GA_ReadOnly
import geopandas as gpd
import pyproj
import rasterio
import rasterio.plot
from rasterio.windows import Window
from shapely.geometry import Point, LineString

# other libraries
import datetime as dt
import glob
from itertools import islice
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.widgets import RangeSlider, Button, Slider, RadioButtons
import numpy as np
import os
import pandas as pd
import pickle
import re
from scipy import interpolate
import sqlalchemy
import subprocess
import sys
import warnings
import geopy as geopy
from geopy.distance import geodesic
from time import mktime
import math
import ipykernel
