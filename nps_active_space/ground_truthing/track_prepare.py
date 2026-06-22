import datetime as dt
import warnings
from typing import TYPE_CHECKING

import geopandas as gpd
import numpy as np
import pandas as pd
from matplotlib.dates import date2num, num2date
from shapely.geometry import Point

from nps_active_space.ground_truthing.segments import AudibleRange
from nps_active_space.utils.computation import audible_time_delay, interpolate_spline, expected_Lp

if TYPE_CHECKING:
    from nps_active_space.utils.models import Nvspl


def load_spectrogram(
    nvspl: "Nvspl",
    points: gpd.GeoDataFrame,
    time_pad: dt.timedelta | None = None,
) -> pd.DataFrame:
    """Load spectrogram window for track points."""
    if time_pad is None:
        time_pad = dt.timedelta(seconds=5*60)
    t_start = str(points.point_dt.iat[0] - time_pad)
    t_end = str(points.point_dt.iat[-1] + time_pad)
    spectro = nvspl.loc[t_start:t_end, '12.5':'20000']
    return spectro


def prepare_spline(points: gpd.GeoDataFrame, mic_point: Point) -> gpd.GeoDataFrame:
    """Sort points, interpolate spline, and compute audibility fields."""
    points.sort_values(by='point_dt', ascending=True, inplace=True)
    spline = interpolate_spline(points)
    spline = audible_time_delay(spline, 'point_dt', mic_point)
    spline = expected_Lp(spline, mic_point)
    return spline


def closest_approach(spline: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, dt.datetime]:
    """Determine the closest spline point to the mic."""
    closest_point = spline[spline.distance_to_target == spline.distance_to_target.min()]
    closest_time = spline.loc[spline.distance_to_target.idxmin()]['time_audible']
    return closest_point, closest_time


def limit_line_bounds(
    closest_time: dt.datetime,
    x_lims: np.ndarray,
) -> tuple[float, float]:
    """Calculate datetime starting points for limit lines."""
    lower_limit_start = max(date2num(closest_time - dt.timedelta(seconds=60)), x_lims[0])
    upper_limit_start = min(date2num(closest_time + dt.timedelta(seconds=60)), x_lims[-1])
    return lower_limit_start, upper_limit_start


def audible_ranges_from_annotations(
    annots: gpd.GeoDataFrame,
    spline: gpd.GeoDataFrame,
    lower_limit_start: float,
    upper_limit_start: float,
) -> list[AudibleRange] | None:
    """
    Load audible ranges from previous annotations or defaults.

    Returns None if annotation times cannot be matched within tolerance (caller should return early).
    """
    if annots.empty:
        # default range slider
        return [[num2date(lower_limit_start), num2date(upper_limit_start)]]

    audible_ranges: list[AudibleRange] = []
    for _, a in annots[annots["valid"] & annots["audible"]].iterrows():
        # note that invalid or fully inaudible annotations will result in no ranges, as desired

        # GEOJSON doesn't save full time precision (up to millisecond it seems),
        # so when trying to match an annotated start_dt or end_dt to a point's point_dt,
        # we have to only check if they're close enough, not exactly equal.
        # Also, if a clock drift file was changed slightly, the annotation start and end time
        # may not match up perfectly (but will be close enough)
        start_diff = (spline["point_dt"] - a["start_dt"]).abs()
        end_diff = (spline["point_dt"] - a["end_dt"]).abs()
        # warn if not close enough
        tol = pd.Timedelta(seconds=1)
        if start_diff.min() > tol:
            warnings.warn(f"Could not find a track point within {tol.total_seconds()} sec of annotation start")
            return None
        if end_diff.min() > tol:
            warnings.warn(f"Could not find a track point within {tol.total_seconds()} sec of annotation end")
            return None
        time_audible_start = spline.loc[start_diff.idxmin(), "time_audible"]
        time_audible_end = spline.loc[end_diff.idxmin(), "time_audible"]
        audible_ranges.append([time_audible_start, time_audible_end])

    return audible_ranges
