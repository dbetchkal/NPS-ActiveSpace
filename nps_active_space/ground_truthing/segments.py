import datetime as dt

import geopandas as gpd
from shapely.geometry import LineString, Point

# Mutable [start, end] pair on the time_audible axis (lists are updated in place).
AudibleRange = list[dt.datetime]


def collapse_audible_ranges(ranges: list[AudibleRange]) -> list[AudibleRange]:
    """Collapse overlapping audible ranges into single audible ranges.
    
    Parameters
    ----------
    ranges: list of [datetime, datetime]
    
    Returns
    -------
    collapsed_ranges: list of [datetime, datetime]
    """
    
    times = []
    for start, end in ranges:
        times.append({"t": start, "type": "start"})
        times.append({"t": end, "type": "end"})
    times.sort(key=lambda x: x["t"])
    
    intervals = []
    n = 0
    start_time = None
    for t in times:
        if t["type"] == "start":
            # if it was previously quiet, start the next interval
            if n == 0:
                start_time = t["t"]
            n += 1
        elif t["type"] == "end":
            n -= 1
            if n == 0:
                # last sound ended, save the interval
                intervals.append([start_time, t["t"]])

    return intervals


def build_annotation_segments(
    track_id: str,
    points: gpd.GeoDataFrame,
    audible_ranges: list[AudibleRange] | None = None,
    valid: bool = True,
    note: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Build annotation segment GeoDataFrame from track points and audible ranges.

    Parameters
    ----------
    track_id : str
        The track unique identifier.
    points: gpd.GeoDataFrame:
        Track and spline points to annotate.
    valid : bool, default True
        If the track was valid.
    audible_ranges: list of [datetime, datetime], default None
        Periods of time when the track was audible. If None or empty, everything was inaudible.
    note: str, default None
        Any note to be added to all points passed for annotation.
    """
    if audible_ranges is None:
        audible_ranges = []

    # Convert points to WGS84 to avoid geopandas bug mentioned in Track model :(
    if 'z' not in points.columns:
        points['z'] = points.geometry.z
        points = points.to_crs('epsg:4326')
        points['geometry'] = points.apply(lambda row: Point(row.geometry.x, row.geometry.y, row.z), axis=1)
        points.drop('z', axis=1, inplace=True)
    

    segments = []

    # Saving invalid tracks or fully inaudible tracks can be done with one line
    if valid is False or len(audible_ranges) == 0:
        segments.append({
            '_id': track_id,
            'start_dt': points.point_dt.iat[0],
            'end_dt': points.point_dt.iat[-1],
            'valid': valid,
            'audible': False,
            'note': note,
            'geometry': points.geometry.iat[0] if points.shape[0] == 1
                            else LineString(points.geometry.tolist())
        })

    else:
        # note that in this conditional block, audible_ranges is not empty

        # simplify overlapping audible ranges
        audible_ranges = collapse_audible_ranges(audible_ranges)

        # remove timezone info to avoid errors comparing tz-naive against tz-aware datetimes
        for r in audible_ranges:
            for i in range(len(r)):
                r[i] = r[i].replace(tzinfo=None)

        segments = []
        
        # add segment for the inaudible range at the beginning, if it exists
        first_inaudible_segment = points[points.time_audible < audible_ranges[0][0]]
        if first_inaudible_segment.shape[0] >= 2:
            segments.append(
                {'_id': track_id,
                'start_dt': first_inaudible_segment.point_dt.iat[0],
                'end_dt': first_inaudible_segment.point_dt.iat[-1],
                'valid': True,
                'audible': False,
                'note': note,
                'geometry': LineString(first_inaudible_segment.geometry.tolist())}
            )

        # add segments for each audible range and the inaudible range following it
        for i, r in enumerate(audible_ranges):                
            audible_segment = points[(points.time_audible >= r[0]) & (points.time_audible < r[1])]
            if audible_segment.shape[0] >= 2:
                segments.append(
                    {'_id': track_id,
                    'start_dt': audible_segment.point_dt.iat[0],
                    'end_dt': audible_segment.point_dt.iat[-1],
                    'valid': True,
                    'audible': True,
                    'note': note,
                    'geometry': LineString(audible_segment.geometry.tolist())}
                )

            if i+1 < len(audible_ranges):
                next_start = audible_ranges[i+1][0]
            else:
                next_start = points.time_audible.iat[-1]
            
            inaudible_segment = points[(points.time_audible >= r[1]) & (points.time_audible < next_start)]
            if inaudible_segment.shape[0] >= 2:
                segments.append(
                    {'_id': track_id,
                    'start_dt': inaudible_segment.point_dt.iat[0],
                    'end_dt': inaudible_segment.point_dt.iat[-1],
                    'valid': True,
                    'audible': False,
                    'note': note,
                    'geometry': LineString(inaudible_segment.geometry.tolist())}
                )

    gdf = gpd.GeoDataFrame(segments, geometry='geometry', crs=points.crs)
    return gdf
