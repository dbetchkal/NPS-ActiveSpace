import datetime as dt

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point

from nps_active_space.utils.geometry import linestring_from_coords

# Mutable [start, end] pair on the time_audible axis (lists are updated in place).
AudibleRange = list[dt.datetime]

# Splines are 1 Hz; storing every vertex makes GeoJSON huge. Reload uses start_dt/end_dt only.
ANNOTATION_MAX_VERTICES = 64
COORD_DECIMALS = 6
Z_DECIMALS = 1


def _round_coord(coord: tuple) -> tuple:
    """Round WGS84 xy and elevation for compact GeoJSON."""
    if len(coord) < 3:
        return (round(coord[0], COORD_DECIMALS), round(coord[1], COORD_DECIMALS))
    return (
        round(coord[0], COORD_DECIMALS),
        round(coord[1], COORD_DECIMALS),
        round(coord[2], Z_DECIMALS),
    )


def downsample_linestring(
    line: LineString,
    max_vertices: int = ANNOTATION_MAX_VERTICES,
) -> LineString:
    """Reduce vertex count for annotation storage and viz (keeps endpoints)."""
    if line.is_empty:
        return line
    coords = list(line.coords)
    if len(coords) <= max_vertices:
        return linestring_from_coords([_round_coord(c) for c in coords])
    indices = np.unique(np.linspace(0, len(coords) - 1, max_vertices, dtype=int))
    sampled = [_round_coord(coords[i]) for i in indices]
    return linestring_from_coords(sampled)


def storage_geometry_from_points(geometries: list) -> Point | LineString:
    """Build a compact Point or downsampled LineString for GeoJSON export."""
    if len(geometries) == 1:
        geom = geometries[0]
        if isinstance(geom, Point):
            return Point(_round_coord(geom.coords[0]))
        if isinstance(geom, LineString):
            return downsample_linestring(linestring_from_coords(list(geom.coords)))
        raise TypeError(f"Unsupported geometry type: {type(geom)!r}")
    return downsample_linestring(
        linestring_from_coords([g.coords[0] for g in geometries])
    )


def compact_geometry(
    geom: Point | LineString,
    max_vertices: int = ANNOTATION_MAX_VERTICES,
) -> Point | LineString:
    """Downsample and round a stored annotation geometry."""
    if isinstance(geom, Point):
        return Point(_round_coord(geom.coords[0]))
    if isinstance(geom, LineString):
        return downsample_linestring(geom, max_vertices=max_vertices)
    return geom


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
    audible_ranges: list[AudibleRange] = [],
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
    audible_ranges: list of [datetime, datetime], default []
        Periods of time when the track was audible. If an empty list, everything was inaudible.        
    note: str, default None
        Any note to be added to all points passed for annotation.
    """
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
            'geometry': storage_geometry_from_points(
                [points.geometry.iat[0]] if points.shape[0] == 1
                else list(points.geometry)
            )
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
                'geometry': storage_geometry_from_points(
                    list(first_inaudible_segment.geometry)
                )}
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
                    'geometry': storage_geometry_from_points(
                        list(audible_segment.geometry)
                    )}
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
                    'geometry': storage_geometry_from_points(
                        list(inaudible_segment.geometry)
                    )}
                )

    gdf = gpd.GeoDataFrame(segments, geometry='geometry', crs=points.crs)
    return gdf
