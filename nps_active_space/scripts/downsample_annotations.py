"""Downsample dense annotation GeoJSON (1 Hz splines) for smaller files."""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString

from nps_active_space.ground_truthing.segments import (
    ANNOTATION_MAX_VERTICES,
    compact_geometry,
)


def _vertex_count(geom) -> int:
    if isinstance(geom, LineString):
        return len(geom.coords)
    return 1


def downsample_annotations_file(
    input_path: Path,
    output_path: Path,
    *,
    max_vertices: int = ANNOTATION_MAX_VERTICES,
) -> tuple[int, int]:
    """Read annotation GeoJSON, compact geometries, write output. Returns (before, after) vertex totals."""
    gdf = gpd.read_file(input_path)
    before = sum(_vertex_count(g) for g in gdf.geometry)
    gdf = gdf.copy()
    gdf["geometry"] = gdf.geometry.map(
        lambda g: compact_geometry(g, max_vertices=max_vertices)
    )
    after = sum(_vertex_count(g) for g in gdf.geometry)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(output_path, driver="GeoJSON")
    return before, after


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Downsample dense ground-truthing annotation GeoJSON."
    )
    parser.add_argument("input", type=Path, help="Input .geojson")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="Output .geojson (default: <input_stem>_compact.geojson)",
    )
    parser.add_argument(
        "--max-vertices",
        type=int,
        default=ANNOTATION_MAX_VERTICES,
        help=f"Max vertices per segment (default: {ANNOTATION_MAX_VERTICES})",
    )
    args = parser.parse_args()

    input_path = args.input.resolve()
    if not input_path.is_file():
        parser.error(f"input not found: {input_path}")

    output_path = (
        args.output.resolve()
        if args.output is not None
        else input_path.with_name(f"{input_path.stem}_compact.geojson")
    )

    before_bytes = input_path.stat().st_size
    before_verts, after_verts = downsample_annotations_file(
        input_path,
        output_path,
        max_vertices=args.max_vertices,
    )
    after_bytes = output_path.stat().st_size

    print(f"Wrote {output_path}")
    print(f"  vertices: {before_verts:,} -> {after_verts:,}")
    print(f"  file size: {before_bytes / 1e6:.1f} MB -> {after_bytes / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
