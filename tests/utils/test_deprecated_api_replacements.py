"""
Tests verifying that the deprecated API replacements in the codebase
produce identical or correct behavior compared to the old calls.

Covers three replacements:
1. pd.read_csv delim_whitespace=True -> sep=r'\\s+'
2. DataFrame.append() -> pd.concat()
3. gpd.sjoin(op='within') -> gpd.sjoin(predicate='within')
"""

import textwrap

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point, Polygon


class TestReadCsvWhitespaceSeparator:
    """Verify that sep=r'\\s+' behaves the same as the old delim_whitespace=True."""

    def test_basic_whitespace_separated(self, tmp_path):
        """Spaces and tabs between columns should parse identically."""
        data_file = tmp_path / "data.txt"
        data_file.write_text(textwrap.dedent("""\
            col_a col_b col_c
            1 2 3
            4 5 6
        """))

        result = pd.read_csv(data_file, sep=r"\s+")

        assert list(result.columns) == ["col_a", "col_b", "col_c"]
        assert result.shape == (2, 3)
        assert result["col_a"].tolist() == [1, 4]
        assert result["col_b"].tolist() == [2, 5]
        assert result["col_c"].tolist() == [3, 6]

    def test_mixed_whitespace(self, tmp_path):
        """Mixed spaces and tabs should be handled by the regex separator."""
        data_file = tmp_path / "mixed.txt"
        data_file.write_text("x\ty\n10\t  20\n30 40\n")

        result = pd.read_csv(data_file, sep=r"\s+")

        assert result.shape == (2, 2)
        assert result["x"].tolist() == [10, 30]
        assert result["y"].tolist() == [20, 40]

    def test_header_file_pattern(self, tmp_path):
        """Mimics the .hdr file reading pattern used in the codebase."""
        hdr_file = tmp_path / "test.hdr"
        hdr_file.write_text(textwrap.dedent("""\
            NROWS 100
            NCOLS 200
            XDIM 30.0
            ULYMAP 5000.0
        """))

        df = pd.read_csv(hdr_file, header=None, sep=r"\s+", index_col=0).T

        assert float(df.NROWS) == 100
        assert float(df.NCOLS) == 200
        assert float(df.XDIM) == 30.0
        assert float(df.ULYMAP) == 5000.0


class TestConcatReplacesAppend:
    """Verify that pd.concat produces the same result as the old .append() call."""

    def test_concat_two_dataframes(self):
        """Basic concatenation with ignore_index=True."""
        df1 = pd.DataFrame({"geometry": ["POINT(0 0)"], "value": [1]})
        df2 = pd.DataFrame({"geometry": ["POINT(1 1)"], "value": [2]})

        result = pd.concat([df1, df2], ignore_index=True)

        assert len(result) == 2
        assert result["value"].tolist() == [1, 2]
        assert result.index.tolist() == [0, 1]

    def test_concat_preserves_columns(self):
        """Columns from both frames should be present in the output."""
        df1 = pd.DataFrame({"a": [1], "b": [2]})
        df2 = pd.DataFrame({"a": [3], "b": [4]})

        result = pd.concat([df1, df2], ignore_index=True)

        assert list(result.columns) == ["a", "b"]
        assert result["a"].tolist() == [1, 3]
        assert result["b"].tolist() == [2, 4]

    def test_concat_accumulation_loop(self):
        """Simulates the accumulation loop in generate_active_space_mesh.py."""
        accumulated = None
        chunks = [
            pd.DataFrame({"val": [10, 20]}),
            pd.DataFrame({"val": [30]}),
            pd.DataFrame({"val": [40, 50, 60]}),
        ]

        for chunk in chunks:
            if accumulated is None:
                accumulated = chunk
            else:
                accumulated = pd.concat([accumulated, chunk], ignore_index=True)

        assert len(accumulated) == 6
        assert accumulated["val"].tolist() == [10, 20, 30, 40, 50, 60]
        assert accumulated.index.tolist() == list(range(6))


class TestSjoinPredicateParameter:
    """Verify that gpd.sjoin with predicate='within' returns the correct spatial join."""

    @pytest.fixture
    def polygon_gdf(self):
        """A single polygon covering (0,0) to (10,10)."""
        poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        return gpd.GeoDataFrame(
            {"zone": ["inside_zone"]},
            geometry=[poly],
            crs="EPSG:4326",
        )

    @pytest.fixture
    def points_gdf(self):
        """Three points: two inside the polygon, one outside."""
        return gpd.GeoDataFrame(
            {"name": ["inside_1", "inside_2", "outside"]},
            geometry=[Point(5, 5), Point(2, 3), Point(20, 20)],
            crs="EPSG:4326",
        )

    def test_sjoin_within_filters_correctly(self, points_gdf, polygon_gdf):
        """Only points within the polygon should appear in the result."""
        result = gpd.sjoin(points_gdf, polygon_gdf, predicate="within")

        assert len(result) == 2
        assert set(result["name"].tolist()) == {"inside_1", "inside_2"}

    def test_sjoin_within_excludes_outside(self, points_gdf, polygon_gdf):
        """The point at (20, 20) should not be in the result."""
        result = gpd.sjoin(points_gdf, polygon_gdf, predicate="within")

        assert "outside" not in result["name"].tolist()

    def test_sjoin_within_column_selection(self, points_gdf, polygon_gdf):
        """Selecting only geometry after sjoin (matches the codebase pattern)."""
        result = gpd.sjoin(points_gdf, polygon_gdf, predicate="within")[["geometry"]]

        assert list(result.columns) == ["geometry"]
        assert len(result) == 2
