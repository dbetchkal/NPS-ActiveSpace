import os

from nps_active_space.utils.helpers import estimate_line_count


class TestEstimateLineCount:
    def test_file_smaller_than_sample_is_counted_exactly(self, tmp_path):
        """A file read in full needs no scaling, so the count is exact."""
        path = tmp_path / "small.tsv"
        path.write_text("".join(f"row {i}\n" for i in range(10)))

        assert estimate_line_count(str(path)) == 10

    def test_file_larger_than_sample_is_estimated(self, tmp_path):
        path = tmp_path / "large.tsv"
        path.write_text("".join(f"row {i} {'x' * 40}\n" for i in range(4000)))
        sample_size = os.path.getsize(path) // 4

        estimate = estimate_line_count(str(path), sample_size=sample_size)

        assert 3800 <= estimate <= 4200

    def test_file_without_newlines_returns_zero(self, tmp_path):
        path = tmp_path / "no_newline.tsv"
        path.write_text("single row with no trailing newline")

        assert estimate_line_count(str(path)) == 0
