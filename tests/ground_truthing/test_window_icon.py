from pathlib import Path

from nps_active_space.ground_truthing.window_icon import nps_logo_paths


class TestNpsLogoPaths:
    def test_packaged_logo_files_exist(self) -> None:
        png_path, ico_path = nps_logo_paths()
        assert png_path.is_file(), f"missing packaged logo: {png_path}"
        assert ico_path.is_file(), f"missing packaged logo: {ico_path}"
        assert png_path.name == "flat-four-color.png"
        assert ico_path.suffix == ".ico"
