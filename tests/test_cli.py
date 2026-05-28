from nps_active_space.cli import _apply_console_argv


class TestApplyConsoleArgv:
    def test_default_leaves_argv_and_hides_warnings(self):
        argv = ["nps-ground-truthing", "-e", "DENA"]
        assert _apply_console_argv(argv) is False
        assert argv == ["nps-ground-truthing", "-e", "DENA"]

    def test_show_warnings_flag(self):
        argv = ["nps-viz", "--show-warnings", "-e", "prod"]
        assert _apply_console_argv(argv) is True
        assert argv == ["nps-viz", "-e", "prod"]

    def test_ignore_warnings_overrides_show(self):
        argv = ["nps-viz", "--show-warnings", "--ignore-warnings", "-e", "prod"]
        assert _apply_console_argv(argv) is False
        assert argv == ["nps-viz", "-e", "prod"]
