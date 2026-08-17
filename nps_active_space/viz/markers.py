"""Scene orientation markers and viewport chrome."""

from __future__ import annotations

import sys
from pathlib import Path

import pyvista as pv

WINDOW_TITLE = "NPS ActiveSpace Visualization"
WINDOWS_APP_USER_MODEL_ID = "gov.nps.activespace.viz"
_ASSETS_DIR = Path(__file__).resolve().parent / "assets"
DEFAULT_WINDOW_ICON = _ASSETS_DIR / "waypoint_icon.png"
_taskbar_identity_registered = False


def default_window_icon_path() -> Path:
    """Bundled waypoint/map-pin PNG for the native window icon."""
    return DEFAULT_WINDOW_ICON


def window_icon_supported() -> bool:
    """VTK implements SetIcon on Windows and Linux OpenGL render windows only."""
    return sys.platform != "darwin"


def register_windows_taskbar_identity(
    app_id: str = WINDOWS_APP_USER_MODEL_ID,
) -> None:
    """Register a Windows AppUserModelID so the taskbar uses our window icon.

    Without this, Windows groups the process under python.exe and keeps the
    Python interpreter icon in the taskbar even when VTK SetIcon updates the
    title-bar icon. Must run before any GUI window is created.
    """
    global _taskbar_identity_registered
    if _taskbar_identity_registered or sys.platform != "win32":
        return
    import ctypes

    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
    _taskbar_identity_registered = True


def apply_window_icon(plotter: pv.Plotter, icon_path: Path | None = None) -> None:
    """Apply the window icon when a PNG is available (no-op if missing or unsupported)."""
    if not window_icon_supported():
        return
    path = DEFAULT_WINDOW_ICON if icon_path is None else icon_path
    if path.is_file():
        set_window_icon(plotter, path)


def utm_orientation_axes_kwargs() -> dict:
    """Orientation-marker kwargs for UTM scenes (+X east, +Y north).

    PyVista/VTK allow only one vtkOrientationMarkerWidget per renderer, so axes
    and north-arrow widgets cannot coexist as separate corner widgets.
    """
    return {
        "interactive": False,
        "line_width": 2,
        "xlabel": "E",
        "ylabel": "N",
        "zlabel": "Z",
        "viewport": (0, 0, 0.2, 0.2),
    }


def set_window_icon(plotter: pv.Plotter, icon_path: str | Path) -> None:
    """Set the native window icon from a PNG file.

    Works on Windows and many Linux window managers. macOS often ignores VTK
    window icons; use a bundled .icns via a Qt wrapper if dock icon matters.
    """
    if not window_icon_supported():
        return

    from vtkmodules.vtkIOImage import vtkPNGReader

    path = Path(icon_path)
    if not path.is_file():
        raise FileNotFoundError(path)
    reader = vtkPNGReader()
    reader.SetFileName(str(path))
    reader.Update()
    # VTK requires an initial render before SetIcon on some platforms (see PyVista #73531355).
    if not plotter.off_screen:
        plotter.render_window.Render()
    plotter.render_window.SetIcon(reader.GetOutput())
