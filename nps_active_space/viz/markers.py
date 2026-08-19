"""Scene orientation markers and viewport chrome."""

from __future__ import annotations

WINDOW_TITLE = "NPS ActiveSpace Visualization"


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
