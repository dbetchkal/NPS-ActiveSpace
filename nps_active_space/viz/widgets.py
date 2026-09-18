from __future__ import annotations

from vtkmodules.vtkRenderingCore import vtkTextActor

from nps_active_space.viz.session import VizSession


class PlotWidgets:
    def __init__(self, session: VizSession) -> None:
        self._session = session

    @property
    def _style(self):
        return self._session.style

    def layer_checkbox_xy(self, row: int) -> tuple[int, int]:
        return (
            self._style.layer_checkbox_x,
            self._style.layer_checkbox_y0 + self._style.layer_widget_dy * row,
        )

    def add_labeled_checkbox(
        self,
        callback,
        *,
        value: bool,
        position: tuple[int, int],
        size: int,
        color_on: str,
        label: str,
    ):
        """Checkbox plus a 2D label in the same VTK display-pixel space."""
        plotter = self._session.plotter
        checkbox = plotter.add_checkbox_button_widget(
            callback=callback,
            value=value,
            position=position,
            size=size,
            color_on=color_on,
        )
        x, y = position
        text = vtkTextActor()
        text.SetInput(label)
        text.SetTextScaleModeToNone()
        text.GetPositionCoordinate().SetCoordinateSystemToDisplay()
        text.SetPosition(x + size + 8, y + max(2, (size - 16) // 2))
        prop = text.GetTextProperty()
        prop.SetFontFamilyToArial()
        prop.SetFontSize(16)
        prop.SetColor(1.0, 1.0, 1.0)
        prop.BoldOn()
        prop.ShadowOn()
        prop.SetJustificationToLeft()
        prop.SetVerticalJustificationToBottom()
        plotter.renderer.AddViewProp(text)
        return checkbox

    def add_color_legend(self, *, compare_models: bool) -> None:
        legend_entries = self._session.legend_models
        if not legend_entries and compare_models:
            legend_entries = [
                ("NMSim", self._style.nmsim_activespace_color),
                ("AAM", self._style.aam_activespace_color),
            ]
        if not legend_entries:
            return
        unique: list[tuple[str, str]] = []
        seen: set[str] = set()
        for name, color in legend_entries:
            if name not in seen:
                unique.append((name, color))
                seen.add(name)
        self._session.plotter.add_legend(
            unique,
            loc="upper right",
            bcolor=None,
            face="r",
        )

    def setup_z_scale(self) -> None:
        plotter = self._session.plotter
        plotter.set_scale(1, 1, 2)

        def toggle_z_scale(flag):
            plotter.set_scale(1, 1, 2 if flag else 1)

        self.add_labeled_checkbox(
            toggle_z_scale,
            value=True,
            position=(10, 140),
            size=25,
            color_on=self._style.z_scale_toggle_color,
            label="2× z-scale",
        )
