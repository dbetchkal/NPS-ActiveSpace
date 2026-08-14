import datetime as dt
from functools import partial
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
import shapely
import tkinter as tk
from matplotlib.dates import date2num, DateFormatter, num2date
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Button

from nps_active_space.ground_truthing.widgets import FastRangeSlider
from nps_active_space.utils.enums import TrackSource

if TYPE_CHECKING:
    from nps_active_space.ground_truthing.annotation_frames import _GroundTruthingFrame


def _track_label_type_lines(
    track_source: TrackSource,
    vehicle_type: str | None,
    vessel_name: str | None,
) -> tuple[str, str]:
    """Return type and name label lines for the track info panel."""
    match track_source:
        case TrackSource.AIS:
            type_line = (
                f"Vessel type: {vehicle_type}\n"
                if vehicle_type is not None else ""
            )
            name_line = (
                f"Vessel: {vessel_name}\n"
                if vessel_name is not None else ""
            )
        case _:
            type_line = (
                f"Aircraft Type: {vehicle_type}\n"
                if vehicle_type is not None else ""
            )
            name_line = ""
    return type_line, name_line


def _cursor_status_text(cursor_time: dt.datetime, altitude_m: float | None) -> str:
    """Format cursor time and optional track altitude for the side panel."""
    text = f"Cursor Time: {cursor_time.strftime('%H:%M:%S')}"
    if altitude_m is not None and np.isfinite(altitude_m):
        text += f"\nAltitude: {altitude_m:.0f} m MSL"
    return text


def _point_altitude_m(point_geometry: Any) -> float | None:
    """Return MSL altitude from a shapely point, or None when unavailable."""
    if not shapely.has_z(point_geometry):
        return None
    z = shapely.get_z(point_geometry)
    if not np.isfinite(z):
        return None
    return float(z)


def _track_altitudes_m(points: Any) -> np.ndarray:
    """Return finite MSL altitudes for raw track points."""
    if "z" in points.columns:
        altitudes = points["z"].to_numpy(dtype=float)
    else:
        altitudes = np.array(
            [_point_altitude_m(geom) for geom in points.geometry],
            dtype=float,
        )
    return altitudes[np.isfinite(altitudes)]


def _track_altitude_summary(points: Any, closest_point_geometry: Any) -> str:
    """Format static track altitude stats for the side panel."""
    altitudes = _track_altitudes_m(points)
    if altitudes.size == 0:
        return ""

    lines = [
        f"Track altitude: {altitudes.min():.0f}–{altitudes.max():.0f} m MSL "
        f"(mean {altitudes.mean():.0f})"
    ]
    closest_alt = _point_altitude_m(closest_point_geometry)
    if closest_alt is not None:
        lines.append(f"Closest approach: {closest_alt:.0f} m MSL")
    return "\n".join(lines) + "\n"


def build_plot(frame: "_GroundTruthingFrame") -> None:
    """
    Build the matplotlib GridSpec plot for a track, using class state set by _load_index().
    """
    # ************************************ Build Plot ************************************#

    # destroy old figure widget, to avoid incorrect background saving for blitting,
    # and probably also avoids app slowing down over time
    figsize = [9, 5]
    if frame.fig_widget is not None:
        figsize[0] = frame.fig_widget.winfo_width() / 100
        figsize[1] = frame.fig_widget.winfo_height() / 100
        frame.fig_widget.destroy()

    # clear any existing figures - matplotlib stores them internally if we don't, taking up tons of memory
    plt.close("all")

    # new figure
    fig = plt.figure(figsize=figsize, layout="constrained", dpi=100)
    fig.canvas.manager.set_window_title(f"Microphone: {frame.master.mic.name}, Track Id: {frame.track_id}")
    canvas = FigureCanvasTkAgg(fig, master=frame)

    height_ratios = [12, 6, 1] + [1 for _ in range(len(frame.audible_ranges))] + [1]
    grid = GridSpec(ncols=2, nrows=4+len(frame.audible_ranges), figure=fig,
                    width_ratios=[10,1], height_ratios=height_ratios)

    map_ax = fig.add_subplot(grid[0, 0])
    spectro_ax = fig.add_subplot(grid[1, 0])
    cue_ax = fig.add_subplot(grid[2, 0])
    slider_axes = []
    rm_range_button_axes = []
    for i in range(len(frame.audible_ranges)):
        slider_axes.append(fig.add_subplot(grid[3+i, 0]))
        rm_range_button_axes.append(fig.add_subplot(grid[3+i, 1]))
    new_range_ax = fig.add_subplot(grid[grid.nrows-1, 0])

    # assign member variables so event handlers and AudibleRangeUI can see these
    frame.fig = fig
    frame.canvas = canvas
    frame.bg = None
    frame.slider_axes = slider_axes
    frame.spectro_ax = spectro_ax
    frame.map_ax = map_ax

    # --------------------------------- Plot Track --------------------------------- #

    # Display the study area, track points, spline points, closest point, and microphone
    map_ax.contour(
        frame.dem_x,
        frame.dem_y,
        frame.dem_z,
        levels=8,
        colors="lightgray",
        linewidths=1,
        zorder=0
    )
    frame.master.study_area.geometry.boundary.plot(
        label='study area',
        ax=map_ax,
        ls="--",
        lw=0.5,
        color="blue"
    )
    frame.spline.plot(
        label='interpolated spline point',
        ax=map_ax,
        color="grey",
        zorder=1,
        markersize=0.1,
        alpha=0.5
    )
    frame.points.plot(
        label='track point',
        ax=map_ax,
        color="blue",
        zorder=1,
        markersize=2,
    )
    frame.closest_point.plot(
        label='closest point',
        ax=map_ax,
        color="lime",
        zorder=1,
        markersize=32,
        edgecolor="blue"
    )

    # Point indicating where on the map corresponds to a certain spot in the spectrogram,
    # when the user hovers their mouse. Init at the microphone, but invisible
    frame.map_mousehover_point = map_ax.plot(
        frame.master.mic.x,
        frame.master.mic.y,
        ms=10,
        marker="o",
        markerfacecolor="none",
        markeredgecolor="black",
        zorder=3,
        animated=True  # don't show until explicitly told to, for use with blitting
    )[0]
    frame.map_mousehover_point.set_visible(False)

    map_ax.plot(
        frame.master.mic.x,
        frame.master.mic.y,
        label='microphone',
        ls="",
        marker="x",
        ms=7,
        color="magenta",
        zorder=10
    )

    # Glean the spatial extent of the points. This will result in a square map.
    xmin, ymin, xmax, ymax = frame.master.study_area.total_bounds
    map_ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    map_ax.set_aspect((xmax - xmin) / (ymax - ymin))
    map_ax.tick_params(axis='both', labelsize=6)
    map_ax.ticklabel_format(style='plain')  # disable scientific notation
    map_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # --------------------------------- Plot Spectrogram --------------------------------- #

    spectro_ax.imshow(
        frame.spectro.T,
        origin="lower",
        aspect=(frame.x_lims[-1] - frame.x_lims[0]) / (8 * 33),
        cmap="plasma",
        extent=[frame.x_lims[0], frame.x_lims[-1], 0, frame.spectro.shape[1]],
        interpolation=None,
        vmin=-10,
        vmax=80
    )

    # add in the time of closest approach in red.
    spectro_ax.axvline(
        date2num(frame.closest_time),
        alpha=0.7,
        color="lime",
        zorder=2,
        linewidth=3,
        label='closest track point'
    )

    spectro_ax.legend(bbox_to_anchor=(0.25, 1.4))
    spectro_ax.set_yticks(np.arange(frame.spectro.shape[1])[::6])
    spectro_ax.set_yticklabels(frame.spectro.columns.astype('float')[::6])
    spectro_ax.set_ylabel("Freq. (Hz)", labelpad=15)
    spectro_ax.axhline(8, lw=1.0, color="white", ls=":", alpha=0.4, zorder=200)
    spectro_ax.xaxis_date()  # tell matplotlib that the numeric axis should be formatted as dates
    spectro_ax.xaxis.set_major_formatter(DateFormatter("%b-%d\n%H:%M"))  # tidy them!

    # --------------------------------- Plot Expected Lp Cue Signal --------------------------------- #

    ys = np.zeros(len(frame.spline["time_audible"]))
    cue_ax.scatter(frame.spline["time_audible"], ys, c=frame.spline["Lp_est"], cmap="plasma")
    cue_ax.set_xlim(frame.x_lims[0], frame.x_lims[-1])
    cue_ax.set_ylabel("Sound Level Cue", rotation=0, ha='right', va='center')
    # remove extra plot bits
    for spine in cue_ax.spines.values():
        spine.set_visible(False)
    cue_ax.set_xticks([])
    cue_ax.set_yticks([])

    # --------------------------------- UI for Selecting Audible Ranges --------------------------------- #

    # note - we need to save references to the AudibleRangeUI objects, so the garbage collector doesn't trash them
    # and also - so we can update them in self.update_plot()
    frame.audible_range_uis = []
    frame.rm_range_buttons = []  # same for storing these references

    for i in range(len(frame.audible_ranges)):

        label = f"Audible Extent {i + 1}"
        ui = AudibleRangeUI(frame, i, label)
        frame.audible_range_uis.append(ui)

        rm_button = Button(rm_range_button_axes[i], "Remove")
        rm_button.on_clicked(partial(remove_audible_range, frame, i))  # use functools.partial to bind the current value of i, since i will change in the for loop
        frame.rm_range_buttons.append(rm_button)
        

    # --------------------------------- New Slider Button --------------------------------- #

    frame.new_slider_button = Button(new_range_ax, "Add Sound Event")
    frame.new_slider_button.on_clicked(lambda _: new_audible_range(frame, _))

    # --------------------------------- Update Track Labels and Event Handling --------------------------------- #

    type_line, name_line = _track_label_type_lines(
        frame.master.track_source,
        frame.aircraft_type,
        getattr(frame, "vessel_name", None),
    )
    altitude_line = _track_altitude_summary(
        frame.points,
        frame.closest_point.geometry.iat[0],
    )

    frame.track_label.config(text=f"Microphone: {frame.master.mic.name}\n\n" + \
                                 f"Track Id: {frame.track_id}\n" + \
                                 (f"{frame.aircraft_help_text}\n" if frame.aircraft_help_text is not None else "") + \
                                 name_line + type_line + \
                                 altitude_line + \
                                 f"\nAnnotated: {frame.track_annotated}" + \
                                 f"\nValid: {frame.valid}")
    frame.progress_label.config(text=f"{frame.i+1}/{frame.master.tracks.track_id.nunique()}")
    frame.submit_button.config(command=lambda: frame._store_annotation(frame.track_id, frame.spline, frame.audible_ranges), state=tk.NORMAL)
    frame.unknown_button.config(command=lambda: frame._store_annotation(frame.track_id, frame.spline, valid=False), state=tk.NORMAL)
    canvas.mpl_connect("motion_notify_event", lambda event: on_mouse_move(frame, event))
    canvas.mpl_connect("button_press_event", lambda event: on_mouse_down(frame, event))
    canvas.mpl_connect("draw_event", lambda event: on_draw(frame, event))

    # --------------------------------- Show Plot --------------------------------- #

    frame.fig_widget = canvas.get_tk_widget()
    frame.fig_widget.grid(row=0, column=0, sticky='nsew', rowspan=100)  # large rowspan so we don't have to update it if adding new rows

    # need to manually call canvas.draw() so that the correct background layout gets saved by self.on_draw()
    canvas.draw()


def on_draw(frame: "_GroundTruthingFrame", event: Any = None) -> None:
    frame.bg = frame.canvas.copy_from_bbox(frame.fig.bbox)
    update_plot(frame)

def update_plot(frame: "_GroundTruthingFrame") -> None:
    if frame.bg is None:
        return
    frame.canvas.restore_region(frame.bg)

    # mouse hover track point
    frame.map_ax.draw_artist(frame.map_mousehover_point)

    # audible range UIs
    for ui in frame.audible_range_uis:
        frame.map_ax.draw_artist(ui.highlight)
        frame.spectro_ax.draw_artist(ui.lower_limit_line)
        frame.spectro_ax.draw_artist(ui.upper_limit_line)
        ui.slider.draw()
    
    frame.canvas.blit(frame.fig.bbox)
    frame.canvas.flush_events()

def on_mouse_down(frame: "_GroundTruthingFrame", event: Any) -> None:
    if event.button == 1 and event.inaxes in frame.slider_axes:
        update_plot(frame)

def on_mouse_move(frame: "_GroundTruthingFrame", event: Any) -> None:
    if event.inaxes == frame.spectro_ax or event.inaxes in frame.slider_axes:
        dt = num2date(event.xdata).replace(tzinfo=None)

        # get closest spline point to the mouse position
        closest_idx = (frame.spline["time_audible"] - dt).abs().idxmin()
        closest_pt = frame.spline.loc[closest_idx]

        frame.time_label.config(text=_cursor_status_text(
            dt,
            _point_altitude_m(closest_pt.geometry),
        ))

        # if the mouse is close enough to a point, display the marker on the map
        if abs(closest_pt["time_audible"] - dt) > frame.typical_t_diff:
            frame.map_mousehover_point.set_visible(False)
        else:
            
            frame.map_mousehover_point.set_data(
                [closest_pt.geometry.x], [closest_pt.geometry.y])
            frame.map_mousehover_point.set_visible(True)

        update_plot(frame)
        

def new_audible_range(frame: "_GroundTruthingFrame", _: Any) -> None:
    """
    Add a new audible range. Since we can't change the layout of axes after making them,
    we need to clear the current plot and remake it, taking into account the new audible range.
    """
    frame.audible_ranges.append([
        num2date(frame.lower_limit_start),
        num2date(frame.upper_limit_start)
    ])
    build_plot(frame)

def remove_audible_range(frame: "_GroundTruthingFrame", i: int, _: Any) -> None:
    """Remove a certain audible range. See new_audible_range() for why we have to replot the figure."""
    del frame.audible_ranges[i]
    build_plot(frame)


class AudibleRangeUI():
    """Class to manage the various UI components of an audible range"""
    def __init__(self, gt_frame: "_GroundTruthingFrame", i: int, label: str) -> None:
        self.range_bounds = gt_frame.audible_ranges[i]  # where we store the range limits, part of the _GroundTruthingFrame.audible_ranges list
        self.slider_ax = gt_frame.slider_axes[i]
        self.gt_frame = gt_frame
        # aliases for easier code
        self.spline = self.gt_frame.spline
        self.bg = self.gt_frame.bg
        self.map_ax = self.gt_frame.map_ax
        self.spectro_ax = self.gt_frame.spectro_ax
        self.fig = self.gt_frame.fig
        self.canvas = self.gt_frame.canvas

        low_init = date2num(self.range_bounds[0])
        high_init = date2num(self.range_bounds[1])

        self.highlight = gt_frame.map_ax.plot(
            gt_frame.spline.geometry.x,
            gt_frame.spline.geometry.y,
            lw=8,
            color='deepskyblue',
            ls='-',
            zorder=1,
            alpha=0.4,
            animated=True
        )[0]
        self.lower_limit_line = gt_frame.spectro_ax.axvline(
            low_init,
            ls="--",
            alpha=0.7,
            color="white",
            zorder=2,
            linewidth=1,
            animated=True
        )
        self.upper_limit_line = gt_frame.spectro_ax.axvline(
            high_init,
            ls="--",
            alpha=0.7,
            color="white",
            zorder=2,
            linewidth=1,
            animated=True
        )

        self.slider = FastRangeSlider(
            self.slider_ax,
            label=label,
            valmin=gt_frame.x_lims[0],
            valmax=gt_frame.x_lims[-1],
            valinit=[low_init, high_init]
        )
        self.slider.drawon = False  # don't call draw_idle() automatically, I want to handle UI updates manually using blitting
        self.slider.valtext.set_visible(False)  # Turn off range slider value label.
        self.slider.on_changed(self._slider_update)
        self._slider_update([low_init, high_init])

    def _slider_update(self, val: list[float]) -> None:
        """
        Update spline highlight and spectrogram lines based on slider values.

        Parameters
        ----------
        val : list[float]
            A two item list with the [min, max] values of the range slider.
        """
        if self.bg is not None:
            self.canvas.restore_region(self.bg)
            
        lower_t, upper_t = val

        # update bounds - this gets propagated to _GroundTruthingFrame.audible_ranges, because self.range bounds list is inside that list
        self.range_bounds[0] = num2date(lower_t)
        self.range_bounds[1] = num2date(upper_t)

        # Update the vertical lines on the spectrogram
        self.lower_limit_line.set_xdata([lower_t, lower_t])
        self.upper_limit_line.set_xdata([upper_t, upper_t])

        # Highlight the section of the track that falls within the date window
        # NOTE: .replace(tzinfo) is required to prevent errors from comparing tz-naive again tz-aware datetimes
        subset_mask = np.all([
            self.gt_frame.spline.time_audible >= num2date(lower_t).replace(tzinfo=None),
            self.gt_frame.spline.time_audible <= num2date(upper_t).replace(tzinfo=None)
        ], axis=0)
        subset = self.gt_frame.spline.loc[subset_mask]
        self.highlight.set_data(subset.geometry.x, subset.geometry.y)
        
        # No need to call update_plot() here, any updates will be paired with a mousemove or mousedown
        # and the mousemove/mousedown event handlers call update_plot().
        # This is also nice because we avoid having to deal with duplicate update_plot() calls
