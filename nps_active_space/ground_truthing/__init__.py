import datetime as dt
import math
import tkinter as tk
import traceback
from abc import ABC
from tkinter import filedialog, messagebox
from typing import Any, Optional, Type

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import date2num, DateFormatter
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RangeSlider
from shapely.geometry import Point

from nps_active_space import _ACTIVE_SPACE_DIR
from nps_active_space.utils import audible_time_delay, interpolate_spline, Microphone, Nvspl, Tracks


_app = None


def launch(*args, **kwargs):
    """A wrapper function to launch the ground truthing application."""
    global _app
    _app = _App(*args, **kwargs)
    _app.mainloop()


class _AppFrame(ABC, tk.Frame):
    """Abstract base class for all application frames."""

    def __init__(self, master):
        super().__init__(master)
        self.master = master


class _App(tk.Tk):

    def __init__(self, outfile: str, tracks: Tracks, nvspl: Nvspl, mic: Microphone, crs,
                 site_shp: Optional[gpd.GeoDataFrame] = None, clip: Optional[bool] = False):
        super().__init__()

        self.outfile = outfile
        self.crs = crs
        self.mic = mic.to_crs(crs)
        self.site_shp = site_shp.to_crs(crs)
        self.tracks = gpd.clip(tracks.to_crs(crs), self.site_shp) if clip else tracks.to_crs(crs)
        self.nvspl = nvspl

        # Set app features.
        self.title('NPS Active Space: Ground Truthing Module')
        self.iconbitmap(f"{_ACTIVE_SPACE_DIR}/img/flat-four-color.ico")
        self.geometry('1200x600')

        # Create app menu.
        self.menu = tk.Menu(self)
        self.file_menu = tk.Menu(self.menu, tearoff=False)
        self.file_menu.add_command(label='Save...', command=self._save)
        self.file_menu.add_command(label='Plot...', command=self._plot)
        self.file_menu.add_separator()
        self.file_menu.add_command(label='Exit', command=self._close)
        self.menu.add_cascade(label='File', menu=self.file_menu)
        self.config(menu=self.menu)

        # Create the starting state.
        self.annotations = pd.DataFrame(columns={'_id', 'valid', 'audible', 'starttime', 'endtime'})
        self._saved = True
        self._frame = None

        self.switch_frame(_WelcomeFrame)

    def run(self):
        """Run the main application frame."""
        self.protocol("WM_DELETE_WINDOW", self._close)
        self.switch_frame(_GroundTruthingFrame) # TODO: try accept here with traceback

    def switch_frame(self, frame_class: Type[_AppFrame]):
        """
        Switch the frame that is being displayed in the application window.

        Parameters
        ----------
        frame_class : _AppFrame
            The frame class of the frame to display.
        """
        new_frame = frame_class(self)
        new_frame['bg'] = 'ivory2'
        if self._frame is not None:
            self._frame.destroy()
        self._frame = new_frame
        self._frame.pack(expand=True, anchor='nw', fill=tk.BOTH)

    def add_annotation(self, id_: Any, valid: bool, audible: bool, starttime: str, endtime: str):
        """
        Add a new annotation.

        Parameters # TODO
        ----------
        id_
        valid: bool
        audible : bool, default None
        starttime: str, default None
        endtime: str, default None
        """
        new_record = pd.DataFrame.from_records([
            {'_id': id_,
             'valid': valid,
             'audible': audible,
             'starttime': starttime,
             'endtime': endtime}
        ])
        self.annotations = pd.concat([self.annotations, new_record], ignore_index=True)

    def load_annotations(self, filename):
        self.annotations = pd.read_csv(filename, usecols=self.annotations.columns)
        print(self.annotations)
    
    def unsave(self):
        self._saved = False

    def _close(self):
        """
        A function to safely close the application. If the user has unsaved changes, they will be warned and asked
        if they would like to proceed before closing the application.
        """
        if self._saved:
            result = 'yes'
        else:
            result = tk.messagebox.askquestion(
                title='Exit',
                message='Are you sure you want to exit without saving?',
                icon='warning',
                default='no'
            )
        if result == 'yes':
            plt.close('all')
            self.destroy()

    def _save(self):
        """Save current annotations to the output file."""
        self.annotations.to_csv(self.outfile, mode='w')
        self._saved = True

    def _plot(self):
        # TODO

        # plotting CRS, lat/long looks nicer here...

        fig, ax = plt.subplots(1, 1, figsize=(6, 9))

        # # show basemap in lat/long # TODO
        # with rasterio.open(r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\DENA Park Brochure Map wi.tif") as raster:
        #     rasterio.plot.show(raster, ax=ax, alpha=0.6)

        # plot study area bounding box
        study_area = self.site_shp.to_crs('epsg:4326')
        study_area.geometry.boundary.plot(ax=ax, ls="--", color="navy")

        # # # plot microphone position
        # ax.plot(site_info["long"], site_info["lat"], ls="", marker="*", ms=6, color="black",
        #         zorder=3, label="microphone\n" + u + s + str(y))

        # # audible vs inaudible points:
        # valid_points = valid_points.to_crs(plot_CRS)
        # valid_points.loc[valid_points.audible].plot(ax=ax, color='deepskyblue', alpha=0.05, markersize=3, zorder=3,
        #                                             label="Audible points")
        # valid_points.loc[~valid_points.audible].plot(ax=ax, color='red', alpha=0.05, markersize=3, zorder=2,
        #                                              label="Inaudible points")
        #
        # # create a legend at the bottom center
        # leg = plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), markerscale=2)
        # for lh in leg.legendHandles:
        #     lh.set_alpha(1)  # set legend handles to opaque for visibility
        #
        # xmin, ymin, xmax, ymax = study_area.total_bounds
        # pad = np.array([(xmax - xmin) * 0.1, (ymax - ymin) * 0.1])
        #
        # # this will result in a square map
        # ax.set(xlim=(xmin - pad[0], xmax + pad[0]), ylim=(ymin - pad[1], ymax + pad[1]))
        # ax.set_title("Annotated audible overflights segments")
        # ax.tick_params(axis='both', labelsize=6)
        #
        # plt.show()


class _WelcomeFrame(_AppFrame):
    """
    The opening frame for the Ground Truthing application welcome the user.

    Parameters
    ----------
    master : tk.Tk
        The tkinter window this frame will be shown in.
    """
    def __init__(self, master):
        super().__init__(master)

        # Define widgets.
        frame_label = tk.Label(
            self,
            text='Welcome to the NPS Active Space Project\nGround Truthing Module!',
            font=('Avenir', 20, 'bold'),
            bg='ivory2'
        )
        continue_button = tk.Button(
            self,
            text='Continue >>',
            width=20,
            font=('Avenir', 8),
            bg='ivory2',
            command=lambda: self.master.switch_frame(_AnnotationLoadFrame)
        )

        # TODO: logo??

        # Place widgets.
        frame_label.place(relx=0.5, rely=0.45, anchor='center')
        continue_button.place(relx=0.9, rely=0.9, anchor='center')


class _AnnotationLoadFrame(_AppFrame):
    """
    A frame to allow the user to decide if they would like previously saved annotations to be loaded.

    Parameters
    ----------
    master : tk.Tk
        The tkinter window this frame will be shown in.
    """
    def __init__(self, master):
        super().__init__(master)

        # Define vars.
        self.load_annotations = tk.BooleanVar(value=False)
        self.annotation_filename = tk.StringVar(value='')

        # Define widgets.
        self.select_file_button = tk.Button(
            self,
            text='Select File',
            bg='ivory2',
            command=lambda: self._select_file()
        )
        self.select_file_label = tk.Label(
            self,
            bg='ivory2'
        )

        self.frame_label = tk.Label(
            self,
            text='Would you like to check for saved annotations?',
            font=('Avenir', 14, 'bold'),
            bg='ivory2'
        )
        self.select_label = tk.Label(
            self,
            text='Select an option:',
            font=('Avenir', 11, 'italic'),
            bg='ivory2'
        )
        self.yes_button = tk.Radiobutton(
            self,
            text='Yes, load annotations from file.',
            font=('Avenir', 10),
            value=True,
            variable=self.load_annotations,
            bg='ivory2',
            command=lambda: self.select_file_button.place(relx=0.6, rely=0.48, anchor='w')
        )
        self.no_button = tk.Radiobutton(
            self,
            text='No, do not load prior annotations.',
            font=('Avenir', 10),
            value=False,
            variable=self.load_annotations,
            bg='ivory2',
            command=lambda: self._clear()
        )
        self.continue_button = tk.Button(
            self,
            text='Continue >>',
            width=20,
            font=('Avenir', 8),
            bg='ivory2',
            command=lambda: self._option_selected()
        )

        # Place widgets.
        self.frame_label.place(relx=0.5, rely=0.3, anchor='center')
        self.select_label.place(relx=0.5, rely=0.4, anchor='center')
        self.yes_button.place(relx=0.41, rely=0.48, anchor='w')
        self.no_button.place(relx=0.41, rely=0.53, anchor='w')
        self.continue_button.place(relx=0.9, rely=0.9, anchor='center')

    def _clear(self):
        self.select_file_button.place_forget()
        self.select_file_label.place_forget()
        self.select_file_label.config(text='')
        self.annotation_filename.set('')

    def _select_file(self):
        filetypes = (('csv files', '*.csv'),)
        filename = filedialog.askopenfilename(
            title='Open file',
            initialdir='/',
            filetypes=filetypes
        )
        if filename:
            self.annotation_filename.set(filename)
            self.select_file_label.config(text=f"{filename[:50]}...")
            self.select_file_label.place(relx=0.66, rely=0.48, anchor='w')

    def _option_selected(self):
        """
        Determine what app frame should be shown next depending on if the user would like saved annotations to
        be loaded or not.
        """
        if self.load_annotations.get() is False:
            self.master.switch_frame(_InstructionsFrame)

        elif self.load_annotations.get() is True and self.annotation_filename.get():
            self.master.load_annotations(self.annotation_filename.get())
            self.master.switch_frame(_InstructionsFrame)


class _InstructionsFrame(_AppFrame):
    """
    Frame describing how the module works to the user.

    Parameters
    ----------
    master : tk.Tk
        The tkinter window this frame will be shown in.
    """
    def __init__(self, master):
        super().__init__(master)

        # Define widgets.
        frame_label = tk.Label(
            self,
            text='Instructions:',
            font=('Avenir', 14, 'bold'),
            bg='ivory2'
        )
        instructions = tk.Label(
            self,
            text='Use the range slider to adjust what section of each track is audible, inaudible, or unknown.',
            font=('Avenir', 12),
            bg='ivory2'
        )
        save_reminder = tk.Label(
            self,
            text='As always, make sure to save intermittently!',
            font=('Avenir', 12),
            bg='ivory2'
        )
        start_button = tk.Button(
            self,
            text='Start',
            font=('Avenir', 8),
            width=20,
            bg='ivory2',
            command=lambda: self.master.run()
        )

        # Place widgets.
        frame_label.place(relx=0.5, rely=0.35, anchor='center')
        instructions.place(relx=0.5, rely=0.45, anchor='center')
        save_reminder.place(relx=0.5, rely=0.5, anchor='center')
        start_button.place(relx=0.9, rely=0.9, anchor='center')


class _CompletionFrame(_AppFrame):
    """
    Frame to be shown after all tracks have been classified,

    Parameters
    ----------
        master : tk.Tk
            The tkinter window this frame will be shown in.
    """
    def __init__(self, master):
        super().__init__(master)

        # Define widgets.
        frame_label = tk.Label(
            self,
            text='Annotations Completed!',
            font=('Avenir', 20, 'bold'),
            bg='ivory2',
        )

        # Place widgets.
        frame_label.place(relx=0.5, rely=0.45, anchor='center')


class _GroundTruthingFrame(_AppFrame):

    def __init__(self, master):
        super().__init__(master)

        self.data = iter(self.master.tracks.groupby(by='track_id'))
        self.canvas = None
        self.i = 1

        # Define widgets.
        self.audible_button = tk.Button(
            self,
            text='Audible >>',
            bg='green',
            fg='white',
            width=10,
            font=('Avenir', 12, 'bold')
        )
        self.inaudible_button = tk.Button(
            self,
            text='Inaudible >>',
            bg='red',
            fg='white',
            width=10,
            font=('Avenir', 12, 'bold')
        )
        self.unknown_button = tk.Button(
            self,
            text='Unknown >>',
            bg='yellow',
            fg='black',
            width=10,
            font=('Avenir', 12, 'bold')
        )
        self.progress_label = tk.Label(
            self,
            bg='ivory2'
        )

        # Place widgets.
        self.grid_columnconfigure(0, weight=5)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=2)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=2)
        self.audible_button.grid(row=1, column=1, sticky='n')
        self.inaudible_button.grid(row=1, column=1)
        self.unknown_button.grid(row=1, column=1, sticky='s')
        self.progress_label.grid(row=0, column=1, sticky='ne', padx=10, pady=5)

        self.next()

    def next(self):

        try:
            idx, points = next(self.data)
        except StopIteration:
            self.master.switch_frames(_CompletionFrame)

        if points.shape[0] >= 3:
            pass
            # TODO

        # TODO: Load prior annotation
        self._build_plot(idx, points)
        self.i += 1

    def _build_plot(self, idx, points):

        points.sort_values(by='point_dt', ascending=True, inplace=True)
        start_time = points.point_dt.iat[0]
        end_time = points.point_dt.iat[-1]
        duration_s = (end_time - start_time).total_seconds()

        spline = interpolate_spline(points)
        spline = audible_time_delay(spline, 'point_dt', Point(self.master.mic.x, self.master.mic.y, self.master.mic.z))
        closest_time = spline.loc[spline.distance_to_target.idxmin()]['point_dt']

        fig = plt.figure(constrained_layout=True)
        spec = GridSpec(ncols=1, nrows=10, figure=fig)
        ax1 = fig.add_subplot(spec[0:6, 0])
        ax2 = fig.add_subplot(spec[6:9, 0])
        ax3 = fig.add_subplot(spec[9, 0])

        # ---- Plot Track ---- #
        # we could plot with geopandas, but using pyplot will allow for a dynamic plot that reacts to the range slider input
        # that is, xy_UTM is a copy of the flight_spline geometry for plotting ONLY
        xy_UTM = spline.geometry.bounds.iloc[:, 0:2]
        xy_UTM['time_audible'] = spline['time_audible']  # when we can actually HEAR the noise, mic time

        #`highlight` is a dynamic line object that will update with the range slider
        highlight, = ax1.plot(xy_UTM.minx, xy_UTM.miny, lw=5, color='deepskyblue', ls='-', zorder=1, alpha=0.4)

        self.master.site_shp.geometry.boundary.plot(ax=ax1, ls="--", lw=0.5, color="blue", label='study area')
        points.plot(ax=ax1, color="blue", zorder=1, markersize=3, label='track point')
        ax1.plot(self.master.mic.x, self.master.mic.y, ls="", marker="x", ms=5, color="red",
                 zorder=10, label='microphone')
        # TODO: rasterio.plot.show(raster, ax=ax[0], alpha=0.6)  # show basemap? - makes annotating much slower

        # glean the spatial extent of the points
        xmin, ymin, xmax, ymax = self.master.site_shp.total_bounds

        # this will result in a square map
        ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        ax1.set_aspect((xmax - xmin) / (ymax - ymin))
        ax1.tick_params(axis='both', labelsize=6)
        ax1.ticklabel_format(style='plain')  # disable scientific notation
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # ---- Plot Slider ---- #

        slider = RangeSlider(ax3, "Audible Window [flight duration, s]", 0, duration_s)
        # slider.set_val([-60, 60] + (closest_time - start_time).total_seconds())
        self._slider_update(None)  # update highlight
        slider.on_changed(self._slider_update)

        # ---- Plot Spectrogram ---- #

        # TODO: make sure spectrogram exists
        spectrogram = self.master.nvspl.loc[str(start_time):str(end_time), '12.5':'20000']  # Time slice of NVSPL record
        x_lims = date2num(spectrogram.index)  # convert the NVSPL's nice datetime axis to numbers
        aspect = (x_lims[-1] - x_lims[0]) / (8 * 33)

        ax2.imshow(spectrogram.T, origin="lower", aspect=aspect, cmap="plasma",
                   extent=[x_lims[0], x_lims[-1], 0, spectrogram.shape[1]],
                   interpolation=None, vmin=-10, vmax=80)

        # add in the time of closest approach in red
        ax2.axvline(date2num(closest_time), alpha=0.7, color="red", zorder=2, linewidth=3)

        # Create the moving vertical lines on the histogram with axvline()
        lower_limit_line = ax2.axvline(date2num(start_time + dt.timedelta(seconds=slider.val[0])),
                                         ls="--", alpha=0.7, color="white", zorder=2, linewidth=1)
        upper_limit_line = ax2.axvline(date2num(start_time + dt.timedelta(seconds=slider.val[1])),
                                         ls="--", alpha=0.7, color="white", zorder=2, linewidth=1)

        ax2.set_yticks(np.arange(spectrogram.shape[1])[::6])
        ax2.set_yticklabels(spectrogram.columns.astype('float')[::6])
        ax2.set_ylabel("Freq. (Hz)", labelpad=15)

        # tell matplotlib that the numeric axis should be formatted as dates
        ax2.xaxis_date()
        ax2.xaxis.set_major_formatter(DateFormatter("%b-%d\n%H:%M"))  # tidy them!

        # ---- Show Plot ---- #
        # TODO: add times
        self.audible_button.config(command=lambda: self._click(idx, True, True), state=tk.NORMAL)
        self.inaudible_button.config(command=lambda: self._click(idx, True, False),  state=tk.NORMAL)
        self.unknown_button.config(command=lambda: self._click(idx, False),  state=tk.NORMAL)
        self.progress_label.config(text=f"{self.i}/{self.master.tracks.track_id.nunique()}") # todo: wrong?

        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew', rowspan=3)

    def _click(self, id_: Any, valid: bool, audible: Optional[bool] = None,
              starttime: Optional[str] = None, endtime: Optional[bool] = None):

        # Disable buttons
        self.audible_button.config(state=tk.DISABLED)
        self.inaudible_button.config(state=tk.DISABLED)
        self.unknown_button.config(state=tk.DISABLED)

        # Record the new annotation and set the window as unsaved.
        self.master.add_annotation(id_, valid, audible, starttime, endtime)
        self.master.unsave()

        # Close the plot and canvas to clear the window for the next one.
        plt.close()
        if self.canvas:
            self.canvas.get_tk_widget().destroy()

        self.next()

    @staticmethod
    def _slider_update(val):
        pass

        # # The val passed to a callback by the RangeSlider will
        # # be a tuple of (min, max)
        #
        # # read the slider values
        # lower_t = start + dt.timedelta(seconds=slider.val[0])
        # upper_t = start + dt.timedelta(seconds=slider.val[1])
        # lower_date = date2num(lower_t)
        # upper_date = date2num(upper_t)
        #
        # # update the vertical lines on the spectrogram
        # lower_limit_line.set_xdata([lower_date, lower_date])
        # upper_limit_line.set_xdata([upper_date, upper_date])
        #
        # subset = xy_UTM.loc[np.all([flight_spline.time_audible >= lower_t,
        #                             flight_spline.time_audible <= upper_t],
        #                            axis=0)]  # subset of entire flight spline
        # highlight.set_data(subset.minx, subset.miny)
        #
        # # Redraw the figure to ensure it updates
        # fig.canvas.draw_idle()
        #
        #
        # slider_update(None)  # update highlight
        # slider.on_changed(slider_update)
        #
        # # add radio buttons (selection) too:
        # ax_buttons = plt.axes([0.125, 0.02, 0.1, 0.18])
        # radio_buttons = RadioButtons(ax_buttons, ['Audible', 'Inaudible', 'Unknown', 'Exit'], active=None)

