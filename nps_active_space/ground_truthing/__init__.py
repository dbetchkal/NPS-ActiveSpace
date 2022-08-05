import datetime as dt
import tkinter as tk
import traceback
from abc import ABC
from tkinter import filedialog, messagebox
from typing import Any, List, Optional, Type, TYPE_CHECKING

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import date2num, DateFormatter, num2date
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RangeSlider
from PIL import Image, ImageTk
from shapely.geometry import Point

from nps_active_space import _ACTIVE_SPACE_DIR
from nps_active_space.utils import audible_time_delay, interpolate_spline

if TYPE_CHECKING:
    from nps_active_space.utils import Microphone, Nvspl, Tracks


_app = None


def launch(*args, **kwargs):
    """A wrapper function to launch the ground truthing application."""
    global _app
    _app = _App(*args, **kwargs)
    _app.mainloop()


class _AppFrame(ABC, tk.Frame):
    """
    Abstract base class for all application frames.

    Parameters
    ----------
    master : tk.Tk
        A tkinter app instance that will display the frame.
    """
    def __init__(self, master: tk.Tk):
        super().__init__(master)
        self.master = master


class _App(tk.Tk):
    """
    Main ground truthing application window.

    Parameters
    ----------
    outfile : str
        Absolute path to the csv file where annotation should be output. Format: /path/to/file.csv
    mic : Microphone
        A Microphone object of the microphone deployment to be used for ground truthing.
    nvspl : Nvspl
        An Nvspl object of sound data record at the input microphone locations.
    tracks : Tracks
        a Tracks object of points to classify as audible, inaudible, or unknown from the microphone location.
    crs : str
        The projected coordinate system to be used for the Tracks, study area, and microphone.
        Format of 'epsg:XXXX...', E.g. 'epsg:32632'
    study_area : gpd.GeoDataFrame
        A gpd.GeoDataFrame of polygon(s) that make up the study area.
    clip : bool, default False
        If True, clip the Tracks to the study area.
    """
    def __init__(self, outfile: str, mic: 'Microphone', nvspl: 'Nvspl', tracks: 'Tracks',
                 crs: str, study_area: gpd.GeoDataFrame, clip: bool = False):
        super().__init__()

        self.outfile = outfile
        self.crs = crs
        self.mic = mic.to_crs(crs)
        self.study_area = study_area.to_crs(crs)
        self.tracks = gpd.clip(tracks.to_crs(crs), self.study_area) if clip else tracks.to_crs(crs)
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

        # Create the application starting state.
        self.annotations = pd.DataFrame(columns={'_id', 'valid', 'audible', 'starttime', 'endtime'})
        self._saved = True
        self._frame = None

        self.switch_frame(_WelcomeFrame)

    def run(self):
        """Run the main application frame."""
        self.protocol("WM_DELETE_WINDOW", self._close)
        self.switch_frame(_GroundTruthingFrame)

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

    def add_annotation(self, id_: Any, valid: bool, audible: bool = False,
                       starttime: Optional[str] = None, endtime: Optional[str] = None):
        """
        Add a new track audibility annotation.

        Parameters
        ----------
        id_ : Any
            The track unique identifier.
        valid : bool
            If the track was valid.
        audible : bool
            If the track was valid, was it audible.
        starttime : str, default None
            If the track was audible, when does audibility start.
        endtime : str, default None
            If the track was audible, when does audibility end.
        """
        new_record = pd.DataFrame.from_records([
            {'_id': id_,
             'valid': valid,
             'audible': audible,
             'starttime': starttime,
             'endtime': endtime}
        ])
        self.annotations = pd.concat([self.annotations, new_record], ignore_index=True)
        self._saved = False

    def load_annotations(self, filename: str):
        """
        Simple function to load existing annotations from a csv file.

        Parameters
        ----------
        filename : str
            Absolute path to csv file to load previous annotations from.
        """
        self.annotations = pd.read_csv(filename, usecols=self.annotations.columns)

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
        # TODO this function....

        # plotting CRS, lat/long looks nicer here...

        fig, ax = plt.subplots(1, 1, figsize=(6, 9))

        # # show basemap in lat/long # TODO
        # with rasterio.open(r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\DENA Park Brochure Map wi.tif") as raster:
        #     rasterio.plot.show(raster, ax=ax, alpha=0.6)

        # plot study area bounding box
        study_area = self.study_area.to_crs('epsg:4326')
        study_area.geometry.boundary.plot(ax=ax, ls="--", color="navy")

        # # plot microphone position
        ax.plot(
            self.mic.lon,
            self.mic.lat,
            ls="",
            marker="*",
            ms=6,
            color="black",
            zorder=3,
            label=self.mic.name
        )

        # audible vs inaudible points:
        valid_points = self.annotations[self.annotations.valid == True].copy().to_crs('epsg:4326')
        valid_points.loc[valid_points.audible].plot(ax=ax, color='deepskyblue', alpha=0.05, markersize=3, zorder=3,
                                                    label="Audible points")
        valid_points.loc[~valid_points.audible].plot(ax=ax, color='red', alpha=0.05, markersize=3, zorder=2,
                                                     label="Inaudible points")

        # create a legend at the bottom center
        leg = plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), markerscale=2)
        for lh in leg.legendHandles:
            lh.set_alpha(1)  # set legend handles to opaque for visibility

        xmin, ymin, xmax, ymax = study_area.total_bounds
        pad = np.array([(xmax - xmin) * 0.1, (ymax - ymin) * 0.1])

        # this will result in a square map
        ax.set(xlim=(xmin - pad[0], xmax + pad[0]), ylim=(ymin - pad[1], ymax + pad[1]))
        ax.set_title("Annotated audible overflights segments")
        ax.tick_params(axis='both', labelsize=6)

        plt.show()


class _WelcomeFrame(_AppFrame):
    """
    The opening frame for the Ground Truthing application that welcomes the user.

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
        im = Image.open(f"{_ACTIVE_SPACE_DIR}/img/flat-four-color.png").resize((138, 181))
        nps_logo = ImageTk.PhotoImage(im)
        label = tk.Label(self, image=nps_logo, bg='ivory2')
        label.image = nps_logo  # NOTE: This re-definition is required for windows machines.

        # Place widgets.
        label.place(relx=0.5, rely=0.3, anchor='center')
        frame_label.place(relx=0.5, rely=0.55, anchor='center')
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
        """Remove the Select File and related widgets if No option is selected."""
        self.select_file_button.place_forget()
        self.select_file_label.place_forget()
        self.select_file_label.config(text='')
        self.annotation_filename.set('')

    def _select_file(self):
        """Open File Dialog and save the selected file."""
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
        """If user wants to load existing annotations, load them before proceeding to the app instructions frame."""
        if self.load_annotations.get() is True and self.annotation_filename.get():
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
        back_button = tk.Button(
            self,
            text='<< Back',
            font=('Avenir', 8),
            width=20,
            bg='ivory2',
            command=lambda: self.master.switch_frame(_AnnotationLoadFrame)
        )

        # Place widgets.
        frame_label.place(relx=0.5, rely=0.35, anchor='center')
        instructions.place(relx=0.5, rely=0.45, anchor='center')
        save_reminder.place(relx=0.5, rely=0.5, anchor='center')
        start_button.place(relx=0.9, rely=0.9, anchor='center')
        back_button.place(relx=0.1, rely=0.9, anchor='center')


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
    """
    Main application frame that allows the users to mark the audibility of each track.

    Parameters
    ----------
        master : tk.Tk
            The tkinter window this frame will be shown in.
    """
    def __init__(self, master):
        super().__init__(master)

        # Set frame variables to starting values.
        self.data = iter(self.master.tracks.groupby(by='track_id'))
        self.slider = None
        self.i = 0

        # Define widgets.
        self.track_label = tk.Label(
            self,
            bg='ivory2',
            font=('Avenir', 10, 'bold')
        )
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
        self.track_label.grid(row=0, column=1, pady=20)
        self.audible_button.grid(row=1, column=1, sticky='n')
        self.inaudible_button.grid(row=1, column=1)
        self.unknown_button.grid(row=1, column=1, sticky='s')
        self.progress_label.grid(row=0, column=1, sticky='ne', padx=10, pady=5)

        mpl.rcParams['toolbar'] = 'None'    # Remove matplotlib toolbar

        self._next()

    def _next(self):
        """
        Get the next track ready for plotting.

        Notes
        -----
        There is a bit of weird recursion logic here. In order to only show one plot at a time, we need to use the
        next() method instead of for _ in iter() to loop through the track data. But, to go to the next track we
        need to call self._next() from itself. To eliminate unwanted effects from this necessary recursion, the
        if/elif/else statements are required.
        """
        try:
            idx, points = next(self.data)
            self.i += 1
            spectro = self.master.nvspl.loc[str(points.point_dt.iat[0]):str(points.point_dt.iat[-1]), '12.5':'20000']

            # If the track is already annotated, move on.
            if idx in self.master.annotations._id.values:
                self._next()

            # If the track does not have enough points for processing, mark it as invalid and move on.
            elif points.shape[0] < 3:
                tk.messagebox.showwarning(
                    title='Data Warning',
                    message=f"Track {idx} has fewer than 3 points and therefore cannot be processed. Skipping...",
                    icon='warning'
                )
                self._click(idx, False, False)

            # If there is no spectrogram data for the track, mark it as invalid and move on.
            elif spectro.empty:
                tk.messagebox.showwarning(
                    title='Data Warning',
                    message=f"Track {idx} has no accompanying spectrogram. Skipping...",
                    icon='warning'
                )
                self._click(idx, False, False)

            # If this is an un-annotated track with 3+ points and corresponding SPL data, display its plot.
            else:
                self._build_plot(idx, points, spectro)

        except StopIteration:
            self.master.switch_frame(_CompletionFrame)

    def _build_plot(self, idx: Any, points: 'Tracks', spectrogram: 'Nvspl'):
        """
        Build the matplotlib GridSpec plot for a track.

        Parameters
        ----------
        idx : Any
            The track id.
        points : Track
            The subset of Track points for the track id.
        spectrogram : Nvspl
            Nvspl data that aligns with the track points times.
        """
        def _slider_update(val: List):
            """
            Update spline highlight and spectrogram lines based on slider values.

            Parameters
            ----------
            val : List
                A two item list with the [min, max] values of the range slider.
            """
            lower_t = val[0]
            upper_t = val[1]

            # Update the vertical lines on the spectrogram
            lower_limit_line.set_xdata([lower_t, lower_t])
            upper_limit_line.set_xdata([upper_t, upper_t])

            # Highlight the section of the track that falls within the date window
            #
            # NOTE: .replace(tzinfo) is required to prevent errors from comparing tz-naive again tz-aware datetimes
            subset = spline.loc[np.all([spline.time_audible >= num2date(lower_t).replace(tzinfo=None),
                                        spline.time_audible <= num2date(upper_t).replace(tzinfo=None)],
                                       axis=0)]
            highlight.set_data(subset.geometry.x, subset.geometry.y)

            # todo: CHECK THIS, valid point spline?
            self.audible_button.config(
                state=tk.NORMAL,
                command=lambda: self._click(idx, True, True,  num2date(lower_t),  num2date(upper_t))
            )

            # Redraw the figure to ensure it updates
            fig.canvas.draw_idle()

        points.sort_values(by='point_dt', ascending=True, inplace=True)
        spline = interpolate_spline(points)
        spline = audible_time_delay(spline, 'point_dt', Point(self.master.mic.x, self.master.mic.y, self.master.mic.z))
        closest_point = spline[spline.distance_to_target == spline.distance_to_target.min()]
        closest_time = spline.loc[spline.distance_to_target.idxmin()]['time_audible']

        # ************************************ Build Plot ************************************#

        fig = plt.figure(figsize=(9, 5), constrained_layout=True)
        fig.canvas.manager.set_window_title(f"Microphone: {self.master.mic.name}, Track Id: {idx}")
        spec = GridSpec(ncols=1, nrows=10, figure=fig)
        ax1 = fig.add_subplot(spec[9, 0])
        ax2 = fig.add_subplot(spec[6:9, 0])
        ax3 = fig.add_subplot(spec[0:6, 0])

        # --------------------------------- Plot Slider --------------------------------- #

        x_lims = date2num(spectrogram.index)  # convert the NVSPL's nice datetime axis to numbers
        lower_limit_start = max(date2num(closest_time - dt.timedelta(seconds=60)), x_lims[0])
        upper_limit_start = min(date2num(closest_time + dt.timedelta(seconds=60)), x_lims[-1])

        self.slider = RangeSlider(
            ax1,
            label="Audible Window",
            valmin=x_lims[0],
            valmax=x_lims[-1],
            valinit=[lower_limit_start, upper_limit_start],
        )

        self.slider.valtext.set_visible(False)  # Turn off range slider value label.
        self.slider.on_changed(_slider_update)

        # --------------------------------- Plot Spectrogram --------------------------------- #

        ax2.imshow(
            spectrogram.T,
            origin="lower",
            aspect=(x_lims[-1] - x_lims[0]) / (8 * 33),
            cmap="plasma",
            extent=[x_lims[0], x_lims[-1], 0, spectrogram.shape[1]],
            interpolation=None,
            vmin=-10,
            vmax=80
        )

        # add in the time of closest approach in red.
        ax2.axvline(
            date2num(closest_time),
            alpha=0.7,
            color="red",
            zorder=2,
            linewidth=3,
            label='closest track point'
        )

        # Create the moving vertical lines on the histogram with axvline()
        lower_limit_line = ax2.axvline(
            lower_limit_start,
            ls="--",
            alpha=0.7,
            color="white",
            zorder=2,
            linewidth=1,
        )
        upper_limit_line = ax2.axvline(
            upper_limit_start,
            ls="--",
            alpha=0.7,
            color="white",
            zorder=2,
            linewidth=1
        )
        ax2.legend(bbox_to_anchor=(0.25, 1.4))
        ax2.set_yticks(np.arange(spectrogram.shape[1])[::6])
        ax2.set_yticklabels(spectrogram.columns.astype('float')[::6])
        ax2.set_ylabel("Freq. (Hz)", labelpad=15)
        ax2.xaxis_date()  # tell matplotlib that the numeric axis should be formatted as dates
        ax2.xaxis.set_major_formatter(DateFormatter("%b-%d\n%H:%M"))  # tidy them!

        # --------------------------------- Plot Track --------------------------------- #

        # Display the study area, track points, spline points, closest point, and microphone
        self.master.study_area.geometry.boundary.plot(
            label='study area',
            ax=ax3,
            ls="--",
            lw=0.5,
            color="blue"
        )
        spline.plot(
            label='interpolated spline point',
            ax=ax3,
            color="grey",
            zorder=1,
            markersize=0.5,
            alpha=0.1
        )
        points.plot(
            label='track point',
            ax=ax3,
            color="blue",
            zorder=1,
            markersize=3,
        )
        closest_point.plot(
            label='closest point',
            ax=ax3,
            color="red",
            zorder=1,
            markersize=3,
        )
        ax3.plot(
            self.master.mic.x,
            self.master.mic.y,
            label='microphone',
            ls="",
            marker="x",
            ms=7,
            color="magenta",
            zorder=10
        )

        # Highlight the points within the currently selected SPL region.
        subset = spline.loc[np.all([spline.time_audible >= num2date(lower_limit_start).replace(tzinfo=None),
                                    spline.time_audible <= num2date(upper_limit_start).replace(tzinfo=None)],
                                   axis=0)]
        highlight, = ax3.plot(
            subset.geometry.x,
            subset.geometry.y,
            lw=5,
            color='deepskyblue',
            ls='-',
            zorder=1,
            alpha=0.4
        )

        # Glean the spatial extent of the points. This will result in a square map.
        xmin, ymin, xmax, ymax = self.master.study_area.total_bounds
        ax3.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        ax3.set_aspect((xmax - xmin) / (ymax - ymin))
        ax3.tick_params(axis='both', labelsize=6)
        ax3.ticklabel_format(style='plain')  # disable scientific notation
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # --------------------------------- Show Plot --------------------------------- #

        self.track_label.config(text=f"Microphone: {self.master.mic.name}\nTrack Id: {idx}")
        self.audible_button.config(command=lambda: self._click(idx, True, True, num2date(lower_limit_start), num2date(upper_limit_start)),
                                   state=tk.NORMAL, )
        self.inaudible_button.config(command=lambda: self._click(idx, True, False), state=tk.NORMAL)
        self.unknown_button.config(command=lambda: self._click(idx, False, False), state=tk.NORMAL)
        self.progress_label.config(text=f"{self.i}/{self.master.tracks.track_id.nunique()}")

        canvas = FigureCanvasTkAgg(fig, master=self)  # A tk.DrawingArea.
        canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew', rowspan=3)

    def _click(self, id_: Any, valid: bool, audible: Optional[bool] = None,
               starttime: Optional[str] = None, endtime: Optional[bool] = None):
        """
        Save an annotation depending on what button what audibility button was clicked and clear
        the frame to be able to show the next plot.

        Parameters
        ----------
        id_ : Any
            The track unique identifier.
        valid : bool
            If the track was valid.
        audible : bool
            If the track was valid, was it audible.
        starttime : str, default None
            If the track was audible, when does audibility start.
        endtime : str, default None
            If the track was audible, when does audibility end.
        """
        self.audible_button.config(state=tk.DISABLED)
        self.inaudible_button.config(state=tk.DISABLED)
        self.unknown_button.config(state=tk.DISABLED)

        self.master.add_annotation(id_, valid, audible, starttime, endtime)

        # Close the plot to clear the window for the next one.
        plt.close()

        self._next()
