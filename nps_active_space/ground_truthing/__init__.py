import datetime as dt
import tkinter as tk
import traceback
from abc import ABC
from tkinter import filedialog, messagebox
from typing import Any, List, Optional, Type, TYPE_CHECKING

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import date2num, DateFormatter, num2date
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RangeSlider
from PIL import Image, ImageTk
from shapely.geometry import LineString, Point

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.utils import Annotations, audible_time_delay, interpolate_spline

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
    mic : Microphone
        A Microphone object of the microphone deployment to be used for ground truthing.
    nvspl : Nvspl
        An Nvspl object of sound data record at the input microphone locations.
    tracks : Tracks
        a Tracks object of points to classify as audible, inaudible, or unknown from the microphone location.
    crs : str
        The PROJECTED coordinate system to be used for the Tracks, study area, and microphone.
        Format of 'epsg:XXXX...', E.g. 'epsg:32632'
    study_area : gpd.GeoDataFrame
        A gpd.GeoDataFrame of polygon(s) that make up the study area.
    clip : bool, default False
        If True, clip the Tracks to the study area.
    """
    def __init__(self, mic: 'Microphone', nvspl: 'Nvspl', tracks: 'Tracks',
                 crs: str, study_area: gpd.GeoDataFrame, clip: bool = False):
        super().__init__()

        self.crs = crs
        self.mic = mic.to_crs(crs)
        self.study_area = study_area.to_crs(crs)
        self.tracks = gpd.clip(tracks.to_crs(crs), self.study_area) if clip else tracks.to_crs(crs)
        self.nvspl = nvspl
        self.outfile = None

        # Set app features.
        self.title('NPS Active Space: Ground Truthing Module')
        self.iconbitmap(f"{ACTIVE_SPACE_DIR}/img/flat-four-color.ico")
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
        self.annotations = Annotations()
        self._saved = True
        self._frame = None

        self.switch_frame(_WelcomeFrame)

    def run(self):
        """Run the main application frame."""
        self.protocol("WM_DELETE_WINDOW", self._close)

        if set(self.tracks.track_id.unique()) - set(self.annotations._id.unique()) == set():
            self.switch_frame(_CompletionFrame)

        else:
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

    def add_annotation(self, annotated_lines: gpd.GeoDataFrame):
        """
        Add new track audibility annotations.

        Parameters
        ----------
        annotated_lines: gpd.GeoDataFrame
            a GeoDataFrame of annotated lines for a track to add to the overall annotations GeoDataFrame.
        """
        if annotated_lines.crs != self.annotations.crs:
            annotated_lines = annotated_lines.to_crs(self.annotations.crs)
        self.annotations = pd.concat([self.annotations, annotated_lines], ignore_index=True)
        self._saved = False

    def load_annotations(self, filename: str):
        """
        Simple function to load existing annotations from a geojson file.

        Parameters
        ----------
        filename : str
            Absolute path to the geojson file to load previous annotations from.
        """
        self.annotations = Annotations(filename)

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
        if self._saved is True:
            return
        try:
            self.annotations.to_file(self.outfile, driver='GeoJSON', mode='w', index=False)
            self._saved = True
            tk.messagebox.showinfo(
                title='Save Status',
                message=f"Saved!",
            )
        except Exception:
            tk.messagebox.showerror(
                title='Save Status',
                message=f"Unable to save.\n\n{traceback.format_exc()}",
            )

    def _plot(self):
        """Plot all annotated tracks and points."""

        fig, ax = plt.subplots(1, 1, figsize=(6, 9))

        # Plot study area.
        study_area = self.study_area.to_crs('epsg:4326')
        study_area.geometry.boundary.plot(ax=ax, ls="--", color="navy", label='study area')

        # Plot track audibility.
        valid_segments = self.annotations[self.annotations.valid]
        valid_segments[valid_segments.audible == True].plot(
            ax=ax,
            color='deepskyblue',
            alpha=0.5,
            markersize=3,
            zorder=3,
            label="Audible segments"
        )
        valid_segments[valid_segments.audible == False].plot(
            ax=ax,
            color='red',
            alpha=0.5,
            markersize=3,
            zorder=2,
            label="Inaudible segments"
        )

        # Plot microphone position.
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
        #cx.add_basemap(ax, crs='epsg:4326', source=cx.providers.OpenStreetMap.Mapnik) # TODO

        # This will result in a square map
        xmin, ymin, xmax, ymax = study_area.total_bounds
        pad = np.array([(xmax - xmin) * 0.1, (ymax - ymin) * 0.1])
        ax.set(xlim=(xmin - pad[0], xmax + pad[0]), ylim=(ymin - pad[1], ymax + pad[1]))
        ax.set_title("Annotated Track Segments")
        ax.tick_params(axis='both', labelsize=6)
        plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.35), markerscale=2)

        fig.show()


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
        im = Image.open(f"{ACTIVE_SPACE_DIR}/img/flat-four-color.png").resize((138, 181))
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

        self.create_file_button = tk.Button(
            self,
            text='Create File',
            bg='ivory2',
            command=lambda: self._create_file()
        )
        self.create_file_label = tk.Label(
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
            command=lambda: self._clear_no()
        )
        self.no_button = tk.Radiobutton(
            self,
            text='No, do not load prior annotations.',
            font=('Avenir', 10),
            value=False,
            variable=self.load_annotations,
            bg='ivory2',
            command=lambda: self._clear_yes()
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
        self.create_file_button.place(relx=0.6, rely=0.53, anchor='w')
        self.continue_button.place(relx=0.9, rely=0.9, anchor='center')

    def _clear_yes(self):
        """Remove the Select File and related widgets if No option is selected."""
        self.select_file_button.place_forget()
        self.select_file_label.place_forget()
        self.select_file_label.config(text='')
        self.annotation_filename.set('')

        self.create_file_button.place(relx=0.6, rely=0.53, anchor='w')

    def _clear_no(self):
        """Remove the Create File and related widgets if Yes option is selected."""
        self.create_file_button.place_forget()
        self.create_file_label.place_forget()
        self.create_file_label.config(text='')
        self.annotation_filename.set('')

        self.select_file_button.place(relx=0.6, rely=0.48, anchor='w')

    def _select_file(self):
        """Open File Dialog and save the existing selected annotation file."""
        filetypes = (('geojson files', '*.geojson'),)
        filename = filedialog.askopenfilename(
            title='Open file',
            initialdir='/',
            filetypes=filetypes
        )
        if filename:
            self.annotation_filename.set(filename)
            self.select_file_label.config(text=f"...{filename[-50:]}")
            self.select_file_label.place(relx=0.66, rely=0.48, anchor='w')

    def _create_file(self):
        """Open File Dialog and save the new annotation file."""
        filetypes = (('geojson files', '*.geojson'),)
        filename = filedialog.asksaveasfilename(
            title='Create annotation file',
            filetypes=filetypes,
            initialdir='/',
            initialfile=f"{self.master.mic.name}_saved_annotations",
            defaultextension=".geojson",
        )

        if filename:
            self.annotation_filename.set(filename)
            self.create_file_label.config(text=f"...{filename[-50:]}")
            self.create_file_label.place(relx=0.66, rely=0.53, anchor='w')

    def _option_selected(self):
        """If user wants to load existing annotations, load them before proceeding to the app instructions frame."""
        if self.annotation_filename.get():

            if self.load_annotations.get() is True:
                self.master.load_annotations(self.annotation_filename.get())

            self.master.outfile = self.annotation_filename.get()
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

            time_pad = dt.timedelta(seconds=5*60)
            spectro = self.master.nvspl.loc[str(points.point_dt.iat[0] - time_pad):str(points.point_dt.iat[-1] + time_pad), '12.5':'20000']

            # If the track is already annotated, move on.
            if str(idx) in self.master.annotations._id.values:
                self._next()

            # If the track does not have enough points for processing, mark it as invalid and move on.
            elif points.shape[0] < 3:
                tk.messagebox.showwarning(
                    title='Data Warning',
                    message=f"Track {idx} has fewer than 3 points and therefore cannot be processed. Skipping...",
                    icon='warning'
                )
                self._click(idx, points, valid=False, audible=False, note='Too few points')

            # If there is no spectrogram data for the track, mark it as invalid and move on.
            elif spectro.empty:
                tk.messagebox.showwarning(
                    title='Data Warning',
                    message=f"Track {idx} has no accompanying spectrogram. Skipping...",
                    icon='warning'
                )
                self._click(idx, points, valid=False, audible=False, note='No SPL data')

            # If this is an un-annotated track with 3+ points and corresponding SPL data, display its plot.
            else:
                self._build_plot(idx, points, spectro)

        except StopIteration:
            self.master.switch_frame(_CompletionFrame)

    def _click(self, id_: Any, points: gpd.GeoDataFrame, valid: bool, audible: bool,
               audibility_start: Optional[dt.datetime] = None, audibility_end: Optional[dt.datetime] = None,
               note: Optional[str] = None):
        """
        Save an annotation depending on what button what audibility button was clicked and clear
        the frame to be able to show the next plot.

        Parameters
        ----------
        id_ : Any
            The track unique identifier.
        points: gpd.GeoDataFrame:
            Track and spline points to annotate.
        valid : bool
            If the track was valid.
        audible : bool
            If the track was valid, was it audible.
        audibility_start : dt.datetime, default None
            If the track was audible, when does audibility start.
        audibility_end : dt.datetime, default None
            If the track was audible, when does audibility end.
        note: str, default None
            Any note to be added to all points passed for annotation.
        """
        # Deactivate the decision buttons.
        self.audible_button.config(state=tk.DISABLED)
        self.inaudible_button.config(state=tk.DISABLED)
        self.unknown_button.config(state=tk.DISABLED)

        # Convert points to WGS84 to avoid geopandas bug mentioned in Track model :(
        if 'z' not in points.columns:
            points['z'] = points.geometry.z
            points = points.to_crs('epsg:4326')
            points['geometry'] = points.apply(lambda row: Point(row.geometry.x, row.geometry.y, row.z), axis=1)
            points.drop('z', axis=1, inplace=True)

        # Unknown and inaudible tracks can be saved as a single line.
        if valid is False or audible is False:
            lines = gpd.GeoDataFrame(
                {
                    '_id': [id_],
                    'start_dt': [points.point_dt.iat[0]],
                    'end_dt': [points.point_dt.iat[-1]],
                    'valid': [valid],
                    'audible': [audible],
                    'note': [note],
                    'geometry': [points.geometry.iat[0] if points.shape[0] == 1
                                 else LineString(points.geometry.tolist())]
                },
                geometry='geometry',
                crs=points.crs
            )

        else:

            audible_segment = points[(points.time_audible >= audibility_start) &
                                     (points.time_audible <= audibility_end)]

            inaudible_segment_1 = points.loc[points.point_dt < audible_segment.point_dt.iat[0]]
            inaudible_segment_2 = points.loc[points.point_dt > audible_segment.point_dt.iat[-1]]

            line_segments = []
            for i, segment in enumerate([inaudible_segment_1, audible_segment, inaudible_segment_2]):
                if segment.shape[0] > 1:
                    line_segments.append(
                        {'_id': id_,
                         'start_dt': segment.point_dt.iat[0],
                         'end_dt': segment.point_dt.iat[-1],
                         'valid': True,
                         'audible': True if i == 1 else False,
                         'note': note,
                         'geometry': LineString(segment.geometry.tolist())}
                    )
            lines = gpd.GeoDataFrame(line_segments, geometry='geometry', crs=points.crs)

        self.master.add_annotation(lines)
        plt.close()
        self._next()

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
            subset = spline.loc[np.all(
                [spline.time_audible >= num2date(lower_t).replace(tzinfo=None),
                 spline.time_audible <= num2date(upper_t).replace(tzinfo=None)],
                axis=0)]
            highlight.set_data(subset.geometry.x, subset.geometry.y)

            self.audible_button.config(
                state=tk.NORMAL,
                command=lambda: self._click(
                    idx,
                    spline,
                    valid=True,
                    audible=True,
                    audibility_start=num2date(lower_t).replace(tzinfo=None),
                    audibility_end=num2date(upper_t).replace(tzinfo=None))
            )

            # Redraw the figure to ensure it updates
            fig.canvas.draw_idle()

        # Interpolate a track spline.
        points.sort_values(by='point_dt', ascending=True, inplace=True)
        spline = interpolate_spline(points)
        spline = audible_time_delay(spline, 'point_dt', Point(float(self.master.mic.x), 
                                                              float(self.master.mic.y), 
                                                              float(self.master.mic.z)))

        # Determine the closest spline point to the mic.
        closest_point = spline[spline.distance_to_target == spline.distance_to_target.min()]
        closest_time = spline.loc[spline.distance_to_target.idxmin()]['time_audible']

        # Calculate some datetime starting points.
        x_lims = date2num(spectrogram.index)  # convert the NVSPL's nice datetime axis to numbers
        lower_limit_start = max(date2num(closest_time - dt.timedelta(seconds=60)), x_lims[0])
        upper_limit_start = min(date2num(closest_time + dt.timedelta(seconds=60)), x_lims[-1])

        if upper_limit_start <= lower_limit_start:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {idx} is a double back path causing the limit lines to cross. Skipping...",
                icon='warning'
            )
            self._click(idx, spline, valid=False, audible=False, note='Crossed limit lines.')

        else:
            # ************************************ Build Plot ************************************#

            fig = plt.figure(figsize=(9, 5), constrained_layout=True)
            fig.canvas.manager.set_window_title(f"Microphone: {self.master.mic.name}, Track Id: {idx}")
            spec = GridSpec(ncols=1, nrows=10, figure=fig)
            ax1 = fig.add_subplot(spec[0:6, 0])
            ax2 = fig.add_subplot(spec[6:9, 0])
            ax3 = fig.add_subplot(spec[9, 0])

            # --------------------------------- Plot Track --------------------------------- #

            # Display the study area, track points, spline points, closest point, and microphone
            self.master.study_area.geometry.boundary.plot(
                label='study area',
                ax=ax1,
                ls="--",
                lw=0.5,
                color="blue"
            )
            spline.plot(
                label='interpolated spline point',
                ax=ax1,
                color="grey",
                zorder=1,
                markersize=0.5,
                alpha=0.1
            )
            points.plot(
                label='track point',
                ax=ax1,
                color="blue",
                zorder=1,
                markersize=2,
            )
            closest_point.plot(
                label='closest point',
                ax=ax1,
                color="red",
                zorder=1,
                markersize=3,
            )
            ax1.plot(
                self.master.mic.x,
                self.master.mic.y,
                label='microphone',
                ls="",
                marker="x",
                ms=7,
                color="magenta",
                zorder=10
            )

            highlight, = ax1.plot(
                spline.geometry.x,
                spline.geometry.y,
                lw=8,
                color='deepskyblue',
                ls='-',
                zorder=1,
                alpha=0.4
            )

            # Glean the spatial extent of the points. This will result in a square map.
            xmin, ymin, xmax, ymax = self.master.study_area.total_bounds
            ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
            ax1.set_aspect((xmax - xmin) / (ymax - ymin))
            ax1.tick_params(axis='both', labelsize=6)
            ax1.ticklabel_format(style='plain')  # disable scientific notation
            ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

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
            ax2.axhline(8, lw=1.0, color="white", ls=":", alpha=0.4, zorder=200)
            ax2.xaxis_date()  # tell matplotlib that the numeric axis should be formatted as dates
            ax2.xaxis.set_major_formatter(DateFormatter("%b-%d\n%H:%M"))  # tidy them!

            # --------------------------------- Plot Slider --------------------------------- #

            self.slider = RangeSlider(
                ax3,
                label="Audible Window",
                valmin=x_lims[0],
                valmax=x_lims[-1],
                valinit=[lower_limit_start, upper_limit_start]
            )

            self.slider.valtext.set_visible(False)  # Turn off range slider value label.
            self.slider.on_changed(_slider_update)
            _slider_update([lower_limit_start, upper_limit_start])

            # --------------------------------- Show Plot --------------------------------- #

            self.track_label.config(text=f"Microphone: {self.master.mic.name}\nTrack Id: {idx}")
            self.progress_label.config(text=f"{self.i}/{self.master.tracks.track_id.nunique()}")
            self.inaudible_button.config(command=lambda: self._click(idx, spline, valid=True, audible=False), state=tk.NORMAL)
            self.unknown_button.config(command=lambda: self._click(idx, spline, valid=False, audible=False), state=tk.NORMAL)

            canvas = FigureCanvasTkAgg(fig, master=self)
            canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew', rowspan=3)
