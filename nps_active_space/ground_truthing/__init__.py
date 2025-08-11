import datetime as dt
import tkinter as tk
import traceback
from abc import ABC
from tkinter import filedialog, messagebox
from typing import Any, List, Optional, Type, TYPE_CHECKING
from functools import partial

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import numpy as np
import pandas as pd
from matplotlib.dates import date2num, DateFormatter, num2date
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RangeSlider, Button
from PIL import Image, ImageTk
from shapely.geometry import LineString, Point
import rasterio
import pyproj

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.utils import Annotations, audible_time_delay, interpolate_spline, expected_relative_Lp, FAAReleasable

if TYPE_CHECKING:
    from nps_active_space.utils import Microphone, Nvspl, Tracks

# mplstyle.use("fast")

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
    database_type: str
        Database type. 'GPS', 'ADSB', or 'AIS'.
    dem
        rasterio Dataset object for reading the DEM
    clip : bool, default False
        If True, clip the Tracks to the study area.
    faa_path : str, default None
        If using aircraft data, path to the FAA Releasable database MASTER.txt file. Leave this as None if using AIS data.
    faa_corrections_path : str
        If using aircraft data, path to the FAA Releasable database corrections json file. Leave this as None if using AIS data.
    """
    def __init__(self, mic: 'Microphone', nvspl: 'Nvspl', tracks: 'Tracks',
                 crs: str, study_area: gpd.GeoDataFrame, database_type: str, dem, clip: bool = False,
                 faa_path: str = None, faa_corrections_path: str = ""):
        super().__init__()

        self.crs = crs
        self.mic = mic.to_crs(crs)
        self.study_area = study_area.to_crs(crs)
        self.tracks = gpd.clip(tracks.to_crs(crs), self.study_area) if clip else tracks.to_crs(crs)
        self.nvspl = nvspl
        self.outfile = None
        self.database_type = database_type
        self.dem = dem
        self.faa_path = faa_path
        self.faa_corrections_path = faa_corrections_path

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
        self.switch_frame(_GroundTruthingFrame)

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

    def set_annotation(self, track_id: str, annotated_lines: gpd.GeoDataFrame):
        """
        Add new audibility annotations for a track, replacing any previous annotations for that track.

        Parameters
        ----------
        track_id: str
            ID of the track being annotated. Convenient for searching out previous annotations for this track ID and removing them,
            since they will be overwritten by this new annotation.
        annotated_lines: gpd.GeoDataFrame
            a GeoDataFrame of annotated lines for a track to add to the overall annotations GeoDataFrame.
        """
        if annotated_lines.crs != self.annotations.crs:
            annotated_lines = annotated_lines.to_crs(self.annotations.crs)

        # remove old annotations for this track
        self.annotations = self.annotations[self.annotations["_id"] != track_id]

        # add new annotation
        self.annotations = pd.concat([self.annotations, annotated_lines], ignore_index=True).infer_objects()
        self._saved = False

        print("n segments saved", (self.annotations["_id"] == track_id).sum())

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
        if self.annotations.empty or self.annotations["valid"].sum() == 0:
            tk.messagebox.showinfo(
                title='Plot Tracks',
                message=f"No valid tracks to plot.",
            )
            return

        fig, ax = plt.subplots(1, 1, figsize=(6, 9))

        # Plot study area.
        study_area = self.study_area.to_crs('epsg:4326')
        study_area.geometry.boundary.plot(ax=ax, ls="--", color="navy", label='study area')

        # Plot track audibility.
        valid_segments = self.annotations[self.annotations.valid]
        if valid_segments.audible.any():
            valid_segments[valid_segments.audible == True].plot(
                ax=ax,
                color='deepskyblue',
                alpha=0.5,
                markersize=3,
                zorder=3,
                label="Audible segments"
            )
        if (~valid_segments.audible).any():
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
        self.load_annotations_bool = tk.BooleanVar(value=False)
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
            variable=self.load_annotations_bool,
            bg='ivory2',
            command=lambda: self._clear_no()
        )
        self.no_button = tk.Radiobutton(
            self,
            text='No, do not load prior annotations.',
            font=('Avenir', 10),
            value=False,
            variable=self.load_annotations_bool,
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

            if self.load_annotations_bool.get() is True:
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
        self.data = list(self.master.tracks.groupby(by='track_id'))
        self.i = -1  # will be incremented by 1 to start at 0, and go up

        self.dem_x, self.dem_y, self.dem_z = self.process_dem(self.master.tracks.crs)

        # load aircraft data, if applicable
        self.faa = None
        if self.master.faa_path is not None and self.master.faa_corrections_path is not None:
            if self.master.database_type == "ADSB":
                icaos = self.master.tracks["ICAO_address"].unique()
                self.faa = FAAReleasable(self.master.faa_path, self.master.faa_corrections_path, icao_addresses=icaos).data
                if self.faa.empty:
                    self.faa = None
            elif self.master.database_type == "GPS":
                track_ids = self.master.tracks["track_id"].unique()
                n_numbers = list(map(lambda s: s.split("_")[0][1:], track_ids))
                self.faa = FAAReleasable(self.master.faa_path, self.master.faa_corrections_path, n_numbers=n_numbers).data
                if self.faa.empty:
                    self.faa = None

        # Define widgets.
        self.progress_label = tk.Label(
            self,
            bg='ivory2'
        )
        self.track_label = tk.Label(
            self,
            bg='ivory2',
            font=('Avenir', 10, 'bold'),
            justify='left'
        )
        self.time_label = tk.Label(
            self,
            bg='ivory2',
            font=('Avenir', 10)
        )
        self.submit_button = tk.Button(
            self,
            text='Submit >>',
            bg='green',
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

        self.nav_buttons = tk.Frame(self)
        self.next_button = tk.Button(
            self.nav_buttons,
            text='Next >',
            bg='ivory2',
            fg='black',
            width=10,
            font=('Avenir', 12),
            command=self._next
        )
        self.back_button = tk.Button(
            self.nav_buttons,
            text='< Back',
            bg='ivory2',
            fg='black',
            width=10,
            font=('Avenir', 12),
            command=self._back
        )
        self.next_unannotated_button = tk.Button(
            self,
            text='Next Unannotated >',
            bg='ivory2',
            fg='black',
            width=15,
            font=('Avenir', 12),
            command=self._next_unannotated
        )
        self.next_identifier_button = tk.Button(
            self,
            text='Next Identifier >',
            bg='ivory2',
            fg='black',
            width=15,
            font=('Avenir', 12),
            command=self._next_identifier
        )

        # Place widgets.
        self.grid_columnconfigure(0, weight=5)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=2)
        self.grid_rowconfigure(3, weight=2)
        self.grid_rowconfigure(4, weight=1)
        self.grid_rowconfigure(5, weight=1)
        self.progress_label.grid(row=0, column=1, sticky='ne', padx=10, pady=5)
        self.track_label.grid(row=0, column=1, pady=10)
        self.time_label.grid(row=1, column=1, pady=10)
        self.submit_button.grid(row=2, column=1, sticky='n')
        self.unknown_button.grid(row=2, column=1, sticky='s')
        
        self.nav_buttons.grid(row=3, column=1)
        self.back_button.pack(side=tk.LEFT, padx=10)
        self.next_button.pack(side=tk.LEFT, padx=10)

        self.next_unannotated_button.grid(row=4, column=1, padx=10, pady=5)
        self.next_identifier_button.grid(row=5, column=1, padx=10, pady=5)

        self._next()
    
    
    def process_dem(self, crs):
        dem = self.master.dem
        data = dem.read(1)

        if dem.nodata is not None:
            data = np.ma.masked_equal(data, dem.nodata)
        
        # convert pixel coords to spatial coords
        x = np.arange(0, data.shape[1])
        y = np.arange(0, data.shape[0])
        x_coords, y_coords = np.meshgrid(x, y)  # pixel coords right now
        x_coords, y_coords = rasterio.transform.xy(dem.transform, y_coords, x_coords, offset="center") # now spatial coords
        x_coords = x_coords.reshape(data.shape)
        y_coords = y_coords.reshape(data.shape)

        # convert x,y to correct CRS
        transformer = pyproj.Transformer.from_crs(dem.crs, crs, always_xy=True)
        x_coords, y_coords = transformer.transform(x_coords, y_coords)

        return x_coords, y_coords, data

    
    def _load_index(self, i):
        """Move to the track with index `i`."""

        self.i = i

        # load points
        track_id, points = self.data[i]
        if points.shape[0] < 3:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {track_id} has fewer than 3 points and therefore cannot be processed. Skipping...",
                icon='warning'
            )
            self._store_annotation(track_id, points, valid=False, note='Too few points')
            # that called _next()
            return

        # load spectrogram
        time_pad = dt.timedelta(seconds=5*60)
        t_start = str(points.point_dt.iat[0] - time_pad)
        t_end = str(points.point_dt.iat[-1] + time_pad)
        spectro = self.master.nvspl.loc[t_start:t_end, '12.5':'20000']
        if spectro.empty:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {track_id} has no accompanying spectrogram. Skipping...",
                icon='warning'
            )
            self._store_annotation(track_id, points, valid=False, note='No SPL data')
            # that called _next()
            return
                
        points.sort_values(by='point_dt', ascending=True, inplace=True)
        spline = interpolate_spline(points)
        mic_point = Point(float(self.master.mic.x), 
                          float(self.master.mic.y), 
                          float(self.master.mic.z))
        spline = audible_time_delay(spline, 'point_dt', mic_point)
        spline = expected_relative_Lp(spline, mic_point)

        # Determine the closest spline point to the mic.
        closest_point = spline[spline.distance_to_target == spline.distance_to_target.min()]
        closest_time = spline.loc[spline.distance_to_target.idxmin()]['time_audible']

        # Calculate some datetime starting points.
        x_lims = date2num(spectro.index)  # convert the NVSPL's nice datetime axis to numbers
        lower_limit_start = max(date2num(closest_time - dt.timedelta(seconds=60)), x_lims[0])
        upper_limit_start = min(date2num(closest_time + dt.timedelta(seconds=60)), x_lims[-1])

        if upper_limit_start <= lower_limit_start:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {track_id} is a double back path causing the limit lines to cross. Skipping...",
                icon='warning'
            )
            self._store_annotation(track_id, spline, valid=False, note='Crossed limit lines.')
            # that called _next()
            return
        
        # load audible ranges - check for previous annotations
        # note: audible_ranges is a list of [datetime, datetime], in the SPECTROGRAM's timeline
        # - i.e. uses the 'time_audible' field, not the 'point_dt' field
        annots = self.master.annotations[self.master.annotations["_id"] == track_id]
        if annots.empty:
            # default range slider
            audible_ranges = [[num2date(lower_limit_start), num2date(upper_limit_start)]]
        else:
            audible_ranges = []
            for _, a in annots[annots["valid"] & annots["audible"]].iterrows():
                # note that invalid or fully inaudible annotations will result in no ranges, as desired
                time_audible_start = spline[spline["point_dt"] == a["start_dt"]].iloc[0]["time_audible"]
                time_audible_end = spline[spline["point_dt"] == a["end_dt"]].iloc[0]["time_audible"]
                audible_ranges.append([time_audible_start, time_audible_end])
        
        # load FAA data if applicable
        faa_row = None
        if self.faa is not None:
            if self.master.database_type == "ADSB":
                icao = points.iloc[0]["ICAO_address"]
                matching = self.faa[self.faa["MODE S CODE HEX"] == icao]
                if not matching.empty:
                    faa_row = matching.iloc[0]
                    aircraft_help_text = f"N-Number: {faa_row['N-NUMBER']}"
            elif self.master.database_type == "GPS":
                n_number = track_id.split("_")[0][1:]
                matching = self.faa[self.faa["N-NUMBER"] == n_number]
                if not matching.empty:
                    faa_row = matching.iloc[0]
                    aircraft_help_text = f"ICAO Address: {faa_row["MODE S CODE HEX"]}"

        # update track-specific state, stored in member variables
        # these are stored as member variables because _build_plot() might be called again later,
        # and we want it to always have access to these without having to pass them around constantly
        # nice to do this all at once to better keep track of these variables and not have only some of them update if we returned early
        self.track_id = track_id
        self.track_annotated = not annots.empty
        self.points = points
        self.spectro = spectro
        self.audible_ranges = audible_ranges
        self.spline = spline
        self.typical_t_diff = spline["time_audible"].diff().median()
        self.closest_point = closest_point
        self.closest_time = closest_time
        self.x_lims = x_lims
        self.lower_limit_start = lower_limit_start
        self.upper_limit_start = upper_limit_start
        self.aircraft_help_text = aircraft_help_text if faa_row is not None else None
        self.aircraft_type = faa_row["TYPE AIRCRAFT"] if faa_row is not None else None

        self._build_plot()

    def _next(self):
        if self.i + 1 < len(self.data):
            self._load_index(self.i + 1)
        else:
            self.master.switch_frame(_CompletionFrame)

    def _back(self):
        if self.i <= 0:
            self.master.switch_frame(_InstructionsFrame)
        else:
            self._load_index(self.i - 1)

    def _next_unannotated(self):
        """iterate self.i until we find a track that hasn't been annotated"""
        while (self.i+1 < len(self.data)):
            self.i += 1
            track_id = self.data[self.i][0]
            if str(track_id) not in self.master.annotations._id.values:
                self._load_index(self.i)
                return
        # if all are annotated, we're done!
        self.master.switch_frame(_CompletionFrame)

    def _next_identifier(self):
        """iterate self.i until we find a different vehicle identifier"""
        current_id = self.track_id.split("_")[0]
        while (self.i+1 < len(self.data)):
            self.i += 1
            track_id = self.data[self.i][0]
            if str(track_id).split("_")[0] != current_id:
                self._load_index(self.i)
                return
        # if all are annotated, we're done!
        self.master.switch_frame(_CompletionFrame)

    def _store_annotation(self, track_id: Any, points: gpd.GeoDataFrame, audible_ranges: Optional[list] = [],
                           valid: bool = True, note: Optional[str] = None):
        """
        Save an annotation depending on what button what audibility button was clicked and clear
        the frame to be able to show the next plot.

        Parameters
        ----------
        track_id : Any
            The track unique identifier.
        points: gpd.GeoDataFrame:
            Track and spline points to annotate.
        valid : bool, default True
            If the track was valid.
        audible_ranges: list of [datetime, datetime], default []
            Periods of time when the track was audible. If an empty list, everything was inaudible.        
        note: str, default None
            Any note to be added to all points passed for annotation.
        """
        # Deactivate the decision buttons.
        self.submit_button.config(state=tk.DISABLED)
        self.unknown_button.config(state=tk.DISABLED)

        # Convert points to WGS84 to avoid geopandas bug mentioned in Track model :(
        if 'z' not in points.columns:
            points['z'] = points.geometry.z
            points = points.to_crs('epsg:4326')
            points['geometry'] = points.apply(lambda row: Point(row.geometry.x, row.geometry.y, row.z), axis=1)
            points.drop('z', axis=1, inplace=True)
        

        segments = []

        # Saving invalid tracks or fully inaudible tracks can be done with one line
        if valid is False or len(audible_ranges) == 0:
            segments.append({
                '_id': track_id,
                'start_dt': points.point_dt.iat[0],
                'end_dt': points.point_dt.iat[-1],
                'valid': valid,
                'audible': False,
                'note': note,
                'geometry': points.geometry.iat[0] if points.shape[0] == 1
                                else LineString(points.geometry.tolist())
            })

        else:
            # note that in this conditional block, audible_ranges is not empty

            # simplify overlapping audible ranges
            audible_ranges = _collapse_audible_ranges(audible_ranges)

            # remove timezone info to avoid errors comparing tz-naive against tz-aware datetimes
            for r in audible_ranges:
                for i in range(len(r)):
                    r[i] = r[i].replace(tzinfo=None)

            segments = []
            
            # add segment for the inaudible range at the beginning, if it exists
            first_inaudible_segment = points[points.time_audible < audible_ranges[0][0]]
            if first_inaudible_segment.shape[0] >= 2:
                segments.append(
                    {'_id': track_id,
                    'start_dt': first_inaudible_segment.point_dt.iat[0],
                    'end_dt': first_inaudible_segment.point_dt.iat[-1],
                    'valid': True,
                    'audible': False,
                    'note': note,
                    'geometry': LineString(first_inaudible_segment.geometry.tolist())}
                )

            # add segments for each audible range and the inaudible range following it
            for i, r in enumerate(audible_ranges):                
                audible_segment = points[(points.time_audible >= r[0]) & (points.time_audible < r[1])]
                if audible_segment.shape[0] >= 2:
                    segments.append(
                        {'_id': track_id,
                        'start_dt': audible_segment.point_dt.iat[0],
                        'end_dt': audible_segment.point_dt.iat[-1],
                        'valid': True,
                        'audible': True,
                        'note': note,
                        'geometry': LineString(audible_segment.geometry.tolist())}
                    )

                if i+1 < len(audible_ranges):
                    next_start = audible_ranges[i+1][0]
                else:
                    next_start = points.time_audible.iat[-1]
                
                inaudible_segment = points[(points.time_audible >= r[1]) & (points.time_audible < next_start)]
                if inaudible_segment.shape[0] >= 2:
                    segments.append(
                        {'_id': track_id,
                        'start_dt': inaudible_segment.point_dt.iat[0],
                        'end_dt': inaudible_segment.point_dt.iat[-1],
                        'valid': True,
                        'audible': False,
                        'note': note,
                        'geometry': LineString(inaudible_segment.geometry.tolist())}
                    )

        gdf = gpd.GeoDataFrame(segments, geometry='geometry', crs=points.crs)
        self.master.set_annotation(track_id, gdf)
        plt.close("all")
        self._next()

    def _build_plot(self):
        """
        Build the matplotlib GridSpec plot for a track, using class state set by _load_index().
        """
        # ************************************ Build Plot ************************************#

        fig = plt.figure(figsize=(9, 5), constrained_layout=True)
        fig.canvas.manager.set_window_title(f"Microphone: {self.master.mic.name}, Track Id: {self.track_id}")

        height_ratios = [12, 6, 1] + [1 for _ in range(len(self.audible_ranges))] + [1]
        grid = GridSpec(ncols=2, nrows=4+len(self.audible_ranges), figure=fig,
                        width_ratios=[10,1], height_ratios=height_ratios)

        map_ax = fig.add_subplot(grid[0, 0])
        spectro_ax = fig.add_subplot(grid[1, 0])
        cue_ax = fig.add_subplot(grid[2, 0])
        slider_axes = []
        rm_range_button_axes = []
        for i in range(len(self.audible_ranges)):
            slider_axes.append(fig.add_subplot(grid[3+i, 0]))
            rm_range_button_axes.append(fig.add_subplot(grid[3+i, 1]))
        new_range_ax = fig.add_subplot(grid[grid.nrows-1, 0])

        # make some axes member variables too so that on_mouse_move() can see them
        self.fig = fig
        self.slider_axes = slider_axes
        self.spectro_ax = spectro_ax

        # --------------------------------- Plot Track --------------------------------- #

        # Display the study area, track points, spline points, closest point, and microphone
        map_ax.contour(
            self.dem_x,
            self.dem_y,
            self.dem_z,
            levels=8,
            colors="lightgray",
            linewidths=1,
            zorder=0
        )
        self.master.study_area.geometry.boundary.plot(
            label='study area',
            ax=map_ax,
            ls="--",
            lw=0.5,
            color="blue"
        )
        self.spline.plot(
            label='interpolated spline point',
            ax=map_ax,
            color="grey",
            zorder=1,
            markersize=0.1,
            alpha=0.5
        )
        self.points.plot(
            label='track point',
            ax=map_ax,
            color="blue",
            zorder=1,
            markersize=2,
        )
        self.closest_point.plot(
            label='closest point',
            ax=map_ax,
            color="lime",
            zorder=1,
            markersize=32,
            edgecolor="blue"
        )

        # Point indicating where on the map corresponds to a certain spot in the spectrogram,
        # when the user hovers their mouse. Init at the microphone, but invisible
        self.map_mousehover_point = map_ax.plot(
            self.master.mic.x,
            self.master.mic.y,
            ms=10,
            marker="o",
            markerfacecolor="none",
            markeredgecolor="black",
            zorder=3
        )[0]
        self.map_mousehover_point.set_visible(False)

        map_ax.plot(
            self.master.mic.x,
            self.master.mic.y,
            label='microphone',
            ls="",
            marker="x",
            ms=7,
            color="magenta",
            zorder=10
        )

        # Glean the spatial extent of the points. This will result in a square map.
        xmin, ymin, xmax, ymax = self.master.study_area.total_bounds
        map_ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        map_ax.set_aspect((xmax - xmin) / (ymax - ymin))
        map_ax.tick_params(axis='both', labelsize=6)
        map_ax.ticklabel_format(style='plain')  # disable scientific notation
        map_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # --------------------------------- Plot Spectrogram --------------------------------- #

        spectro_ax.imshow(
            self.spectro.T,
            origin="lower",
            aspect=(self.x_lims[-1] - self.x_lims[0]) / (8 * 33),
            cmap="plasma",
            extent=[self.x_lims[0], self.x_lims[-1], 0, self.spectro.shape[1]],
            interpolation=None,
            vmin=-10,
            vmax=80
        )

        # add in the time of closest approach in red.
        spectro_ax.axvline(
            date2num(self.closest_time),
            alpha=0.7,
            color="lime",
            zorder=2,
            linewidth=3,
            label='closest track point'
        )

        spectro_ax.legend(bbox_to_anchor=(0.25, 1.4))
        spectro_ax.set_yticks(np.arange(self.spectro.shape[1])[::6])
        spectro_ax.set_yticklabels(self.spectro.columns.astype('float')[::6])
        spectro_ax.set_ylabel("Freq. (Hz)", labelpad=15)
        spectro_ax.axhline(8, lw=1.0, color="white", ls=":", alpha=0.4, zorder=200)
        spectro_ax.xaxis_date()  # tell matplotlib that the numeric axis should be formatted as dates
        spectro_ax.xaxis.set_major_formatter(DateFormatter("%b-%d\n%H:%M"))  # tidy them!

        # --------------------------------- Plot Expected Lp Cue Signal --------------------------------- #

        ys = np.zeros(len(self.spline["time_audible"]))
        cue_ax.scatter(self.spline["time_audible"], ys, c=self.spline["Lp_est"], cmap="plasma")
        cue_ax.set_xlim(self.x_lims[0], self.x_lims[-1])
        cue_ax.set_ylabel("Sound Level Cue", rotation=0, ha='right', va='center')
        # remove extra plot bits
        for spine in cue_ax.spines.values():
            spine.set_visible(False)
        cue_ax.set_xticks([])
        cue_ax.set_yticks([])

        # --------------------------------- UI for Selecting Audible Ranges --------------------------------- #

        # note - we need to save references to the AudibleRangeUI objects, so the garbage collector doesn't trash them
        self.audible_range_uis = []
        self.rm_range_buttons = []  # same for storing these references

        for i, r in enumerate(self.audible_ranges):

            label = f"Audible Extent {i + 1}"
            ui = AudibleRangeUI(label, r, fig, map_ax, spectro_ax, slider_axes[i], self.spline, self.x_lims)
            self.audible_range_uis.append(ui)

            rm_button = Button(rm_range_button_axes[i], "Remove")
            rm_button.on_clicked(partial(self.remove_audible_range, i))  # use functools.partial to bind the current value of i, since i will change in the for loop
            self.rm_range_buttons.append(rm_button)
            

        # --------------------------------- New Slider Button --------------------------------- #

        self.new_slider_button = Button(new_range_ax, "Add Sound Event")
        self.new_slider_button.on_clicked(self.new_audible_range)

        # --------------------------------- Update Track Labels and Event Handling --------------------------------- #

        self.track_label.config(text=f"Microphone: {self.master.mic.name}\n\n" + \
                                     f"Track Id: {self.track_id}\n" + \
                                     (f"{self.aircraft_help_text}\n" if self.aircraft_help_text is not None else "") + \
                                     (f"Aircraft Type: {self.aircraft_type}\n" if self.aircraft_type is not None else "") + \
                                     f"\nAnnotated: {self.track_annotated}")
        self.progress_label.config(text=f"{self.i+1}/{self.master.tracks.track_id.nunique()}")
        self.submit_button.config(command=lambda: self._store_annotation(self.track_id, self.spline, self.audible_ranges), state=tk.NORMAL)
        self.unknown_button.config(command=lambda: self._store_annotation(self.track_id, self.spline, valid=False), state=tk.NORMAL)
        fig.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)

        # --------------------------------- Show Plot --------------------------------- #

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew', rowspan=100)  # large rowspan so we don't have to update it if adding new rows

    
    def new_audible_range(self, _):
        """
        Add a new audible range. Since we can't change the layout of axes after making them,
        we need to clear the current plot and remake it, taking into account the new audible range.
        """
        plt.close("all")
        self.audible_ranges.append([
            num2date(self.lower_limit_start),
            num2date(self.upper_limit_start)
        ])
        self._build_plot()
    
    def remove_audible_range(self, i, _):
        """Remove a certain audible range. See new_audible_range() for why we have to replot the figure."""
        plt.close("all")
        del self.audible_ranges[i]
        self._build_plot()
    
    def on_mouse_move(self, event):
        if event.inaxes == self.spectro_ax or event.inaxes in self.slider_axes:
            dt = num2date(event.xdata).replace(tzinfo=None)

            # update time display
            self.time_label.config(text=f"Cursor Time: {dt.strftime('%H:%M:%S')}")

            # get closest spline point to the mouse position
            closest_idx = (self.spline["time_audible"] - dt).abs().idxmin()
            closest_pt = self.spline.loc[closest_idx]

            # if the mouse is close enough to a point, display the marker on the map
            if abs(closest_pt["time_audible"] - dt) > self.typical_t_diff:
                self.map_mousehover_point.set_visible(False)
            else:
                self.map_mousehover_point.set_data(
                    [closest_pt.geometry.x], [closest_pt.geometry.y])
                self.map_mousehover_point.set_visible(True)

            self.fig.canvas.draw_idle()
            
    

class AudibleRangeUI():
    """Class to manage the various UI components of an audible range"""
    def __init__(self, label, range_bounds, fig, map_ax, spectro_ax, slider_ax, spline, x_lims):
        self.range_bounds = range_bounds  # where we store the range limits, part of the _GroundTruthingFrame.audible_ranges list
        self.fig = fig
        self.map_ax = map_ax
        self.spectro_ax = spectro_ax
        self.slider_ax = slider_ax
        self.spline = spline
        self.x_lims = x_lims

        low_init = date2num(range_bounds[0])
        high_init = date2num(range_bounds[1])

        self.highlight = map_ax.plot(
            spline.geometry.x,
            spline.geometry.y,
            lw=8,
            color='deepskyblue',
            ls='-',
            zorder=1,
            alpha=0.4
        )[0]
        self.lower_limit_line = spectro_ax.axvline(
            low_init,
            ls="--",
            alpha=0.7,
            color="white",
            zorder=2,
            linewidth=1,
        )
        self.upper_limit_line = spectro_ax.axvline(
            high_init,
            ls="--",
            alpha=0.7,
            color="white",
            zorder=2,
            linewidth=1
        )

        self.slider = RangeSlider(
            slider_ax,
            label=label,
            valmin=x_lims[0],
            valmax=x_lims[-1],
            valinit=[low_init, high_init]
        )
        self.slider.valtext.set_visible(False)  # Turn off range slider value label.
        self.slider.on_changed(self._slider_update)
        self._slider_update([low_init, high_init])

    def _slider_update(self, val: List):
        """
        Update spline highlight and spectrogram lines based on slider values.

        Parameters
        ----------
        val : List
            A two item list with the [min, max] values of the range slider.
        """
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
            self.spline.time_audible >= num2date(lower_t).replace(tzinfo=None),
            self.spline.time_audible <= num2date(upper_t).replace(tzinfo=None)
        ], axis=0)
        subset = self.spline.loc[subset_mask]
        self.highlight.set_data(subset.geometry.x, subset.geometry.y)


def _collapse_audible_ranges(ranges: list):
    """Collapse overlapping audible ranges into single audible ranges.
    
    Parameters
    ----------
    ranges: list of [datetime, datetime]
    
    Returns
    -------
    collapsed_ranges: list of [datetime, datetime]
    """
    
    times = []
    for start, end in ranges:
        times.append({"t": start, "type": "start"})
        times.append({"t": end, "type": "end"})
    times.sort(key=lambda x: x["t"])
    
    intervals = []
    n = 0
    start_time = None
    for t in times:
        if t["type"] == "start":
            # if it was previously quiet, start the next interval
            if n == 0:
                start_time = t["t"]
            n += 1
        elif t["type"] == "end":
            n -= 1
            if n == 0:
                # last sound ended, save the interval
                intervals.append([start_time, t["t"]])

    return intervals