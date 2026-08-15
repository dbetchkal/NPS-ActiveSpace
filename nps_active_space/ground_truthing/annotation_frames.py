import geopandas as gpd
import tkinter as tk
from tkinter import filedialog, messagebox
from matplotlib.dates import date2num
from shapely.geometry import Point

from nps_active_space.ground_truthing.frame_base import _AppFrame
from nps_active_space.ground_truthing.segments import AudibleRange
from nps_active_space.ground_truthing import annotation_plot
from nps_active_space.ground_truthing.annotation_plot import snap_threshold_days
from nps_active_space.ground_truthing import dem
from nps_active_space.ground_truthing import segments
from nps_active_space.ground_truthing import track_context
from nps_active_space.ground_truthing import vehicle_info
from nps_active_space.utils.enums import TrackSource


class _AnnotationLoadFrame(_AppFrame):
    """
    A frame to allow the user to decide if they would like previously saved annotations to be loaded.

    Parameters
    ----------
    master : tk.Tk
        The tkinter window this frame will be shown in.
    """
    def __init__(self, master: tk.Tk) -> None:
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

    def _clear_yes(self) -> None:
        """Remove the Select File and related widgets if No option is selected."""
        self.select_file_button.place_forget()
        self.select_file_label.place_forget()
        self.select_file_label.config(text='')
        self.annotation_filename.set('')

        self.create_file_button.place(relx=0.6, rely=0.53, anchor='w')

    def _clear_no(self) -> None:
        """Remove the Create File and related widgets if Yes option is selected."""
        self.create_file_button.place_forget()
        self.create_file_label.place_forget()
        self.create_file_label.config(text='')
        self.annotation_filename.set('')

        self.select_file_button.place(relx=0.6, rely=0.48, anchor='w')

    def _select_file(self) -> None:
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

    def _create_file(self) -> None:
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

    def _option_selected(self) -> None:
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
    def __init__(self, master: tk.Tk) -> None:
        super().__init__(master)

        content = tk.Frame(self, bg='ivory2')

        # Define widgets.
        frame_label = tk.Label(
            content,
            text='Instructions:',
            font=('Avenir', 14, 'bold'),
            bg='ivory2'
        )
        instructions = tk.Label(
            content,
            text=(
                'Annotate audible extent with the range sliders under the spectrogram.\n'
                'Use "Add Sound Event" to annotate another extent on the same track.\n'
                'If no periods are audible, click "Remove" to remove all extents.\n'
                'When finished annotating a track, click "Submit" to move to the next track.\n'
                'Click "Ignore" to dismiss a track from analysis for any reason.'
            ),
            font=('Avenir', 14),
            bg='ivory2',
            justify='center',
        )
        save_reminder = tk.Label(
            content,
            text='As always, make sure to save intermittently!',
            font=('Avenir', 14),
            bg='ivory2'
        )
        start_button = tk.Button(
            self,
            text='Start',
            font=('Avenir', 8),
            width=20,
            bg='ivory2',
            command=lambda: self.master.switch_frame(_GroundTruthingFrame)
        )
        back_button = tk.Button(
            self,
            text='<< Back',
            font=('Avenir', 8),
            width=20,
            bg='ivory2',
            command=lambda: self.master.switch_frame(_AnnotationLoadFrame)
        )

        frame_label.pack(pady=(0, 8))
        instructions.pack(pady=(0, 16))
        save_reminder.pack()

        # Place widgets.
        content.place(relx=0.5, rely=0.45, anchor='center')
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
    def __init__(self, master: tk.Tk) -> None:
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
    def __init__(self, master: tk.Tk) -> None:
        super().__init__(master)

        # Set frame variables to starting values.
        self.data = list(self.master.tracks.groupby(by='track_id'))
        self.i = 0

        self.dem_x, self.dem_y, self.dem_z = dem.process_dem(self.master.dem, self.master.tracks.crs)

        # load aircraft data, if applicable
        self.faa = vehicle_info.load_faa(
            self.master.faa_path,
            self.master.faa_corrections_path,
            self.master.track_source,
            self.master.tracks,
        )

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
        self.ignore_button = tk.Button(
            self,
            text='Ignore >>',
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
        self.next_annotated_button = tk.Button(
            self,
            text='Next Annotated >',
            bg='ivory2',
            fg='black',
            width=15,
            font=('Avenir', 12),
            command=self._next_annotated
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
        self.last_annotated_button = tk.Button(
            self,
            text='Last Annotated',
            bg='ivory2',
            fg='black',
            width=15,
            font=('Avenir', 12),
            command=self._to_last_annotated
        )

        self.fig_widget = None

        # Place widgets.
        self.grid_columnconfigure(0, weight=5)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=2)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=3)
        self.grid_rowconfigure(3, weight=3)
        self.grid_rowconfigure(4, weight=1)
        self.grid_rowconfigure(5, weight=1)
        self.grid_rowconfigure(6, weight=1)
        self.grid_rowconfigure(7, weight=1)
        self.progress_label.grid(row=0, column=1, sticky='ne', padx=10, pady=5)
        self.track_label.grid(row=0, column=1, pady=10)
        self.time_label.grid(row=1, column=1, pady=10)
        self.submit_button.grid(row=2, column=1, sticky='n')
        self.ignore_button.grid(row=2, column=1, sticky='s')
        
        self.nav_buttons.grid(row=3, column=1)
        self.back_button.pack(side=tk.LEFT, padx=10)
        self.next_button.pack(side=tk.LEFT, padx=10)

        self.next_annotated_button.grid(row=4, column=1, padx=10, pady=1)
        self.next_unannotated_button.grid(row=5, column=1, padx=10, pady=1)
        self.next_identifier_button.grid(row=6, column=1, padx=10, pady=1)
        self.last_annotated_button.grid(row=7, column=1, padx=10, pady=1)

        self._load_index(0)
    
    
    def _load_index(self, i: int) -> None:
        """Move to the track with index `i`."""

        # set index
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
            # self._store_annotation() advances to the next track
            return

        # load spectrogram
        spectro = track_context.load_spectrogram(self.master.nvspl, points)
        if spectro.empty:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {track_id} has no accompanying spectrogram. Skipping...",
                icon='warning'
            )
            self._store_annotation(track_id, points, valid=False, note='No SPL data')
            # self._store_annotation() advances to the next track
            return
                
        mic_point = Point(float(self.master.mic.x), 
                          float(self.master.mic.y), 
                          float(self.master.mic.z))
        try:
            spline = track_context.prepare_spline(points, mic_point)
        except ValueError as exc:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {track_id} could not be interpolated ({exc}). Skipping...",
                icon='warning'
            )
            self._store_annotation(track_id, points, valid=False, note='Spline failed')
            return

        # Determine the closest spline point to the mic.
        closest_point, closest_time = track_context.closest_approach(spline)

        # Calculate some datetime starting points.
        x_lims = date2num(spectro.index)  # convert the NVSPL's nice datetime axis to numbers
        lower_limit_start, upper_limit_start = track_context.limit_line_bounds(closest_time, x_lims)

        if upper_limit_start <= lower_limit_start:
            tk.messagebox.showwarning(
                title='Data Warning',
                message=f"Track {track_id} is a double back path causing the limit lines to cross. Skipping...",
                icon='warning'
            )
            self._store_annotation(track_id, spline, valid=False, note='Crossed limit lines.')
            # self._store_annotation() advances to the next track
            return
        
        # load audible ranges - check for previous annotations
        # note: audible_ranges is a list of [datetime, datetime], in the SPECTROGRAM's timeline
        # - i.e. uses the 'time_audible' field, not the 'point_dt' field
        annots = self.master.annotations[self.master.annotations["_id"] == track_id]
        audible_ranges = track_context.audible_ranges_from_annotations(
            annots, spline, lower_limit_start, upper_limit_start
        )
        if audible_ranges is None:
            return
        
        vehicle = vehicle_info.lookup_track_vehicle(
            self.faa, self.master.track_source, track_id, points
        )
        self.vessel_name = vehicle.vessel_name

        # update track-specific state, stored in member variables
        # these are stored as member variables because _build_plot() might be called again later,
        # and we want it to always have access to these without having to pass them around constantly
        # nice to do this all at once to better keep track of these variables and not have only some of them update if we returned early
        self.track_id = track_id
        self.track_annotated = not annots.empty
        self.valid = annots["valid"].all()  # if invalid should only be one row anyways
        self.points = points
        self.spectro = spectro
        self.audible_ranges = audible_ranges
        self.spline = spline
        self.spline_time_num = date2num(spline["time_audible"].to_numpy())
        self.snap_threshold_days = snap_threshold_days(
            spline["time_audible"].diff().median()
        )
        self.closest_point = closest_point
        self.closest_time = closest_time
        self.x_lims = x_lims
        self.lower_limit_start = lower_limit_start
        self.upper_limit_start = upper_limit_start
        if vehicle.faa_row is not None or self.master.track_source is TrackSource.AIS:
            self.aircraft_help_text = vehicle.help_text
        else:
            self.aircraft_help_text = None
        self.aircraft_type = vehicle.vehicle_type

        self._build_plot()

    def _next(self) -> None:
        if self.i + 1 < len(self.data):
            self._load_index(self.i + 1)
        else:
            self.master.switch_frame(_CompletionFrame)

    def _back(self) -> None:
        if self.i <= 0:
            self.master.switch_frame(_InstructionsFrame)
        else:
            self._load_index(self.i - 1)

    def _next_unannotated(self) -> None:
        """iterate self.i until we find a track that hasn't been annotated"""
        while (self.i+1 < len(self.data)):
            self.i += 1
            track_id = self.data[self.i][0]
            if str(track_id) not in self.master.annotations._id.values:
                self._load_index(self.i)
                return
        # if all are annotated, we're done!
        self.master.switch_frame(_CompletionFrame)
    
    def _next_annotated(self) -> None:
        i = self.i
        while (i+1 < len(self.data)):
            i += 1
            track_id = self.data[i][0]
            if str(track_id) in self.master.annotations._id.values:
                self._load_index(i)
                return
        # if none are annotated, don't go anywhere
        tk.messagebox.showinfo(
            title='Info',
            message=f"None annotated."
        )

    def _next_identifier(self) -> None:
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

    def _to_last_annotated(self) -> None:
        """Search from the end until we find an annotated record"""
        for i in range(len(self.data)-1, -1, -1):
            track_id = self.data[i][0]
            if str(track_id) in self.master.annotations._id.values:
                self._load_index(i)
                return
        # if didn't find it, then presumably none annotated, warn the user
        tk.messagebox.showinfo(
            title='Status',
            message=f"Failed to find last annotated track. This is probably because no tracks have been annotated."
        )


    def _store_annotation(
        self,
        track_id: str,
        points: gpd.GeoDataFrame,
        audible_ranges: list[AudibleRange] = [],
        valid: bool = True,
        note: str | None = None,
    ) -> None:
        """
        Save an annotation depending on what button what audibility button was clicked and clear
        the frame to be able to show the next plot.

        Parameters
        ----------
        track_id : str
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
        self.ignore_button.config(state=tk.DISABLED)

        gdf = segments.build_annotation_segments(
            track_id, points, audible_ranges=audible_ranges, valid=valid, note=note
        )
        self.master.set_annotation(track_id, gdf)
        self._next()

    def _build_plot(self) -> None:
        annotation_plot.build_plot(self)
