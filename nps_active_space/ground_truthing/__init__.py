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

from nps_active_space.utils import calculate_spline, Microphone, Nvspl, Tracks


_app = None


def launch(*args, **kwargs):
    """A wrapper function to launch the ground truthing application."""
    global _app
    _app = _App(*args, **kwargs)
    _app.mainloop()


class _AppFrame(tk.Frame):
    """Abstract base class for all application frames."""

    def __init__(self, master):
        super().__init__(master)
        self.master = master


class _App(ABC, tk.Tk):

    def __init__(self, outfile: str,
                 # tracks: Tracks, nvspl: Nvspl, mic: Microphone, crs,
                 # site_shp: Optional[gpd.GeoDataFrame] = None, tracks_datetime = None,
                 # clip: Optional[bool] = None
                 ):
        super().__init__()

        self.outfile = outfile
        # self.crs = crs
        # self.mic = mic.to_crs(crs)
        # self.site_shp = site_shp.to_crs(crs)
        # self.tracks = gpd.clip(tracks.to_crs(crs), self.site_shp) if clip else tracks.to_crs(crs)
        # self.nvspl = nvspl
        # self.tracks_datetime = tracks_datetime

        # Set app features.
        self.title('NPS Active Space: Ground Truthing Module')
        self.iconbitmap()
        self.geometry('1200x600')

        # Create app menu.
        self.menu = tk.Menu(self)
        self.file_menu = tk.Menu(self.menu, tearoff=False)
        self.file_menu.add_command(label='Save...', command=self._save)
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
        self.switch_frame(_GroundTruthingFrame)

    def switch_frame(self, frame_class: Type[_AppFrame]):
        """Switch the frame that is being displayed in the application window.

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

    def add(self, id_: Any, valid: bool, audible: bool, starttime: str, endtime: str):
        """Add a new annotation.

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

    def _close(self):
        """A function to safely close the application. If the user has unsaved changes, they will be warned and asked
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
        """Determine what app frame should be shown next depending on if the user would like saved annotations to
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
    """Frame to be shown after all tracks have been classified,

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
        continue_button = tk.Button(
            self,
            text='Continue >>',
            width=20,
            font=('Avenir', 8),
            bg='ivory2',
            command=lambda: self.master.switch_frame(_PlotResultsFrame)
        )

        # Place widgets.
        frame_label.place(relx=0.5, rely=0.45, anchor='center')
        continue_button.place(relx=0.9, rely=0.9, anchor='center')


class _PlotResultsFrame(_AppFrame):

    def __init__(self, master):
        super().__init__(master)


class _GroundTruthingFrame(_AppFrame):

    def __init__(self = None, master = None):
        super().__init__(master)
#         self.data = iter(self.master.tracks.groupby(self.master.track_id, 0, **('by', 'axis')))
#         self.canvas = None
#         self.i = 0
#         self.audible_button = tk.Button(self, 'Audible >>', 'green', 'white', 10, ('Avenir', 12, 'bold'), **('text', 'bg', 'fg', 'width', 'font'))
#         self.inaudible_button = tk.Button(self, 'Inaudible >>', 'red', 'white', 10, ('Avenir', 12, 'bold'), **('text', 'bg', 'fg', 'width', 'font'))
#         self.unknown_button = tk.Button(self, 'Unknown >>', 'yellow', 12, ('Avenir', 11, 'bold'), **('text', 'bg', 'width', 'font'))
#         self.progress_label = tk.Label(self, 'ivory2', **('bg',))
#         self.grid_columnconfigure(0, 5, **('weight',))
#         self.grid_columnconfigure(1, 1, **('weight',))
#         self.grid_rowconfigure(0, 2, **('weight',))
#         self.grid_rowconfigure(1, 1, **('weight',))
#         self.grid_rowconfigure(2, 2, **('weight',))
#         self.audible_button.grid(1, 1, 'n', **('row', 'column', 'sticky'))
#         self.inaudible_button.grid(1, 1, **('row', 'column'))
#         self.unknown_button.grid(1, 1, 's', **('row', 'column', 'sticky'))
#         self.progress_label.grid(0, 1, 'ne', 10, 5, **('row', 'column', 'sticky', 'padx', 'pady'))
#         self.next_graph()
#
#
#     def clicked(self, id_ = None, valid = None, audible = None, starttime = (False, None, None), endtime = {
#         'valid': bool,
#         'audible': bool,
#         'starttime': str,
#         'endtime': str }):
#         self.master.add(id_, valid, audible, starttime, endtime)
#         self.master._saved = False
#         plt.close()
#         if self.canvas:
#             self.canvas.get_tk_widget().destroy()
#         self.next_graph()
#
#
#     def next_graph(self):