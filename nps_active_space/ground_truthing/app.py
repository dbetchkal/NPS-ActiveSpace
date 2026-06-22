import traceback
from typing import TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tkinter as tk
from rasterio.io import DatasetReader
from tkinter import messagebox

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.models import Annotations
from nps_active_space.ground_truthing.frame_base import _AppFrame
from nps_active_space.ground_truthing.welcome_frame import _WelcomeFrame

if TYPE_CHECKING:
    from nps_active_space.utils.models import Microphone, Nvspl, Tracks

_app: "_App | None" = None


def launch(
    mic: "Microphone",
    nvspl: "Nvspl",
    tracks: "Tracks",
    crs: str,
    study_area: gpd.GeoDataFrame,
    database_type: TrackSource,
    dem: DatasetReader,
    clip: bool = False,
    faa_path: str | None = None,
    faa_corrections_path: str = "",
) -> None:
    """A wrapper function to launch the ground truthing application."""
    global _app
    _app = _App(
        mic,
        nvspl,
        tracks,
        crs,
        study_area,
        database_type,
        dem,
        clip=clip,
        faa_path=faa_path,
        faa_corrections_path=faa_corrections_path,
    )
    _app.mainloop()


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
    database_type: TrackSource
        Track feed for causal data (GPS, ADSB, or AIS).
    dem
        rasterio Dataset object for reading the DEM
    clip : bool, default False
        If True, clip the Tracks to the study area.
    faa_path : str, default None
        If using aircraft data, path to the FAA Releasable database MASTER.txt file. Leave this as None if using AIS data.
    faa_corrections_path : str
        If using aircraft data, path to the FAA Releasable database corrections json file. Leave this as None if using AIS data.
    """
    def __init__(
        self,
        mic: "Microphone",
        nvspl: "Nvspl",
        tracks: "Tracks",
        crs: str,
        study_area: gpd.GeoDataFrame,
        database_type: TrackSource,
        dem: DatasetReader,
        clip: bool = False,
        faa_path: str | None = None,
        faa_corrections_path: str = "",
    ) -> None:
        super().__init__()

        self.crs = crs
        self.mic = mic.to_crs(crs)
        self.study_area = study_area.to_crs(crs)
        self.nvspl = nvspl
        self.outfile = None
        self.database_type = database_type
        self.dem = dem
        self.faa_path = faa_path
        self.faa_corrections_path = faa_corrections_path

        # clip and project tracks, and then add z coordinate into geometry so it gets saved to GEOJSON properly later
        self.tracks = gpd.clip(tracks.to_crs(crs), self.study_area) if clip else tracks.to_crs(crs)
        self.tracks.geometry = gpd.points_from_xy(
            self.tracks.geometry.x, self.tracks.geometry.y, self.tracks.z)

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

        # Handle window closing
        self.protocol("WM_DELETE_WINDOW", self._close)

        self.switch_frame(_WelcomeFrame)


    def switch_frame(self, frame_class: type[_AppFrame]) -> None:
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

    def set_annotation(self, track_id: str, annotated_lines: gpd.GeoDataFrame) -> None:
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

    def load_annotations(self, filename: str) -> None:
        """
        Simple function to load existing annotations from a geojson file.

        Parameters
        ----------
        filename : str
            Absolute path to the geojson file to load previous annotations from.
        """
        self.annotations = Annotations(filename)

    def _close(self) -> None:
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
            plt.close("all")
            self.destroy()

    def _save(self) -> None:
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

    def _plot(self) -> None:
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
