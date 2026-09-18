from __future__ import annotations

import configparser
import glob
import os
from collections.abc import Callable

from nps_active_space.ground_truthing.load_tracks import load_tracks
from nps_active_space.scripts.run_audible_transits import AudibleTransits
from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.helpers import get_deployment
from nps_active_space.viz.geometry import create_polyline_3d, track_points_to_linestring
from nps_active_space.viz.polyline_builders import TrackPolylineBuilder
from nps_active_space.viz.scene import ScenePlotter
from nps_active_space.viz.session import VizSession
from nps_active_space.viz.widgets import PlotWidgets


class TracksPlotter:
    def __init__(
        self,
        session: VizSession,
        widgets: PlotWidgets,
        scene: ScenePlotter,
        polylines: TrackPolylineBuilder,
    ) -> None:
        self._session = session
        self._widgets = widgets
        self._scene = scene
        self._polylines = polylines

    def plot_audible_transits(self, audible_transits_pkl: str | None = None) -> None:
        session = self._session
        style = session.style
        if audible_transits_pkl:
            session.status(f"Loading audible transits from {audible_transits_pkl}")
            listener = AudibleTransits.from_pickle(audible_transits_pkl)
        else:
            session.status("Loading audible transits")
            matches = glob.glob(
                os.path.join(
                    session.project_dir,
                    session.unit + session.site,
                    "Output_Data",
                    "AUDIBLE_TRANSITS",
                    f"3D*{session.year}-01-01*Active Space {session.year}*",
                    "AudibleTransits_object.pkl",
                )
            )
            if len(matches) == 0:
                session.status("No audible transits pkl file found")
                return
            listener = AudibleTransits.from_pickle(matches[0])

        tracks = listener.tracks
        if tracks.empty:
            session.status("Audible transits is empty, skipping")
            return

        session.status(f"{len(tracks)} audible transits")
        if len(tracks) > session.max_tracks:
            session.status("Too many, sampling")
            tracks = tracks.sample(session.max_tracks, random_state=4)
            session.status(f"Showing {len(tracks)} transits")

        tracks = tracks.to_crs(session.crs)

        actors = []
        for _, track in tracks.iterrows():
            polyline = create_polyline_3d(track["interp_geometry"])
            actor = session.plotter.add_mesh(
                polyline, color=style.audible_transits_color, point_size=2, line_width=2
            )
            actors.append(actor)

        def toggle(flag):
            for actor in actors:
                actor.SetVisibility(flag)

        self._widgets.add_labeled_checkbox(
            toggle,
            value=True,
            position=(10, 180),
            size=25,
            color_on="purple",
            label="transits",
        )

    def plot_tracks(
        self,
        source: TrackSource,
        start_date: str | None = None,
        end_date: str | None = None,
        *,
        annotation_polyline: Callable | None = None,
        sea_surface_polyline: Callable | None = None,
        add_track_line: Callable | None = None,
    ) -> None:
        """Plot causal tracks from GPS, ADSB, or MXAK AIS for the study window."""
        to_annotation = annotation_polyline or self._polylines.annotation_polyline
        to_sea_surface = sea_surface_polyline or self._polylines.sea_surface_polyline
        add_line = add_track_line or self._scene.add_track_line

        session = self._session
        style = session.style
        start_date = start_date or f"{session.year}-01-01"
        end_date = end_date or f"{session.year}-12-31"
        microphone = get_deployment(session.project_dir, session.unit, session.site, session.year)

        session.status(f"Querying {source} tracks from {start_date} to {end_date}")
        try:
            loaded = load_tracks(
                source,
                start_date=start_date,
                end_date=end_date,
                study_area=session.study_area,
                microphone=microphone,
                include_faa_paths=False,
            )
        except (configparser.NoSectionError, configparser.NoOptionError) as exc:
            session.status(f"No {source} tracks loaded: missing config option {exc!r}")
            return
        except KeyError as exc:
            session.status(f"No {source} tracks loaded: {exc}")
            return
        except (AssertionError, ValueError) as exc:
            session.status(f"No {source} tracks loaded: {exc}")
            return

        tracks = loaded.tracks.to_crs(session.crs)
        if tracks.empty:
            session.status(f"No {source} tracks loaded.")
            return

        track_ids = tracks["track_id"].drop_duplicates()
        session.status(f"{len(track_ids)} {source} tracks ({len(tracks)} points)")
        if len(track_ids) > session.max_tracks:
            session.status(f"More than {session.max_tracks}, sampling")
            selected = track_ids.sample(session.max_tracks, replace=False, random_state=3)
            tracks = tracks[tracks["track_id"].isin(selected)]
            session.status(f"Showing {selected.nunique()} tracks")

        color = (
            style.vessel_track_color
            if source is TrackSource.AIS
            else style.flight_track_color
        )
        actors = []
        for _track_id, group in tracks.groupby("track_id", sort=False):
            group = group.sort_values("point_dt")
            line = track_points_to_linestring(
                group,
                include_z=source is not TrackSource.AIS,
            )
            if line.is_empty or line.length == 0:
                continue
            match source:
                case TrackSource.AIS:
                    polyline = to_sea_surface(line)
                case TrackSource.ADSB | TrackSource.GPS:
                    polyline = to_annotation(line)
                case _:
                    raise ValueError(f"Unknown track source: {source}")
            actor = add_line(polyline, color=color)
            actors.append(actor)

        if not actors:
            session.status(f"No {source} tracks to plot.")
            return

        def toggle(flag):
            for actor in actors:
                actor.SetVisibility(flag)

        self._widgets.add_labeled_checkbox(
            toggle,
            value=True,
            position=(10, 20),
            size=25,
            color_on=color,
            label="tracks",
        )
