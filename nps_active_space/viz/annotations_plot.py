from __future__ import annotations

import pyvista as pv
from shapely.geometry import box
from tqdm import tqdm

from nps_active_space.utils.helpers import load_annotations
from nps_active_space.utils.models import Annotations
from nps_active_space.viz.annotations import format_annotation_summary
from nps_active_space.viz.geometry import iter_plot_linestrings
from nps_active_space.viz.polyline_builders import TrackPolylineBuilder
from nps_active_space.viz.scene import ScenePlotter
from nps_active_space.viz.session import VizSession
from nps_active_space.viz.widgets import PlotWidgets


class AnnotationsPlotter:
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

    def plot(self, annotation_file: str | None = None) -> None:
        session = self._session
        style = session.style
        if annotation_file is None:
            session.status(
                f"Loading annotations from project dir for {session.deployment} (valid only)"
            )
            annotations = load_annotations(
                session.project_dir, session.unit, session.site, session.year, only_valid=True
            )
        else:
            session.status(f"Loading annotations from {annotation_file} (valid only)")
            annotations = Annotations(annotation_file, only_valid=True)

        session.status(f"Parsed annotations: {format_annotation_summary(annotations)}")
        if annotations.empty:
            session.status("No annotations found, skipping.")
            return

        track_ids = annotations["_id"].drop_duplicates()
        if len(track_ids) > session.max_tracks:
            session.status(
                f"Sampling {session.max_tracks} of {len(track_ids)} tracks "
                f"({len(annotations)} segments before sample)"
            )
            selected_track_ids = track_ids.sample(session.max_tracks, replace=False, random_state=2)
            annotations = annotations[annotations["_id"].isin(selected_track_ids)]
            session.status(f"After sample: {format_annotation_summary(annotations)}")

        n_loaded = len(annotations)
        session.status(f"Reprojecting annotations to {session.crs} and clipping to study area")
        annotations = annotations.to_crs(session.crs)
        ann_bounds = annotations.total_bounds
        study_bounds = session.study_area.total_bounds
        session.status(f"Annotation bounds: {ann_bounds}")
        session.status(f"Study area bounds: {study_bounds}")

        n_empty = int(annotations.geometry.is_empty.sum())
        if n_empty:
            session.status(f"Dropping {n_empty} empty geometries before clip")
        annotations = annotations[~annotations.geometry.is_empty]
        annotations = annotations.clip(box(*study_bounds))
        annotations = annotations[~annotations.geometry.is_empty].explode(ignore_index=True)
        session.status(
            f"After clip: {len(annotations)} segments "
            f"({n_loaded - len(annotations)} removed outside study area)"
        )

        audible_actors = []
        inaudible_actors = []
        audible_polylines: list[pv.PolyData] = []
        inaudible_polylines: list[pv.PolyData] = []
        n_plotted = 0
        n_skipped = 0
        for _, annot in tqdm(
            annotations.iterrows(),
            total=len(annotations),
            desc="Building annotation lines",
            unit="segment",
        ):
            for line in iter_plot_linestrings(annot["geometry"]):
                polyline = self._polylines.annotation_polyline(line)
                if polyline.n_points < 2:
                    n_skipped += 1
                    continue
                n_plotted += 1
                if annot["audible"]:
                    audible_polylines.append(polyline)
                else:
                    inaudible_polylines.append(polyline)

        audible_actor = self._scene.add_annotation_lines(
            audible_polylines, color=style.audible_annotation_color
        )
        if audible_actor is not None:
            audible_actors.append(audible_actor)
        inaudible_actor = self._scene.add_annotation_lines(
            inaudible_polylines, color=style.inaudible_annotation_color
        )
        if inaudible_actor is not None:
            inaudible_actors.append(inaudible_actor)

        if n_plotted == 0:
            session.status(
                f"No annotation segments plotted ({n_loaded} loaded; "
                f"{n_skipped} degenerate lines skipped). "
                "Check CRS and study-area overlap."
            )
            return
        session.status(
            f"Plotted {n_plotted} annotation line(s) in "
            f"{len(audible_actors) + len(inaudible_actors)} mesh(es) "
            f"({len(audible_polylines)} audible, {len(inaudible_polylines)} inaudible; "
            f"{n_skipped} skipped)"
        )

        def toggle_audible(flag):
            for actor in audible_actors:
                actor.SetVisibility(flag)

        def toggle_inaudible(flag):
            for actor in inaudible_actors:
                actor.SetVisibility(flag)

        self._widgets.add_labeled_checkbox(
            toggle_audible,
            value=True,
            position=(10, 90),
            size=25,
            color_on="deepskyblue",
            label="audible",
        )
        self._widgets.add_labeled_checkbox(
            toggle_inaudible,
            value=True,
            position=(10, 55),
            size=25,
            color_on="red",
            label="inaudible",
        )
