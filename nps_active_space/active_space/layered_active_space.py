from __future__ import annotations

import logging
import numpy as np
import os
import geopandas as gpd
import glob
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from shapely.geometry import Point
from nps_active_space.utils.computation import normalize_point_density
from nps_active_space.utils.legacy_nmsim_paths import find_layer_geojson

logger = logging.getLogger(__name__)


class LayeredActiveSpace():
    def __init__(self, designator: str, layer_dirs: dict[int, str], study_area: gpd.GeoDataFrame,
                 gain: float | None = None, crs: str = "epsg:4326") -> None:
        self.designator: str = designator
        self.layer_dirs: dict[int, str] = dict(sorted(layer_dirs.items()))
        self.study_area: gpd.GeoDataFrame = study_area
        self.activespaces: dict[int, gpd.GeoDataFrame] = {}
        self.all_activespaces: dict[float, dict[int, gpd.GeoDataFrame] | None] = {}  # set by self.preload_all_activespaces()
        self.crs: str = crs
        self.gain: float | None = None
        if gain is not None:
            self.set_gain(gain)
        self.fit_pbar: tqdm | None = None

        if not self.layer_dirs:
            self.gain_values: list[float] = []
            self.min_gain = 0.0
            self.max_gain = 0.0
            self.activespaces = None
            if gain is not None:
                self.gain = gain
            return

        self.gain_values = self._discover_gain_values()
        if self.gain_values:
            self.min_gain: float = self.gain_values[0]
            self.max_gain: float = self.gain_values[-1]
        else:
            self.min_gain = 0.0
            self.max_gain = 0.0

    @staticmethod
    def _gains_in_layer_dir(layer_dir: str) -> set[float]:
        from nps_active_space.utils.helpers import omni_to_gain

        active_names = glob.glob(os.path.join(layer_dir, "*_O_*.geojson"))
        return {omni_to_gain(path) for path in active_names}

    def _discover_gain_values(self) -> list[float]:
        """Omni gains present on every altitude layer (intersection across layers)."""
        if not self.layer_dirs:
            return []
        common: set[float] | None = None
        for layer_dir in self.layer_dirs.values():
            layer_gains = self._gains_in_layer_dir(layer_dir)
            common = layer_gains if common is None else common & layer_gains
        return sorted(common or [])
    
    def load_activespaces(self, gain: float) -> dict[int, gpd.GeoDataFrame] | None:
        if gain in self.all_activespaces:
            return self.all_activespaces[gain]

        activespaces: dict[int, gpd.GeoDataFrame] = {}
        for altitude, layer_dir in self.layer_dirs.items():
            activespace_file = find_layer_geojson(layer_dir, gain)
            if activespace_file is None:
                logger.warning(
                    "No active space for gain %s in %s; skipping altitude %s m",
                    gain, layer_dir, altitude,
                )
                continue
            activespaces[altitude] = gpd.read_file(activespace_file).to_crs(self.crs)

        if not activespaces:
            return None
        return activespaces

    def preload_all_activespaces(self) -> None:
        self.all_activespaces = {}
        for gain in tqdm(self.gain_values, desc="Loading all active spaces"):
            self.all_activespaces[gain] = self.load_activespaces(gain)

    def set_gain(self, gain: float) -> None:
        self.activespaces = self.load_activespaces(gain)
        self.gain = gain

    def fit(self, annotations: gpd.GeoDataFrame, beta: float = 1., plot: bool = True,
            plot_savepath: str | None = None) -> pd.Series:

        # Extract all valid points from their LineStrings. These will be needed for calculating fbeta scores later.
        valid_points_lst = []
        for idx, row in tqdm(annotations.iterrows(), total=annotations.shape[0], desc='Extracting valid points', unit=' valid track', colour='white'):
            valid_points_lst.extend([{'annotation_idx': idx, 'audible': row.audible, 'geometry': Point(coords)} for coords in row.geometry.coords])
        valid_points = gpd.GeoDataFrame(data=valid_points_lst, geometry='geometry', crs=annotations.crs)

        result = self.fit_points(valid_points, self.min_gain, self.max_gain, beta, plot, plot_savepath)
        result["Number of valid annotated segments"] = len(annotations)
        return result

    def fit_points(self, points: gpd.GeoDataFrame, min_gain: float = -10., max_gain: float = 40.,
                   beta: float = 1., plot: bool = True, plot_savepath: str | None = None) -> pd.Series:
        assert min_gain <= max_gain

        # Reduce point density to median density, so very dense areas (e.g. airports) don't skew the fit
        points_before_kde = len(points)
        points = normalize_point_density(points, self.study_area, random_seed=679)
        points_after_kde = len(points)

        # Convert to activespace CRS
        if points.crs != self.crs:
            points = points.to_crs(self.crs)
        
        logger.info("Assigning points to their closest layer")
        points = self.assign_layers(points)

        gain_values = self._resolve_fit_gain_values(min_gain, max_gain)
        if not gain_values:
            raise ValueError(
                f"No omni gain geojsons found on all layers for {self.designator} "
                f"(discovered {self.gain_values!r})",
            )

        # Compute precision, recall, and F-Beta for each gain on disk
        logger.info(
            "Computing performance for gains %s dB (%d values)",
            gain_values,
            len(gain_values),
        )
        detection_results = pd.DataFrame([])
        self.fit_pbar = tqdm(gain_values, unit=" Gain Values")
        for gain in self.fit_pbar:
            desc = f"{gain}dB, loading actives"
            self.fit_pbar.set_description(f"{desc:<25}")

            self.set_gain(gain)
            if not self.activespaces:
                logger.warning("Skipping gain %s dB: no active-space layers loaded", gain)
                continue

            in_AS = self.predict(points)
            audible = points["audible"]
            TP = np.all([in_AS, audible], axis=0).sum()
            FP = np.all([in_AS, ~audible], axis=0).sum()
            FN = np.all([~in_AS, audible], axis=0).sum()
            if TP == 0:
                fbeta = 0
                precision = 0
                recall = 0
            else:
                precision = TP / (TP + FP)  # specificity... if a flight enters the active space, is it actually audible?
                recall = TP / (TP + FN)  # sensitivity... if a flight is audible, does it enter the active space?
                fbeta = (1 + np.power(beta, 2)) * ((precision * recall) / ((np.power(beta, 2) * precision) + recall))

            detection_results.loc[gain, f"F{beta}"] = fbeta
            detection_results.loc[gain, "Precision"] = precision
            detection_results.loc[gain, "Recall"] = recall
        self.fit_pbar.close()

        if detection_results.empty:
            raise ValueError(
                f"No F{beta} scores computed for {self.designator}; "
                "check ACTIVESPACES geojson naming matches omni stems (e.g. O_-100 = -10 dB)",
            )

        best_gain = detection_results[f"F{beta}"].idxmax()
        best = detection_results.loc[best_gain]
        logger.info(f"Best gain: {best_gain}dB, F{beta} = {best[f"F{beta}"]}")

        self.set_gain(best_gain)

        if plot:
            # create Precision-Recall Plot.
            fig, ax = plt.subplots()
            ax.plot(detection_results["Recall"], detection_results["Precision"], ls="-", lw=0.2, marker="o", ms=2, color="k")
            ax.plot(best["Recall"], best["Precision"], ls="", marker="o", ms=5, color="None", markeredgecolor="limegreen")
            ax.text(best["Recall"], best["Precision"], f"  {best_gain}dB\n  F{beta}={best[f"F{beta}"]:.3f}", fontsize=8)
            ax.set_title(f'Precision-Recall Curve, F{beta:.1f}', loc="left")
            ax.set_ylabel('Precision')
            ax.set_xlabel('Recall')

            if plot_savepath is None:
                plt.show()
            else:
                os.makedirs(os.path.dirname(plot_savepath), exist_ok=True)
                plt.savefig(plot_savepath)
        
        return pd.Series({
            "Designator": self.designator,
            "KDE reduction (%)": f"{100 * (1 - (points_after_kde / points_before_kde))}%",
            "1/3rd Octave Gain (F1)": best_gain,
            f"F{beta}": best[f"F{beta}"]
        })

    def _resolve_fit_gain_values(
        self,
        min_gain: float,
        max_gain: float,
    ) -> list[float]:
        """Gains to evaluate during fit: on-disk ladder clipped to [min_gain, max_gain]."""
        return [
            g for g in self.gain_values
            if min_gain <= g <= max_gain
        ]

    def assign_layers(self, points: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Returns a copy of points with a new column added representing the active space layer
        each point belongs to (which layer is closest to the point's z value)
        """
        points = points.copy()
        altitudes = list(self.layer_dirs.keys())
        def closest_layer(z: float) -> int:
            return min(altitudes, key=lambda alt: abs(alt - z))
        points["layer"] = points.geometry.z.apply(closest_layer)
        return points


    def predict(self, points: gpd.GeoDataFrame) -> pd.Series:
        """Given a GeoDataFrame of 3D points, predict whether they are audible.
        
        Returns
        -------
        predictions: pd.Series
            Boolean series representing whether the model predicts each point is audible.
        """
        if not self.activespaces:
            raise ValueError("Can't predict with no activespaces")
        assert points.crs == self.crs

        if "layer" not in points.columns:
            points = self.assign_layers(points)
        
        # Determine if points are inside their layer's activespace
        points["in_AS"] = False
        for altitude, activespace in self.activespaces.items():
            desc = f"{self.gain}dB, {altitude}m"
            if self.fit_pbar is not None:
                self.fit_pbar.set_description(f"{desc:<25}")
            layer_mask = points["layer"] == altitude
            in_AS_gdf = gpd.clip(points[layer_mask], activespace)
            points.loc[in_AS_gdf.index, "in_AS"] = True

        return points["in_AS"]
