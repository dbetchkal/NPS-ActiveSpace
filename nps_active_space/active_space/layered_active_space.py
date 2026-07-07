from typing import Dict, Optional

import numpy as np
import os
import geopandas as gpd
import glob
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from shapely.geometry import Point
from nps_active_space.utils.computation import normalize_point_density


class LayeredActiveSpace():
    def __init__(self, designator: str, layer_dirs: Dict[int, str], study_area: gpd.GeoDataFrame,
                 gain: Optional[float] = None, crs: str = "epsg:4326") -> None:
        self.designator: str = designator
        self.layer_dirs: Dict[int, str] = dict(sorted(layer_dirs.items()))
        self.study_area: gpd.GeoDataFrame = study_area
        self.activespaces: Dict[int, gpd.GeoDataFrame] = {}
        self.all_activespaces: Dict[float, Optional[Dict[int, gpd.GeoDataFrame]]] = {}  # set by self.preload_all_activespaces()
        self.crs: str = crs
        self.gain: Optional[float] = None
        if gain is not None:
            self.set_gain(gain)
        self.fit_pbar: Optional[tqdm] = None

        # determine min and max gain - import here to avoid circular import
        from nps_active_space.utils.helpers import omni_to_gain
        first_layer_dir = list(self.layer_dirs.values())[0]
        active_names = glob.glob(os.path.join(first_layer_dir, "*_O_*.geojson"))
        gains = list(map(lambda f: omni_to_gain(f), active_names))
        self.min_gain: float = min(gains)
        self.max_gain: float = max(gains)
    
    def load_activespaces(self, gain: float) -> Optional[Dict[int, gpd.GeoDataFrame]]:
        if gain in self.all_activespaces:
            return self.all_activespaces[gain]
        
        sign = "-" if gain < 0 else "+"
        gain_string = str(np.abs(int(10*gain))).zfill(3)

        activespaces: Dict[int, gpd.GeoDataFrame] = {}
        for altitude, dir in self.layer_dirs.items():
            glob_result = glob.glob(os.path.join(dir, f"*_O_{sign}{gain_string}.geojson"))
            if len(glob_result) == 0:
                print(f"Couldn't find active space for gain {gain} in {dir}")
                return None
            activespace_file = glob_result[0]
            activespaces[altitude] = gpd.read_file(activespace_file).to_crs(self.crs)
        
        return activespaces

    def preload_all_activespaces(self) -> None:
        self.all_activespaces = {}
        for gain in tqdm(np.arange(self.min_gain, self.max_gain + 0.5, 0.5), desc="Loading all active spaces"):
            self.all_activespaces[gain] = self.load_activespaces(gain)

    def set_gain(self, gain: float) -> None:
        self.activespaces = self.load_activespaces(gain)
        self.gain = gain

    def fit(self, annotations: gpd.GeoDataFrame, beta: float = 1., plot: bool = True,
            plot_savepath: Optional[str] = None) -> pd.Series:

        # Extract all valid points from their LineStrings. These will be needed for calculating fbeta scores later.
        valid_points_lst = []
        for idx, row in tqdm(annotations.iterrows(), total=annotations.shape[0], desc='Extracting valid points', unit=' valid track', colour='white'):
            valid_points_lst.extend([{'annotation_idx': idx, 'audible': row.audible, 'geometry': Point(coords)} for coords in row.geometry.coords])
        valid_points = gpd.GeoDataFrame(data=valid_points_lst, geometry='geometry', crs=annotations.crs)

        result = self.fit_points(valid_points, self.min_gain, self.max_gain, beta, plot, plot_savepath)
        result["Number of valid annotated segments"] = len(annotations)
        return result

    def fit_points(self, points: gpd.GeoDataFrame, min_gain: float = -10., max_gain: float = 40.,
                   beta: float = 1., plot: bool = True, plot_savepath: Optional[str] = None) -> pd.Series:
        assert min_gain <= max_gain

        # Reduce point density to median density, so very dense areas (e.g. airports) don't skew the fit
        points_before_kde = len(points)
        points = normalize_point_density(points, self.study_area, random_seed=679)
        points_after_kde = len(points)

        # Convert to activespace CRS
        if points.crs != self.crs:
            points = points.to_crs(self.crs)
        
        print("Assigning points to their closest layer")
        points = self.assign_layers(points)

        # Compute precision, recall, and F-Beta for each gain
        print(f"Computing performance for each gain {min_gain}dB to {max_gain}dB")
        detection_results = pd.DataFrame([])
        self.fit_pbar = tqdm(np.arange(min_gain, max_gain + 0.5, 0.5), unit=" Gain Values")
        for gain in self.fit_pbar:
            desc = f"{gain}dB, loading actives"
            self.fit_pbar.set_description(f"{desc:<25}")

            self.set_gain(gain)

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

        best_gain = detection_results[f"F{beta}"].idxmax()
        best = detection_results.loc[best_gain]
        print(f"Best gain: {best_gain}dB, F{beta} = {best[f"F{beta}"]}")

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
        assert len(self.activespaces) > 0, "Can't predict with no activespaces"
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
