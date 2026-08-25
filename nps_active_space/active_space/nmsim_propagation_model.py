"""NMSim implementation of :class:`PropagationModel`."""

from __future__ import annotations

import glob
import logging
import os
import subprocess
from dataclasses import dataclass
from typing import Optional
from uuid import uuid4

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.setup.site_writer import write_listener_site_file
from nps_active_space.utils.constants import IS_WINDOWS
from nps_active_space.utils.computation import coords_to_utm, project_raster
from nps_active_space.active_space.propagation_model import (
    NMSIM_PREDICTIONS_SUBDIR,
    NMSIM_SCRATCH_SUBDIR,
)

logger = logging.getLogger(__name__)


def nmsim_control_path(path: str) -> str:
    """Format a filesystem path for an NMSIM control file (Wine backslashes on Linux/Mac)."""
    if IS_WINDOWS:
        return path
    return "Z:" + os.path.abspath(path).replace("/", "\\")


@dataclass(frozen=True)
class NmsimSiteContext:
    dem_file: str
    flt_file: str
    site_file: str


class NmsimPropagationModel:
    max_points_per_run = 4000
    predictions_subdir = NMSIM_PREDICTIONS_SUBDIR

    def __init__(self, nmsim_exe: str, root_dir: str) -> None:
        assert os.path.exists(nmsim_exe), "NMSIM not found"
        self.nmsim_exe = nmsim_exe
        self.root_dir = root_dir

    def prepare_site(
        self,
        dem_src: str,
        study_area: gpd.GeoDataFrame,
        mic: Microphone,
        *,
        project_dem: bool = True,
        suffix: str = "",
    ) -> NmsimSiteContext:
        dem_file = self._mask_dem_file(dem_src, study_area, project=project_dem, suffix=suffix)
        flt_file = self._create_dem_flt(dem_file)
        site_file = self._create_site_file(mic, flt_file)
        return NmsimSiteContext(dem_file=dem_file, flt_file=flt_file, site_file=site_file)

    def predict(
        self,
        site: NmsimSiteContext,
        source_pts: gpd.GeoDataFrame,
        omni_source: str,
        altitude_m: int,
        job_name: str,
        heading: int | None = None,
    ) -> pd.DataFrame:
        crs = source_pts.crs.to_string().lower()
        trajectory_filename = self._create_trajectory_file(source_pts, crs, job_name, heading)
        batch_file = self._create_instruction_files(
            site.flt_file, site.site_file, trajectory_filename, omni_source,
        )
        process = subprocess.Popen(
            [self.nmsim_exe, batch_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        stdout, stderr = process.communicate()
        if stderr:
            for line in stderr.decode("utf-8").splitlines():
                logger.error(line.strip())
        return self._postprocess_trj_tis(
            trajectory_filename,
            f"{self.root_dir}/{NMSIM_SCRATCH_SUBDIR}/{job_name}.tis",
            cleanup=True,
        )

    def _mask_dem_file(
        self,
        dem_src: str,
        study_area: gpd.GeoDataFrame,
        project: bool = False,
        buffer: Optional[int] = None,
        suffix: str = "",
    ) -> str:
        if project:
            dem_projected_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation{suffix}.tif"
            project_raster(dem_src, dem_projected_filename, study_area.crs)
            dem_src = dem_projected_filename

        study_area_filename_prefix = (
            f"{self.root_dir}/Input_Data/01_ELEVATION/study_area{suffix}_{uuid4()}"
        )
        study_area_filename = f"{study_area_filename_prefix}.shp"
        if buffer:
            equal_area_crs, _ = coords_to_utm(
                study_area.centroid.iat[0].y, study_area.centroid.iat[0].x,
            )
            study_area_m = study_area.to_crs(equal_area_crs)
            study_area_m = study_area_m.buffer(buffer * 1000)
            study_area = study_area_m.to_crs(study_area.crs)

        if "FID" in study_area.columns:
            study_area = study_area.drop(columns=["FID"])

        study_area.to_file(study_area_filename)

        dem_masked_filename = f"{self.root_dir}/Input_Data/01_ELEVATION/elevation{suffix}_masked.tif"
        # Pass the path to Warp: GDAL Dataset is not a context manager on the
        # versions we ship in Docker (osgeo.gdal < 3.8).
        gdal.Warp(
            dem_masked_filename,
            dem_src,
            cutlineDSName=study_area_filename,
            cropToCutline=True,
            dstNodata=-9999,
        )

        for filename in glob.glob(f"{study_area_filename_prefix}*"):
            os.remove(filename)

        return dem_masked_filename

    def _create_dem_flt(self, dem_file: str) -> str:
        flt_filename = dem_file.replace(".tif", ".flt")
        flt_header_filename = dem_file.replace(".tif", ".hdr")

        gdal.Translate(flt_filename, dem_file, options="-ot Float32 -of ehdr -a_nodata -9999")

        old_hdr = pd.read_csv(flt_header_filename, header=None, sep=r"\s+", index_col=0).T
        yllcorner = float(old_hdr.ULYMAP) - float(old_hdr.NROWS) * float(old_hdr.XDIM)

        with open(flt_header_filename, "w") as header:
            header.write("{:14}{:}\n".format("ncols", old_hdr.NCOLS.values[0]))
            header.write("{:14}{:}\n".format("nrows", old_hdr.NROWS.values[0]))
            header.write("{:14}{:}\n".format("xllcorner", old_hdr.ULXMAP.values[0]))
            header.write("{:14}{:}\n".format("yllcorner", yllcorner))
            header.write("{:14}{:}\n".format("cellsize", old_hdr.XDIM.values[0]))
            header.write("{:14}{:}\n".format("NODATA_value", old_hdr.NODATA.values[0]))
            header.write("{:14}{:}".format("byteorder", "LSBFIRST"))

        return flt_filename

    def _create_site_file(self, mic: Microphone, dem_file: str) -> str:
        return str(write_listener_site_file(
            self.root_dir, mic.name, mic.x, mic.y, mic.z, dem_file,
        ))

    def _create_trajectory_file(
        self,
        points: gpd.GeoDataFrame,
        crs: str,
        filename: str,
        heading: Optional[int] = None,
    ) -> str:
        trajectory_filename = f"{self.root_dir}/Input_Data/03_TRAJECTORY/{filename}.trj"

        trajectory = points.to_crs(crs)
        trajectory["heading"] = (
            heading if heading is not None
            else np.random.choice(range(0, 360), size=len(points), replace=True)
        )
        trajectory["climb_angle"] = 0
        trajectory["power"] = 95
        trajectory["rol"] = 0

        velocity = 70
        trajectory["knots"] = 1.94384 * velocity

        dist = np.diff(
            np.array([
                [x, y, z] for x, y, z in zip(
                    trajectory.geometry.x, trajectory.geometry.y, trajectory.geometry.z,
                )
            ]),
            axis=0,
        )
        time_elapsed = np.cumsum(np.array([np.linalg.norm(d) for d in dist]) / velocity)
        time_elapsed = np.append(time_elapsed, np.nan)
        trajectory["time_elapsed"] = time_elapsed

        with open(trajectory_filename, "w") as trajectory_file:
            trajectory_file.write("Flight track trajectory variable description:\n")
            trajectory_file.write(" time - time in seconds from the reference time\n")
            trajectory_file.write(" Xpos - x coordinate (UTM)\n")
            trajectory_file.write(" Ypos - y coordinate (UTM)\n")
            trajectory_file.write(" UTM Zone  " + crs + "\n")
            trajectory_file.write(" Zpos - z coordinate in meters MSL\n")
            trajectory_file.write(" heading - aircraft compass bearing in degrees\n")
            trajectory_file.write(" climbANG - aircraft climb angle in degrees\n")
            trajectory_file.write(" vel - aircraft velocity in knots\n")
            trajectory_file.write(" power - % engine power\n")
            trajectory_file.write(" roll - bank angle (right wing down), degrees\n")
            trajectory_file.write("FLIGHT " + filename + "\n")
            trajectory_file.write("TEMP.  59.0\n")
            trajectory_file.write("Humid.  70.0\n")
            trajectory_file.write("\n")
            trajectory_file.write(
                "         time(s)        Xpos           Ypos           Zpos         heading        climbANG       Vel            power          rol\n",
            )

            for _, point in trajectory.iterrows():
                trajectory_file.write(
                    "{0:15.3f}".format(point.time_elapsed)
                    + "{0:15.3f}".format(point.geometry.x)
                    + "{0:15.3f}".format(point.geometry.y)
                    + "{0:15.3f}".format(point.geometry.z)
                    + "{0:15.3f}".format(point.heading)
                    + "{0:15.3f}".format(point.climb_angle)
                    + "{0:15.3f}".format(point.knots)
                    + "{0:15.3f}".format(point.power)
                    + "{0:15.3f}".format(point.rol) + "\n"
                )

        return trajectory_filename

    def _create_instruction_files(
        self,
        flt_file: str,
        site_file: str,
        trajectory_file: str,
        omni_source_file: str,
    ) -> str:
        control_file = (
            f"{self.root_dir}/control_{os.path.basename(trajectory_file).replace('.trj', '')}.nms"
        )
        batch_file = (
            f"{self.root_dir}/batch_{os.path.basename(trajectory_file).replace('.trj', '')}.txt"
        )
        tis_directory = f"{self.root_dir}/{NMSIM_SCRATCH_SUBDIR}"

        with open(control_file, "w") as nms:
            nms.write(nmsim_control_path(flt_file) + "\n")
            nms.write("-\n")
            nms.write(nmsim_control_path(site_file) + "\n")
            nms.write(nmsim_control_path(trajectory_file) + "\n")
            nms.write(nmsim_control_path(f"{ACTIVE_SPACE_DIR}/data/default.wea") + "\n")
            nms.write("-\n")
            nms.write(nmsim_control_path(omni_source_file) + "\n")
            nms.write("{0:11.4f}   \n".format(500.0000))
            nms.write("-\n")
            nms.write("-")

        with open(batch_file, "w") as batch:
            batch.write("open\n")
            batch.write(nmsim_control_path(control_file) + "\n")
            batch.write("site\n")
            batch.write(
                nmsim_control_path(
                    f"{tis_directory}/{os.path.basename(trajectory_file)[:-4]}",
                ) + "\n",
            )
            batch.write("dbf: no\n")
            batch.write("hrs: 0\n")
            batch.write("min: 0\n")
            batch.write("sec: 0.0")

        return batch_file

    def _postprocess_trj_tis(
        self,
        trajectory_file: str,
        tis_file: str,
        cleanup: bool = True,
    ) -> pd.DataFrame:
        traj_df = pd.read_fwf(trajectory_file, header=14, widths=[16, 14] + [15] * 7)
        traj_df = traj_df.drop(
            ["time(s)", "heading", "climbANG", "Vel", "power", "rol"], axis=1,
        )

        tis_df = pd.read_fwf(tis_file, header=15, skipfooter=1, widths=[4, 12] + [5] * 35)
        tis_df.drop(["SP#", "F", "42"], axis=1, inplace=True)
        tis_df.drop(tis_df.head(2).index, inplace=True)
        tis_df.reset_index(drop=True, inplace=True)
        tis_df.columns = [
            "TIME", "A", "10", "12.5", "15.8", "20", "25", "31.5", "40", "50", "63",
            "80", "100", "125", "160", "200", "250", "315", "400", "500", "630", "800", "1000",
            "1250", "1600", "2000", "2500", "3150", "4000", "5000", "6300", "8000", "10000", "12500",
        ]
        tis_df = tis_df.drop("TIME", axis=1)
        tis_df = (tis_df.astype(float) * 0.1).round(6)

        assert not tis_df.empty, (
            "NMSIM didn't run, try reducing max_pts in _preprocess_source_pts() and try again"
        )
        assert len(traj_df) == len(tis_df), (
            f"# trajectory points ({len(traj_df)}) is not equal to # tis points ({len(tis_df)})"
        )
        new_rows = pd.concat([traj_df, tis_df], axis=1)

        if cleanup:
            os.remove(trajectory_file)
            os.remove(tis_file)

        return new_rows
