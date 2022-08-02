# -*- coding: utf-8 -*-
"""
Initializes a 'CreateActiveSpace' object from either an existing NMSim project folder or a location in a temporary
directory. Runs NMSim from the batch script to identify the boundary between audibility and inaudibility and returns
a GeoDataFrame shape.

Audibility is decided from either:
    1) Lasting NVSPL record for sites from the inventory monitoring project with an audio record
    2) Dan Mennitt's broadband ambience model

Flight altitude is calculated from either:
    1) A user-specified value
    2) A default value (12,000 ft)
    3) The average flight altitude from flight tracks collected by overflights monitoring

Created 2021-07-06 11:30

@author: Kirby Heck
"""

# geospatial libraries
import shapely.errors
import tempfile
import warnings
import os
from shapely.ops import triangulate
from shapely.geometry import Polygon, MultiPoint

# other libraries
from matplotlib.tri import Triangulation
import sys

# DENA RDS computer
RDS = r"\\inpdenards\overflights"

sys.path.append(os.path.join(RDS, "scripts"))
from query_tracks import query_tracks

# specialized NPS libraries
cd = os.getcwd()
sys.path.append(os.path.join(cd, "iyore-master"))
import iyore

sys.path.append(os.path.join(cd, "soundDB-master"))
from soundDB import nvspl, srcid

sys.path.append(os.path.dirname(cd))  # is this necessary?
from active_space_utils import *


class CreateActiveSpace(object):

    def __init__(self, unit_name='DENA', site_name=None, yr=None,
                 root_dir=r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_ProjectDirectories",
                 project_dir=None,
                 source_path=os.path.join(cd, "NMSIM-Python", "NMSIM", "Sources", "AirTourFixedWingSources", "C207.src"),
                 NMSIM_path=os.path.join(cd, "NMSIM-Python", "NMSIM", "Nord2000batch.exe"),
                 ras_path=r"T:\ResMgmt\WAGS\Sound\GIS\AlaskaRegion\AK Soundscape Model\DENA_L50_Existing.tif",
                 save_pkl=False, printout=False, nvspl_max=None,
                 ambience_source='NVSPL', freq_bands=False, lxx=50,
                 coord=None, crs=None):
        """
        Initializes an instance of CreateActiveSpace object that can run the NMSim model.

        Parameters
        ----------
        unit_name (str): optional, four character park code ('DENA'). Default DENA if none provided
        site_name (str): optional, four character site name ('KAHP'). Default `None` (if none provided, this will fail
            to look up nvspl files and use broadband instead)
        yr (int): optional, default `None` (if none provided, this will fail to look up nvspl files and use broadband
            instead)
        root_dir (str): optional, directory for NMSIM project directories
        project_dir (str): optional, directory to the site-specific NMSim project directory. If included, overrides the
            root directory
        source_path (str): optional, file path for the aircraft source file
        NMSIM_path (str): optional, file path for the NMSIM batch analysis .exe
        ras_path (str): optional, file path for the parkwide ambience raster
        save_pkl (bool): optional, saves the active space to a pickle binary file. Default False
        printout (bool): optional, shows graphs and prints additional statements in the console while creating the
            active space. Default False
        nvspl_max (int): optional, maximum number of NVSPL files to load. Default None
        ambience_source (str): optional, default 'NVSPL'. Selects which kind of audibility test to assign to the output
            .tis files. Either 'NVSPL' or 'Mennitt'
        freq_bands (bool): optional, default False. If input as True, looks for a third octave band spectrum in the
            NVSPL record; requires ambience_source to be 'NVSPL'
        lxx (float): optional, quantile L_x level (between 0-100) for 3rd octave bands if comparing to NVSPL ambience.
            Default 50; median value
        coord (tuple): optional, coordinates of the site (long, lat) or (x, y)
        crs (str): required for reading in long/lat coordinates. Formatted like 'epsg:4326' for WGS84. It is
            recommended that the coordinate system is the same as the basemap and other raster files.
        """

        # ======== INSTANCE VARIABLES ========================

        if project_dir in {'tmp', 'temp', 'TMP', 'TEMP'}:
            self.tmp = tempfile.TemporaryDirectory()  # make a temporary project directory
            project_dir = self.tmp.name
        else:
            self.tmp = None

        if crs is not None:  # no site, assume long/lat are given
            self.in_crs = crs
            if coord is not None:
                self.coord = coord
                self.project_dir = project_dir

                # set the root directory based on the project directory unless none is provided
                if project_dir is None:
                    self.root_dir = None
                else:
                    self.root_dir = os.path.dirname(self.project_dir)

                self.yr = 9999  # filler
                self.site_name = '0000'  # filler
            else:
                raise ValueError('Longitude (x) and latitude (y) must be given in a respective coordinate system')

        elif site_name is not None:  # this is an inventory site in the park
            self.site_name = site_name
            self.unit_name = unit_name

            # also set lat/long of the site:
            site_info = get_site_info(unit_name=self.unit_name, site_name=site_name, year=yr)
            self.coord = (site_info.x.values[0], site_info.y.values[0])
            self.in_crs = site_info.utm_zone.values[0]

            # the year may change if get_site_info() cannot find a match
            self.yr = site_info.year.values[0]

            # only proj_dir OR root_dir are required for the inventoried sound station sites
            if project_dir is not None:
                self.set_proj_dir(project_dir)
            else:
                self.root_dir = root_dir
                self.project_dir = self.root_dir + os.sep + self.unit_name + self.site_name

        else:
            raise ValueError('Either a coordinate system with long/lat coordinates '
                             'or a site name and year must be given. ')

        self.source_path = source_path
        self.NMSIM_path = NMSIM_path
        self.ras_path = ras_path
        self.printout = printout  # show progress, text, and [many] plots
        self.save_pkl = save_pkl
        self.bands = freq_bands
        self.ambience_source = ambience_source
        self.nvspl_max = nvspl_max
        self.Lxx = lxx
        self.study_area = None  # set later in create_loc() if needed

        self.Lx = self._set_ambience()

    def __str__(self):
        """
        Overrides the string method to print out a description of the CreateActiveSpace object.
        """

        if self.site_name == '0000':
            to_str = "for (long, lat) ({:.2f}, {:.2f})".format(self.coord[0], self.coord[1])
        else:
            to_str = "for " + self.project_dir
        return "CreateActiveSpace object " + to_str

    def _set_ambience(self):
        """
        Loads an ambience level or spectrum into self.Lx from either the NVSPL record or Dan Mennitt's statewide
        L50 ambient levels

        Returns
        -------
        Lx (either dataframe or float): for third octave band ambience, returns a pandas dataframe with headers for the
            frequency bands and values [dB] for the Lxx quantile. For broadband, returns the broadband dBA.
        """

        if self.ambience_source == 'NVSPL':  # Load NVSPL dataset
            if self.yr == 2021:  # TODO: This will inevitably change
                archive = iyore.Dataset(r'V:\Data 2021')
            else:
                archive = iyore.Dataset(r"E:")

            nv = nvspl(archive, unit=self.unit_name, site=self.site_name, year=self.yr, n=self.nvspl_max).combine()

            if len(nv.index.names) > 1:
                nv.index = nv.index.droplevel()

            if self.bands:  # load third octave bands
                Lx = nv.loc[:, "12.5":"20000"].quantile(1 - (self.Lxx / 100))

                if len(Lx) != 33:
                    print("\t_set_ambience(): Something went wrong loading the frequency bands... \
                          try removing '∞' strings? ")
                    # in some sites (TEK4...?) there are '∞' strings in the NVSPL data! remove these...
                    nv[nv == '∞'] = np.nan  # try and get rid of strings
                    Lx = nv.loc[:, "12.5":"20000"].astype(float).quantile(1 - (self.Lxx / 100))

            else:  # load broadband
                Lx = nv.loc[:, 'dbA'].quantile(1 - (self.Lxx/100))

        else:  # not using NVSPL, read the broadband L50 from an input raster
            if self.ambience_source != 'Mennitt':
                warnings.warn("ambience_source should be 'Mennitt' or 'NVSPL'. Loading Mennitt model")

            if self.Lxx != 50:
                warnings.warn("Switching to L50 median sound levels for Mennitt's ambience model.")
                self.Lxx = 50  # Mennitt's model only has L50 levels

            if self.bands:
                warnings.warn("Switching to broadband sound level for Mennitt's ambience model.")
                self.bands = False

            # open the Mennitt raster and pick out the location's broadband L50
            with rasterio.open(self.ras_path) as ras:
                if self.in_crs != ras.crs:
                    in_proj = pyproj.CRS.from_string(self.in_crs)

                    # convert into the raster's coordinate system
                    input_to_m = pyproj.Transformer.from_crs(in_proj, ras.crs)
                    x, y = input_to_m.transform(self.coord[0], self.coord[1])  # apply projection
                else:
                    x, y = self.coord[0], self.coord[1]

                # now, look up the broadband level at that location
                band1 = ras.read(1)
                Lx = band1[ras.index(x, y)]

        return Lx

    def set_ambience(self, nvspl_max=None, ambience_source=None, freq_bands=None, lxx=None):
        """
        Updates the ambience level (Lx) with the given input arguments. Calls _set_ambience()

        Parameters
        ----------
        nvspl_max (int): optional, maximum number of NVSPL files to load. Default None
        ambience_source (str): optional, default 'NVSPL'. Selects which kind of audibility test to assign to the output
            .tis files. Either 'NVSPL' or 'Mennitt'
        freq_bands (bool): optional, default None. If input as True, looks for a third octave band spectrum in the
            NVSPL record; requires ambience_source to be 'NVSPL'
        lxx (float): optional, quantile L_x level (between 0-100) for 3rd octave bands if comparing to NVSPL ambience.
            Default 50; median value
        """

        if ambience_source is not None:
            self.ambience_source = ambience_source

        if nvspl_max is not None:
            self.nvspl_max = nvspl_max

        if freq_bands is not None:
            self.bands = freq_bands

        if lxx is not None:
            self.Lxx = lxx

        self.Lx = self._set_ambience()

    def create_loc(self, **kwargs):
        """
        Creates a DEM for the site coordinates and calls set_loc() to set the location of the site. Also sets
        self.study_area to the appropriate GeoDataFrame

        Parameters
        ----------
        see CreateActiveSpace.set_loc()
        """

        # Creating a location will require the correct computer infrastructure; try to create a directory
        make_NMSIM_project_dir(self.project_dir)  # create directory

        if len(kwargs) > 0:
            self.set_loc(**kwargs)
        self._write_DEM()

    def _write_DEM(self):
        """
        Internal function; writes a DEM of the newly clipped study area in the ModelManager project directory.

        Uses self.coord for positional data

        Returns
        -------
        mask (geodataframe): rectangular mask that fits the clipped DEM
        """

        # if we're in a temporary directory, pass that into clip_DEM()
        if self.tmp is not None:
            tmp_dir = self.project_dir
        else:
            tmp_dir = None

        # write DEM in .tif and .flt formats, requires gdal_transform.exe
        out_path = clip_DEM(self.coord, self.in_crs, proj_dir=self.project_dir, tmp_dir=tmp_dir)
        create_NMSIM_flt(out_path)

        self.study_area = get_raster_extent(out_path)

    def set_loc(self, coord=None, crs=None, unit_name=None, site_name=None, yr=None):
        """
        Updates the location based on longitude/latitude/coordinate system and calls _set_ambience() to update the
        ambient levels as well.

        Parameters
        ----------
        coord (tuple): coordinate (long, lat) or (x,y) for the site
        crs (str): required for reading in long/lat coordinates. Formatted like 'epsg:4326' for WGS84. It is
            recommended that the coordinate system is the same as the basemap and other raster files.
        """

        if coord is not None and crs is not None:
            self.coord = coord
            self.in_crs = crs
        elif site_name is not None and yr is not None:
            self.site_name = site_name
            self.yr = yr

            # also set lat/long of the site:
            if unit_name is None:
                unit_name = self.unit_name
            else:
                self.unit_name = unit_name

            site_info = get_site_info(unit_name=unit_name, site_name=site_name, year=yr)
            self.coord = (site_info.x.values[0], site_info.y.values[0])
            self.in_crs = site_info.utm_zone.values[0]
        else:
            warnings.warn('Insufficient information provided (need coordinate+crs OR [unit]+site+year). \
                        Location unchanged')
            return None

        # important! also set the ambience levels; location changed.
        self.Lx = self._set_ambience()

    def set_proj_dir(self, project_dir):
        self.project_dir = project_dir
        self.root_dir = os.path.dirname(self.project_dir)

    def set_source(self, source_str):
        """
        Changes the NMSim source filepath.

        Inputs
        ------
        source_str (string): new (absolute) path for the source file
        """
        if os.path.exists(source_str):
            self.source_path = source_str

        else:
            warnings.warn('Could not find source path; source unchanged.')

    def _get_flight_alt(self):
        """
        Queries overflight tracks and returns the average flight altitude over the point's study area.

        Returns
        -------
        altitude_ft (float) average flight altitude in [ft]

        Throws
        ------
        ValueError if there is no study area for the site
        IndexError if there are no flights in the study area
        """

        if self.study_area is None:
            raise ValueError('No study area associated with the site.')

        area_tracks = query_tracks(connection_txt=r"\\inpdenards\overflights\config\connection_info.txt",
                                   start_date='2017-01-01', end_date=dt.date.today(),
                                   mask=self.study_area, aircraft_info=False)

        if len(area_tracks) > 0:
            alt_ft = area_tracks.altitude_ft.mean()
            return alt_ft
        else:
            raise ValueError('No tracks in study area. ')

    def find_audible(self, prev_space=None):
        """
        Loads the .tis results from NMSim in the project directory and appends audible points to
        a GeoDataFrame

        Inputs
        ------
        prev_space (GeoDataFrame): previous active space gdf to append to, if it exists

        Outputs
        -------
        total_space (GeoDataFrame): GeoDataFrame of the current active + inactive space
        """

        trj_tis_paths = pair_trj_to_tis_results(self.project_dir)

        audible_points = []
        inaudible_points = []
        for trj_path, tis_path in trj_tis_paths:

            # load the .trj file as written
            current_trj = trj_reader(trj_path)

            # load the .tis file as written; format to decibels
            current_tis = tis_resampler(tis_path, utc_offset=-8)

            # select the frequency bands we have Lx data for...
            # Check to see if any of the frequency bands are louder than ambient
            if self.bands:
                audible_times = (current_tis.loc[:, "12.5":"12500"] > self.Lx["12.5":"12500"].values).sum(axis=1)
            else:
                audible_times = current_tis.loc[:, "A"] > self.Lx

            # look up the positions of the aircraft that produce audible sound
            # add these to a list of dataframes
            audible_points.append(current_trj.loc[current_tis[audible_times > 0].index, ["Xpos", "Ypos"]])
            inaudible_points.append(current_trj.loc[current_tis[audible_times == 0].index, ["Xpos", "Ypos"]])

            if self.printout:
                # prints the heading to the console
                curr_heading = current_trj["heading"]
                print("Heading: ", curr_heading.iloc[0])

                # format the site files into spectrograms
                result_spect = current_tis.loc[:, "12.5":"12500"]
                # visualize the results "on the fly"
                plt.figure(figsize=(10, 1))
                plt.imshow(result_spect.T, origin="lower", aspect="auto",
                           cmap="plasma", vmin=-10, vmax=50)
                plt.show()
                print((current_tis.loc[audible_times.index, "12.5":"12500"] > 0).sum())
                print("\n")

        # converts lists to DataFrames and strings to floats
        audible_points = pd.concat(audible_points).astype(float)
        audible_points["audible"] = 1
        inaudible_points = pd.concat(inaudible_points).astype(float)
        inaudible_points["audible"] = 0

        _, utm_zone = get_extent(self.project_dir)

        # convert audible points result into a geodataframe
        # TODO: how to handle duplicate points with different results?  (if n_tracks ≠ 1)
        active_space = gpd.GeoDataFrame(audible_points, crs=utm_zone,
                                        geometry=[Point(x, y) for x, y in zip(audible_points["Xpos"].astype('float'),
                                                                              audible_points["Ypos"].astype('float'))])
        null_space = gpd.GeoDataFrame(inaudible_points, crs=utm_zone,
                                      geometry=[Point(x, y) for x, y in zip(inaudible_points["Xpos"].astype('float'),
                                                                            inaudible_points["Ypos"].astype('float'))])

        # the goal is to build a total space (inaudible + audible)
        total_space = active_space.append(null_space, ignore_index=True)
        total_space = total_space.append(prev_space, ignore_index=True)

        return total_space

    def run_model(self, altitude_ft=12000, n_tracks=1, density=48, out_crs='epsg:4326', n_iter=3, simplify=100):
        """
        Compute the active space for a given flight altitude and source amplification over the location described in
        the CreateActiveSpace instance. Begins by testing a uniform grid of points spanning the study area; after 1 or 2
        iterations, switches to a Triangulated Irregular Network (TIN) approach of points to test.

        Parameters
        ----------
        altitude_ft (float): optional, flight altitude for NMSim trajectories. Default 12,000 ft.
            Set to 'flight' (str) to calculate the average flight altitude over that area (loads overflights)
        n_tracks (int): optional, number of trajectory tracks to test the model on; more tracks provides a robust active
            space for directional (anisotropic) sound sources. Default 1
        density (int): optional, point density of the grid-based meshing scheme. Default 48 [creates 48^2 = 2304 points]
        out_crs (str): optional, coordinate reference frame to output final geometry. Default 'epsg:4326' (WGS84)
        n_iter (int): optional, number of iterations to run. Default 3
        simplify (int): optional, smoothing done by Shapely to reduce the final filesize. Default is 100 [m]. Set to 0
            or False to turn off geometric simplification

        Returns
        -------
        final_poly (geodataframe): multipolygon geodataframe of the final active space shape
        OR None, if triangulation or contouring fails. This can occur if the source is too quiet (or loud).
        """

        # =============== Use NMSim to Compile Active Space =====================

        # try to calculate the average flight altitude if possible
        if altitude_ft == 'calc':
            try:
                altitude_ft = self._get_flight_alt()  # may throw an OperationalError or also Timeout
            except IndexError as e:
                print(e)
                print('Using default value for flight altitude: 12,000 ft.')
                altitude_ft = 12000

        altitude = int(np.round(altitude_ft / 3.28084))  # ft to meters
        velocity = 70  # m/s, this doesn't matter actually (?)

        print('run_model(): Beginning active space computations for ' + self.__str__() + ' at {:d} m\n'.format(altitude))

        # ================== let 'er rip ==================

        active_space = None
        total_space = None
        xy = None

        for i in range(0, n_iter):
            full_extent, utm_zone = get_extent(self.project_dir)

            # if this isn't the first iteration, active_space exists
            if active_space is not None:
                xmax, ymax, _ = active_space.max().astype(float)
                xmin, ymin, _ = active_space.min().astype(float)
                xpad = 0.2 * (xmax - xmin)  # pad extents by 20% on each side, 40% total
                ypad = 0.2 * (ymax - ymin)

                extent = ([xmin - xpad, xmax + xpad], [ymin - ypad, ymax + ypad])

                # how much smaller is the extent compared to the original study area?
                shrinkage = np.divide(np.diff(extent) - np.diff(full_extent), np.diff(full_extent))
                if (min(shrinkage) > -0.30) or (i > 1):
                    # if the largest lat/long change is more conservative than -30%, then move to contour refinement
                    print("run_model(): MESHING: Switching to contouring point selection...")
                    try:
                        xy = get_contour_xy(total_space)
                    except RuntimeError as e:
                        print(e)
                        print(self.__str__() + ": Triangulation failed")
                        plt.close('all')
                        return None

            else:
                extent = full_extent

            # create_trajectories
            try:
                track = create_trajectories(self.project_dir, self.yr, n_tracks, altitude, velocity,
                                            id_string="C-206_", coord=self.coord, crs=self.in_crs,
                                            density=density, heading="set", extent=extent, xy=xy)
            except IndexError as e:
                print(e)
                print(self.__str__() + ": create_trajectories was given an empty track")
                plt.close('all')
                return None

            # run model
            NMSIM_create_tis(self.project_dir, self.source_path,
                             NMSIMpath=self.NMSIM_path, printout=self.printout)

            # check points for audibility
            total_space = self.find_audible(prev_space=total_space)
            active_space = total_space.loc[total_space.audible == 1]

            print(f"run_model(): Complete: {i + 1} iteration(s) of active space. ")

            if self.printout:
                print(f"\tActive space after {i + 1} iteration(s): ")
                active_space.plot(zorder=10, markersize=3)
                plt.show()

        # ========== finish up with a nicer plot of the final geometry ============

        if not self.printout:
            # the proper way to do this is to call ax.tricontour's source code
            plt.close('all')  # cheap work-around for the contour plotting issue is to close all plots

        else:  # plot the triangulated active space grid
            as_multipoint = MultiPoint(active_space.geometry.values)
            c_hull = as_multipoint.convex_hull
            fig, ax = plt.subplots(figsize=(6, 6))

            p_hull = gpd.GeoSeries(c_hull)
            p_hull.plot(ax=ax, alpha=0.3)

            triangles = triangulate(active_space.geometry.unary_union)
            triangles_gdf = gpd.GeoDataFrame()
            triangles_gdf.geometry = triangles

            triangles_gdf.plot(ax=ax, zorder=5, edgecolor="red", alpha=0.2)

            active_space.plot(ax=ax, zorder=10, markersize=3)
            plt.show()
            # end plotting

        # ================== Create Multipolygon ======================

        total_space.sort_index(inplace=True)
        contour_data = total_space.values[:, :-1]  # slice Xpos, Ypos, audible columns
        fig, ax = plt.subplots()  # need an axis to call tricontour function

        levels = np.linspace(0, 1, 10, endpoint=False)
        x_tri = contour_data[:, 0]
        y_tri = contour_data[:, 1]
        z_tri = contour_data[:, 2]  # this is either a 1 or 0 ('audible' or 'inaudible', respectively)

        # create the triangulated irregular network and contour lines
        tri = mpl.tri.Triangulation(list(x_tri), list(y_tri))
        cs = ax.tricontour(tri, list(z_tri), levels=levels)

        # pick some contour level to choose as the active space boundary... 0.5 is somewhat arbitrary
        level_ind = np.where(cs.levels == 0.5)[0][0]
        level_val = cs.levels[level_ind]
        if self.printout:
            print("run_model(): using contour level =", level_val, "at level index", level_ind)

        else:
            plt.close('all')  # close triangulation figure

        # iterate through all contour paths in the line collection at level_ind
        for i, contour_path in enumerate(cs.collections[level_ind].get_paths()):
            x = contour_path.vertices[:, 0]
            y = contour_path.vertices[:, 1]
            new_poly = Polygon([(i[0], i[1]) for i in zip(x, y)])

            if i == 0:  # first polygon
                poly = new_poly
            else:  # append new polygons to the conglomerate if they are larger than some threshold area
                if new_poly.area > 50000:  # TODO: this area threshold may need to be changed!
                    try:
                        poly = poly.symmetric_difference(new_poly)

                    except shapely.errors.TopologicalError as e:
                        # this occurs if the 'active space' contours extend past the extents of the study area, which
                        # results in segments of polygons (but are indistinguishable from polygons)...
                        # the 'proper' solution would be to add a ring of zero entries before calculating contours to
                        # ensure each contour is a closed polygon
                        print(e)
                        print('Active space possibly too large (or the study area is too small).')
                        plt.close('all')
                        return None

        # put the final Multipolygon in a geodataframe
        final_poly = gpd.GeoSeries(poly, crs=utm_zone)

        if simplify:  # simplify geometry
            final_poly = final_poly.simplify(simplify)

        if final_poly.crs != out_crs:  # convert coordinate project if necessary
            final_poly = final_poly.to_crs(out_crs)

        # write the final polygon into a binary format that can be acquired by a pickle
        src_name = self.source_path.split(os.sep)[-1][:-4]  # gleans source filename and hacks off the extension
        if self.save_pkl:
            filename = self.project_dir + os.sep + self.unit_name + self.site_name + str(self.yr) + \
                       "_src" + src_name + "_" + str(altitude_ft) + "ft_active_space.pkl"
            filewriter = open(filename, 'wb')
            pickle.dump(final_poly, filewriter)
            filewriter.close()

        if self.printout:
            fig, ax = plt.subplots()
            final_poly.plot(ax=ax)
            plt.show()

        print("\nActive space complete for " + self.unit_name + self.site_name + str(self.yr))

        final_gdf = gpd.GeoDataFrame(geometry=final_poly)
        final_gdf = final_gdf.assign(altitude=altitude)  # keep the altitude info with the geodataframe
        return final_gdf

    def run_once(self, **kwargs):
        """
        Runs the model once and closes the object. Only runs if self.project_dir is a temporary directory; else, see
        CreateActiveSpace.run_model()

        Parameters
        ----------
        kwargs (dict): arguments to feed run_model()

        Returns
        -------
        final_poly (geodataframe): multipolygon geodataframe of the final active space shape
        OR None, if triangulation or contouring fails. This can occur if the source is too quiet (or loud).
        """

        self.create_loc()
        ret = self.run_model(**kwargs)
        self.close()

        return ret

    def plot(self, crs=None, ax=None, **kwargs):
        """
        Plot the current site location on an existing axis, or create one to return.

        Parameters
        ----------
        crs (str): optional, coordinate system of the point to return. Default: self.in_crs
        ax (Axes obj): optional, axis to plot on.
        **kwargs (dict): optional line spec plotting for pyplot

        Returns
        -------
        ax (Axes obj): Matplotlib axis object
        """

        pt = gpd.GeoDataFrame(geometry=[Point(self.coord)], crs=self.in_crs)
        if crs is not None:
            pt.to_crs(crs, inplace=True)

        if ax is None:
            _, ax = plt.subplots()

        if kwargs is None:  # default line specifications
            kwargs = {'marker': '*', 'color': 'black', 'markersize': 20, 'legend': 'Site Location'}

        pt.plot(ax=ax, **kwargs)
        return ax

    def close(self):
        """
        Deletes the project directory associated with the CreateActiveSpace object if it is located in the temporary
        directory.

        Returns None
        """

        try:
            self.tmp.cleanup()
            self.tmp = None  # throw an error the next time
        except AttributeError:
            pass  # already closed
