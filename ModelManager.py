# -*- coding: utf-8 -*-
"""
Manages CreateActiveSpace object and runs independently of any existing directory infrastructure.

Created 2021-08-12 11:29

@author: Kirby Heck
"""

import multiprocessing
from multiprocessing.dummy import Pool as DummyPool
import sys
from sqlalchemy.exc import OperationalError

from active_space_utils import *
from NestablePool import NestablePool

sys.path.append(os.path.join(os.getcwd(), "NMSIM-Python"))
from ActiveSpace import CreateActiveSpace


def timeout_wrapper(func, *args, **kwargs):
    """
    Wraps self.run() [or any other function `func`] that returns None if a run times out.

    Parameters
    ----------
    func (method): method call, static or otherwise
    *args: function arguments for `func`
    **kwargs (dict): keyword arguments for timeout_wrapper(). Currently only uses 'timeout'

    Returns
    -------
    output of `func` if the function completes before `timeout` seconds
    None otherwise
    """
    timeout = kwargs.get('timeout', None)
    p = DummyPool(1)  # create a dummy pool for this thread
    res = p.apply_async(func, args=args)  # call the wrapped function

    try:
        # attempt get a return from the function `func` before timing out
        out = res.get(timeout=timeout)
        return out  # return function results

    except multiprocessing.context.TimeoutError:
        print('Process {}: Timed out after {} seconds, terminated.'.format(os.getpid(), timeout))
        return None


class ModelManager:

    def __init__(self, coords_ls=None, crs=None, n_jobs=1, timeout=45., init_args=None, run_args=None):
        """
        Initializes a ModelManager object that wraps the CreateActiveSpace object.

        Parameters
        ----------
        coords (list): list of lists [id, (x, y)] to run through the model from get_coords_list()
        crs (str): coordinate projection of the given coordinates
        n_jobs (int OR float): optional, number of multiprocessing processes. If n_jobs is between 0 and 1, will use
            that fraction of total computational cores.
        timeout (float): optional, length of time to wait for each Active Space computation before terminating the
            process. Default 45 seconds
        init_args (dict): initialization arguments, see CreateActiveSpace in ActiveSpace.py. Omit `coord=` and `crs=`
            arguments (will be filled in by the first two parameters)
        run_args (dict): model arguments for CreateActiveSpace.run_model()
        """

        if coords_ls is None or crs is None:
            raise ValueError('Missing arguments: coords list or crs.')

        self.coords_ls = coords_ls  # list of [id, coordinate tuples] lists
        self.crs = crs
        self.timeout = timeout

        # ok, let's start fresh
        if n_jobs > 0:
            if n_jobs < 1:  # get fraction of total CPUs (cannot use all cores)
                n_jobs = os.cpu_count() * n_jobs

            if n_jobs > os.cpu_count():
                n_jobs = os.cpu_count()-1

            self.n_jobs = int(n_jobs)
        else:
            raise ValueError('Number of jobs must be positive')

        # default CreateActiveSpace init args; this does NOT include coord= or crs=
        if init_args is None or len(init_args) == 0:
            # if this is none, choose a more realistic ambience source than the default Cessna 206
            sound_source_dir = r'T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_TuningSources'
            sound_sources = get_omni_sources(sound_source_dir, upper=23, lower=23)
            sound_src = sound_sources.full_path.values[0]

            init_args = {'project_dir': 'tmp',
                         'ambience_source': 'Mennitt',
                         'source_path': sound_src}
        self.init_args = init_args

        # fill in default run_model arguments
        if run_args is None or len(run_args) == 0:
            run_args = {'altitude_ft': 'calc',
                        'n_tracks': 1,
                        'out_crs': crs}
        self.run_args = run_args

    def run_all(self, **run_args):
        """
        Runs all CreateActiveSpace models with an average flight altitude based on runtime arguments provided upon
        ModelManager initialization, or new ones provided.

        Parameters
        ----------
        run_args correspond with CreateActiveSpace.run_model() key word arguments

        Returns
        -------
        GeoDataFrame of Active Space results
        """

        print("\n========== ModelManager: Running {} active space locations ==========".format(len(self.coords_ls)))
        print("\tInitiating {} process(es)...\n\n".format(self.n_jobs))

        kwargs = {'timeout': self.timeout}  # timeout limit

        if run_args is not None and len(run_args) > 0:
            self.run_args = run_args

        if self.n_jobs == 1:
            gdf_list = []

            for input_args in self.coords_ls:
                # call the helper run function with the timeout wrapper
                gdf_list.append(timeout_wrapper(self._run, input_args, **kwargs))

        else:  # start multiprocessing
            # create a pool to multiprocess this task
            with NestablePool(processes=self.n_jobs, maxtasksperchild=1) as p:
                # call the timeout_wrapper function to run the CreateActiveSpace model with each coordinate
                res_ls = [p.apply_async(timeout_wrapper, args=(self._run, coord, ), kwds=kwargs)
                          for coord in self.coords_ls]

                # wait for results to come in asynchronously
                gdf_list = [res.get() for res in res_ls]

                p.close()  # not sure what these do, but I'm hopping on the bandwagon
                p.join()

        return pd.concat(gdf_list)  # returns the list as one GeoDataFrame

    def _run(self, input_args):
        """
        Actually runs the CreateActiveSpace model.
        """

        id = input_args[0]
        coord = input_args[1]
        print('\n{}: Running active space computations at {}: '.format(id, coord))

        ret = None
        count = 1
        while ret is None:
            try:
                cas = CreateActiveSpace(coord=coord, crs=self.crs, **self.init_args)
                ret = cas.run_once(**self.run_args)

            except OperationalError:
                # this could quickly get dangerous if there's a permanent/long-term error at play...
                # so stop the attempts after 5
                if count > 5:
                    ret = gpd.GeoDataFrame(index=[id])
                    return ret  # return a blank geodataframe

                count += 1
                warnings.warn("{}: Couldn't connect to the server... "
                              "Recomputing (beginning attempt #{})".format(id, count))

        ret = ret.assign(id=id)
        ret.set_index('id', inplace=True)
        return ret
