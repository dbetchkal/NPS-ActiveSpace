# -*- coding: utf-8 -*-
"""
Runs evaluate_truthing.py with the subprocess module.

Created 2021-08-13 16:20

@author: Kirby Heck
"""

import multiprocessing
import sys

# need the active space scripts as well
sys.path.append(r"C:\Users\kheck\PythonScripts\NMSIM-Python")
from legacy_code.active_space_utils import *


def main():
    multi = False

    u = 'DENA'
    yr = 2021  # this is not true for all the sites, but xyz_UTM() will fix the year for us.

    # find all project directories
    paths = glob.glob(r"T:\ResMgmt\WAGS\Sound\Users\Kirby_Heck\NMSIM_ProjectDirectories\*")
    sites = [os.path.basename(f)[4:8] for f in paths]

    sites = ['UWBT']

    if multi:
        with multiprocessing.Pool(int(os.cpu_count()*0.8)) as p:
            p.map(run, sites)
    else:
        for s in sites:
            run(s)

    print('====== FINISHED ALL SITES =======')


def run(s):
    """
    Runs eval_truthing.py with site `s`
    """
    print("========== run_eval_truthing.py: running " + u + s + str(yr) + " =============")
    subprocess.call(['python', 'evaluate_truthing.py', u, s, str(yr)])


if __name__ == '__main__':
    main()
