#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""
sami2py
-----------

Functions
---------------------------------------------------------------------------
run_model(year, day, lat=0, lon=0, alt=300,
          f107=120, f107a=120, ap=0,
          rmin=100, rmax=2000, gams=3, gamp=3, altmin=85.,
          dthr=0.25, hrinit=0., hrpr=24., hrmax=48.,
          dt0=30., maxstep=100000000, denmin=1.e-6,
          nion1=1, nion2=7, mmass=48, h_scale=1, o_scale=1,
          no_scale=1, o2_scale=1, he_scale=1, n2_scale=1, n_scale=1,
          Tinf_scale=1, Tn_scale=1., euv_scale=1,
          wind_scale=1, hwm_model=14,
          fejer=True, ExB_drifts=np.zeros((10,2)), ve01=0., exb_scale=1,
          alt_crit=150., cqe=7.e-14,
          tag='test', clean=False, test=False)

    Initializes a run of the SAMI2 model and archives the data.
---------------------------------------------------------------------------

Classes
---------------------------------------------------------------------------
model
    Loads, reshapes, and holds SAMI2 output for a given model run
    specified by the user.
---------------------------------------------------------------------------
"""
import logging
import sys
import os

__version__ = str('0.1.2')

# get home directory
env_dir = sys.prefix
# set sami2py directory path in home directory
sami2py_dir = os.path.join(env_dir, '.sami2py')
# make sure a sami2py directory for model output exists
if not os.path.isdir(sami2py_dir):
    # create directory
    os.mkdir(sami2py_dir)
    print(''.join(('Created .sami2py directory in ' + env_dir + ' to',
                   'store settings.')))


archive_path = os.path.join(sami2py_dir, 'archive_path.txt')
if os.path.isfile(archive_path):
    # load up stored data path
    with open(archive_path, 'r') as f:
        archive_dir = f.readline()
else:
    # create file
    with open(archive_path, 'w+') as f:
        f.write('')
    archive_dir = ''
    print(''.join(('Run sami2py.utils.set_archive_dir to set the path to',
                   ' top-level directory for model outputs.')))

# load fortran directory
with open(os.path.join(sami2py_dir, 'fortran_path.txt'), 'r') as f:
    fortran_dir = f.readline()
# load test_data directory
with open(os.path.join(sami2py_dir, 'test_data_path.txt'), 'r') as f:
    test_data_dir = f.readline()


# import main functions
try:
    from sami2py import _core, _core_class, utils
    from sami2py._core import run_model
    from sami2py._core_class import Model
except ImportError as errstr:
    logging.exception('problem importing sami2py: ' + str(errstr))
