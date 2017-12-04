#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
#-----------------------------------------------------------------------------
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

Modules
---------------------------------------------------------------------------

---------------------------------------------------------------------------
"""
import logging
import os

__version__ = str('0.1a1')

# Imports
#---------------------------------------------------------------------

# get home directory
home_dir = os.path.expanduser('~')
# set sami2py directory path in home directory
sami2py_dir = os.path.join(home_dir, '.sami2py')
# make sure a sami2py directory exists
if not os.path.isdir(sami2py_dir):
    # create directory
    os.mkdir(sami2py_dir)
    # create file
    with open(os.path.join(sami2py_dir, 'model_path.txt'),'w') as f:
        f.write('')
    print('Created .sami2py directory in user home directory to store settings.')
    model_dir=''
else:
    # load up stored data path
    with open(os.path.join(sami2py_dir, 'model_path.txt'),'r') as f:
        model_dir = f.readline()

if model_dir == '':
    print(''.join(('Run sami2py.utils.set_model_dir to set the path',
        ' to top-level directory that will/does contain model outputs.')))

# import main functions
try:
    from sami2py import (run_model, utils, model)
    from sami2py.run_model import (run_model)
    from sami2py.model import (model)
except ImportError as e:
    logging.exception('problem importing sami2py: ' + str(e))
