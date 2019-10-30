#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""
sami2py - sami2py is another model of the ionosphere python style
=================================================================

Sami2py is a python module that runs the SAMI2 model, as well as archives,
loads and plots the resulting modeled values. SAMI2 is a model developed
by the Naval Research Laboratory to simulate the motions of plasma in a
2D ionospheric environment along a dipole magnetic field [Huba et al, 2000].
SAMI2 solves for the chemical and dynamical evolution of seven ion species
in this environment (H+, He+, N+, O+, N2+, NO+, and O2+).

"""
from __future__ import print_function
import logging
import os

__version__ = str('0.2.0-dev')

# get home directory
env_dir = os.path.expanduser('~')
# set sami2py directory path in home directory
sami2py_dir = os.path.join(env_dir, '.sami2py')
# make sure a sami2py directory for model output exists
if not os.path.isdir(sami2py_dir):
    # create directory
    os.mkdir(sami2py_dir)
    print('Created {} directory to store settings.'.format(sami2py_dir))


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
    print('Run sami2py.utils.set_archive_dir to set the path to'
          ' top-level directory for model outputs.')

# load fortran directory
with open(os.path.join(sami2py_dir, 'fortran_path.txt'), 'r') as f:
    fortran_dir = f.readline()
# load test_data directory
with open(os.path.join(sami2py_dir, 'test_data_path.txt'), 'r') as f:
    test_data_dir = f.readline()


# import main functions
try:
    from sami2py import _core, _core_class, utils  # noqa: F401
    from sami2py._core import run_model  # noqa: F401
    from sami2py._core_class import Model  # noqa: F401
except ImportError as errstr:
    logging.exception('problem importing sami2py: ' + str(errstr))
