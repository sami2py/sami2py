#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""sami2py - sami2py is another model of the ionosphere python style.

Sami2py is a python module that runs the SAMI2 model, as well as archives,
loads and plots the resulting modeled values. SAMI2 is a model developed
by the Naval Research Laboratory to simulate the motions of plasma in a
2D ionospheric environment along a dipole magnetic field [Huba et al, 2000].
SAMI2 solves for the chemical and dynamical evolution of seven ion species
in this environment (H+, He+, N+, O+, N2+, NO+, and O2+).

"""
import logging
import os
import sys

# set the version
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'version.txt')) as version_file:
    __version__ = version_file.read().strip()
del here, version_file

# get home directory
home_dir = os.path.expanduser('~')
# get virtual environment directory
env_name = os.path.split(sys.prefix)[-1]
# set sami2py directory path in home directory with environment subdirectory
sami2py_dir = os.path.join(home_dir, '.sami2py', env_name)
# make sure a sami2py directory for model output exists
if not os.path.isdir(sami2py_dir):
    # create directory
    os.makedirs(sami2py_dir)
    print('Created {} directory to store settings.'.format(sami2py_dir))

_archive_path = os.path.join(sami2py_dir, 'archive_path.txt')
if os.path.isfile(_archive_path):
    # load up stored data path
    with open(_archive_path, 'r') as fin:
        archive_dir = fin.readline()
        del fin
else:
    # create file
    with open(_archive_path, 'w+') as fout:
        fout.write('')
        del fout
    archive_dir = ''
    print('Run sami2py.utils.set_archive_dir to set the path to'
          ' top-level directory for model outputs.')

# flag, True if on readthedocs
on_rtd = os.environ.get('READTHEDOCS') == 'True'

if not on_rtd:
    # load fortran directory
    with open(os.path.join(sami2py_dir, 'fortran_path.txt'), 'r') as fin:
        fortran_dir = fin.readline()
        del fin
    # load test_data directory
    with open(os.path.join(sami2py_dir, 'test_data_path.txt'), 'r') as fin:
        test_data_dir = fin.readline()
        del fin

del home_dir, env_name, on_rtd

# import main functions
try:
    from sami2py import _core, _core_class, utils  # noqa: F401
    from sami2py._core import run_model  # noqa: F401
    from sami2py._core_class import Model  # noqa: F401
except ImportError as errstr:
    logging.exception('problem importing sami2py: ' + str(errstr))
