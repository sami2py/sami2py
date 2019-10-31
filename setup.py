#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK
# Full license can be found in License.md
# -----------------------------------------------------------------------------

from __future__ import print_function
import os
import sys
from setuptools import setup
import subprocess

# get home directory
home_dir = os.path.expanduser('~')
# get name of virtual environment
env_name = os.path.split(sys.prefix)[-1]

# generate path for fortran model files
here = os.path.abspath(os.path.dirname(__file__))
fortran_path = os.path.join(here, 'sami2py', 'fortran')
# generate path for test files
test_data_path = os.path.join(here, 'sami2py', 'tests', 'test_data')
# generate path to store test and fortran file paths
file_path = os.path.join(home_dir, '.sami2py', env_name)

# %% build


if not os.path.isfile(os.path.join(fortran_path, 'sami2py.x')):
    try:  # py27 does not have shutil.which()
        cmd = ['gfortran', '-fno-range-check', '-fno-automatic',
               '-ffixed-line-length-none', '-o', 'sami2py.x']
        src = ['nrlmsise00_modified.f', 'grid-1.00.f', 'sami2py-1.00.f',
               'hwm93.f', 'hwm07e_modified.f90', 'apexcord.f90', 'hwm14.f90']
        subprocess.call(cmd + src, cwd=fortran_path)
    except OSError:
        pass

if not os.path.isfile(os.path.join(fortran_path, 'sami2py.x')):
    print('\nYou will need to compile the fortran files.  Try\n'
          '$  make -C sami2py/fortran compile\n', file=sys.stderr)

if not os.path.isdir(file_path):
    os.mkdir(file_path)
    print('Created {} directory to store settings.'.format(file_path))

with open(os.path.join(file_path, 'fortran_path.txt'), 'w+') as f:
    f.write(fortran_path)
with open(os.path.join(file_path, 'test_data_path.txt'), 'w+') as f:
    f.write(test_data_path)

setup()
