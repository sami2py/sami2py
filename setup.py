#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK
# Full license can be found in License.md
# -----------------------------------------------------------------------------

from __future__ import print_function
import sys
from os import path, mkdir
from setuptools import setup, find_packages
import subprocess

HOME = os.path.expanduser('~')


# Define a read function for using README for long_description
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# generate path for fortran model files
here = os.path.abspath(os.path.dirname(__file__))
fortran_path = os.path.join(here, 'sami2py', 'fortran')
test_data_path = os.path.join(here, 'sami2py', 'tests', 'test_data')
file_path = os.path.join(HOME, '.sami2py')

# %% build


if not path.isfile(path.join(fortran_path, 'sami2py.x')):
    try:  # py27 does not have shutil.which()
        cmd = ['gfortran', '-fno-range-check', '-fno-automatic', '-ffixed-line-length-none',
               '-o', 'sami2py.x']
        src = ['nrlmsise00_modified.f', 'grid-1.00.f', 'sami2py-1.00.f', 'hwm93.f', 'hwm07e_modified.f90',
               'apexcord.f90', 'hwm14.f90']
        subprocess.call(cmd + src, cwd=fortran_path)
    except OSError:
        pass

if not path.isfile(path.join(fortran_path, 'sami2py.x')):
    print('\nYou will need to compile the fortran files.  Try\n'
          '$  make -C sami2py/fortran compile\n', file=sys.stderr)

if not os.path.isdir(file_path):
    os.mkdir(file_path)
    print('Created {} directory to store settings.'.format(file_path))

with open(os.path.join(file_path, 'fortran_path.txt'), 'w+') as f:
    f.write(fortran_path)
with open(os.path.join(file_path, 'test_data_path.txt'), 'w+') as f:
    f.write(test_data_path)

# Define a test suite

# def sami2_test_suite():
#     import unittest
#
#     test_loader = unittest.TestLoader()
#     test_path = path.join(path.dirname(__file__), 'sami2py/tests')
#     test_suite = test_loader.discover(test_path, pattern='test_*.py')
#     return test_suite

# Run setup

setup(test_suite='setup.sami2py_test_suite')
