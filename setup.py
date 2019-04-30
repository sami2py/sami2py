#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK
# Full license can be found in License.md
# -----------------------------------------------------------------------------

from os import path, mkdir
from setuptools import setup, find_packages


# Define a read function for using README for long_description
def read(fname):
    return open(path.join(path.dirname(__file__), fname)).read()


# generate path for fortran model files
here = path.abspath(path.dirname(__file__))
fortran_path = path.join(here, 'sami2py', 'fortran')
test_data_path = path.join(here, 'sami2py', 'tests', 'test_data')
file_path = path.join(path.expanduser('~'), '.sami2py')

if not path.isfile(fortran_path + '/sami2py.x'):
    print('\n'.join(['\nYou will need to compile the fortran files.  Try',
                     '$  cd sami2py/fortran',
                     '$  make compile\n']))

if not path.isdir(file_path):
    mkdir(file_path)

with open(path.join(file_path, 'fortran_path.txt'), 'w+') as f:
    f.write(fortran_path)
with open(path.join(file_path, 'test_data_path.txt'), 'w+') as f:
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

setup(name='sami2py',
      version='0.1a2',
      url='gitlab.com/jklenzing/sami2py',
      author='Jeff Klenzing',
      author_email='jeffrey.klenzing@nasa.gov',
      description='Generate, read, and plot sami2 model runs',
      long_description=read('README.md'),
      packages=find_packages(),
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Topic :: Scientific/Engineering :: Physics",
          "Intended Audience :: Science/Research",
          "License :: BSD",
          "Natural Language :: English",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.4",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Operating System :: MacOS :: MacOS X",
      ],
      include_package_data=True,
      zip_safe=False,
      test_suite='setup.sami2py_test_suite',
      )
