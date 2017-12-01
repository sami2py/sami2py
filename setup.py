#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Copyright (C) 2017, JK
# Full license can be found in License.md
#-----------------------------------------------------------------------------

from os import path
from setuptools import setup, find_packages

# Define a read function for using README for long_description

def read(fname):
    return open(path.join(path.dirname(__file__), fname)).read()

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
      version='0.1a1',
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
          "Programming Language :: Python :: 3.3",
          "Programming Language :: Python :: 3.4",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Operating System :: MacOS :: MacOS X",
      ],
      include_package_data=True,
      zip_safe=False,
      test_suite='setup.sami2py_test_suite',
)
