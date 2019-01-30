<<<<<<< HEAD
'''
   Unit tests for utils.py
'''

import os
import sami2py
from nose.tools import raises
from sami2py.utils import generate_path, set_archive_dir


class TestGeneratePath():
    def test_generate_path(self):
        '''test the basic functionality of generating a path
        '''
        from sami2py import test_data_dir
        out = generate_path(tag='test', lon=0, year=1900, day=0, test=True)
        assert out == test_data_dir + 'test/lon000/1900_000/'

    @raises(NameError)
    def test_generate_path_exception(self):
        '''if no archive directory has been specified a name error should be
           produced
        '''
        tmp_archive_dir = sami2py.archive_dir
        sami2py.archive_dir = ''
        out = generate_path(tag='test', lon=0, year=0, day=0)
        sami2py.archive_dir = tmp_archive_dir


class TestArchiveDir():
    def test_set_archive_dir(self):
        '''test that set_archive_dir has set and stored the archive directory
        '''
        from sami2py import test_data_dir
        tmp_archive_dir = sami2py.archive_dir
        set_archive_dir(path=test_data_dir)
        home_dir = os.path.expanduser('~')
        sami2py_dir = os.path.join(home_dir, '.sami2py')
        archive_path = os.path.join(sami2py_dir, 'archive_path.txt')
        with open(archive_path, 'r') as f:
            archive_dir = f.readline()
        assert archive_dir == test_data_dir
        # return the archive dir to its previous value
        set_archive_dir(path=tmp_archive_dir)

    @raises(ValueError)
    def test_set_archive_dir_exception(self):
        '''if the provided path is invalid a value error should be produced
        '''
        set_archive_dir('dummy_invalid_path')
=======
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & AGB
# Full license can be found in License.md
# -----------------------------------------------------------------------------
""" Tests the utilities functions
"""

from __future__ import (print_function)
import sami2py
from nose.tools import assert_raises, raises
import nose.tools
import numpy as np


class TestUtils():

    def setUp(self):
        """Runs before every method to create a clean testing setup."""
        sami2py.archive_dir = 'test'

    @raises(NameError)
    def test_generate_path_w_blank_archive_dir(self):
        """Tests generation of a path without archive_dir set"""
        sami2py.archive_dir = ''
        sami2py.utils.generate_path(tag='test', lon=0, year=2012, day=277)

    @raises(TypeError)
    def test_generate_path_w_numeric_tag(self):
        """Tests generation of a path with a numeric tag"""

        sami2py.utils.generate_path(tag=7, lon=0, year=2012, day=277)

    @raises(TypeError)
    def test_generate_path_w_nonnumeric_lon(self):
        """Tests generation of a path with a nonnumeric longitude"""

        sami2py.utils.generate_path(tag='test', lon='0', year=2012, day=277)

    @raises(TypeError)
    def test_generate_path_w_nonnumeric_year(self):
        """Tests generation of a path with a nonnumeric year"""

        sami2py.utils.generate_path(tag='test', lon=0, year='2012', day=277)

    @raises(TypeError)
    def test_generate_path_w_nonnumeric_day(self):
        """Tests generation of a path with a nonnumeric longitude"""

        sami2py.utils.generate_path(tag='test', lon=0, year=2012, day='277')
>>>>>>> be7a2291ba5f1baa1742f1e135faec2a1718c641
