#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & AGB
# Full license can be found in License.md
# -----------------------------------------------------------------------------
""" Tests the utilities functions
"""

from __future__ import (print_function)
import os
import numpy as np
import sami2py
import pytest


class TestGeneratePath():
    """Test basic functionality of the generate_path function
    """
    def setup(self):
        """Runs before every method to create a clean testing setup."""
        sami2py.archive_dir = 'test'

    def test_successful_path_generation(self):
        """Tests generation of a path that is successful"""
        out_path = sami2py.utils.generate_path(tag='test', lon=0, year=2012,
                                               day=277, test=True)
        assert out_path == sami2py.test_data_dir + '/test/lon000/2012_277/'

    def test_generate_path_w_blank_archive_dir(self):
        """Tests generation of a path without archive_dir set"""
        with pytest.raises(NameError):
            sami2py.archive_dir = ''
            sami2py.utils.generate_path(tag='test', lon=0, year=2012, day=277)

    def test_generate_path_w_numeric_tag(self):
        """Tests generation of a path with a numeric tag"""
        with pytest.raises(TypeError):
            sami2py.utils.generate_path(tag=7, lon=0, year=2012, day=277)

    def test_generate_path_w_nonnumeric_lon(self):
        """Tests generation of a path with a nonnumeric longitude"""
        with pytest.raises(ValueError):
            sami2py.utils.generate_path(tag='test', lon='0', year=2012,
                                        day=277)

    def test_generate_path_w_nonnumeric_year(self):
        """Tests generation of a path with a nonnumeric year"""
        with pytest.raises(ValueError):
            sami2py.utils.generate_path(tag='test', lon=0, year='2012',
                                        day=277)

    def test_generate_path_w_nonnumeric_day(self):
        """Tests generation of a path with a nonnumeric day"""
        with pytest.raises(ValueError):
            sami2py.utils.generate_path(tag='test', lon=0, year=2012,
                                        day='277')


class TestArchiveDir():
    """Test basic functionality of the set_archive_dir function"""
    def test_set_archive_dir(self):
        """Test that set_archive_dir has set and stored the archive directory

           To leave the archive directory unchanged it must be gathered and
           reset after the test is complete
        """
        tmp_archive_dir = sami2py.archive_dir

        from sami2py import test_data_dir
        sami2py.utils.set_archive_dir(path=test_data_dir)

        with open(sami2py.archive_path, 'r') as archive_file:
            archive_dir = archive_file.readline()
        assert archive_dir == test_data_dir

        if os.path.isdir(tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=tmp_archive_dir)
        else:
            with open(sami2py.archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

    def test_set_archive_dir_exception(self):
        """if the provided path is invalid a value error should be produced"""
        with pytest.raises(ValueError):
            sami2py.utils.set_archive_dir('dummy_invalid_path')


class testGetUnformattedData():
    """Test basic functionality of the get_unformatted_data function"""
    def setup(self):
        """setup the model_path variable for accessing unformatted data"""
        self.model_pathF = sami2py.utils.generate_path('test', 256, 1999, 256,
                                                       test=True)
        self.model_pathU = sami2py.utils.generate_path('test', 256, 1999, 257,
                                                       test=True)

    def test_successful_get(self):
        """Test a successful get of unformatted data"""
        ret_data = sami2py.utils.get_unformatted_data(self.model_pathU, 'glat')
        glat = np.loadtxt(self.model_pathF + 'glatf.dat')
        assert ret_data.size == glat.size

    def test_get_with_reshape_true(self):
        """Test a successful get of unformatted data with the reshape flag
        set to True
        """
        dim0 = 98*101*7 + 2  # nf*nz*ni + 2
        dim1 = 2             # nt
        ret_data = sami2py.utils.get_unformatted_data(self.model_pathU, 'deni',
                                                      dim0=dim0, dim1=dim1,
                                                      reshape=True)
        glat = np.loadtxt(self.model_pathF + 'denif.dat')
        assert ret_data.size == glat.size

    def test_reshape_exception(self):
        """Reshape should raise an error if invalid dimensions are provided"""
        with pytest.raises(ValueError):
            dim0 = 2
            dim1 = 2
            sami2py.utils.get_unformatted_data(self.model_pathU, 'deni',
                                               dim0=dim0, dim1=dim1,
                                               reshape=True)

    def file_open_error(self):
        """File open should raise an error if invalid file path provided"""
        with pytest.raises(IOError):
            sami2py.utils.get_unformatted_data(self.model_pathU, 'glat')
