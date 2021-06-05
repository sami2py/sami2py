# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & AGB
# Full license can be found in License.md
# -----------------------------------------------------------------------------
""" Tests the utilities functions
"""

from __future__ import (print_function)
import os
import numpy as np
import pytest
import sami2py


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
        assert out_path == os.path.join(sami2py.test_data_dir,
                                        'test', 'lon000', '2012_277')

    def test_generate_path_w_blank_archive_dir(self):
        """Tests generation of a path without archive_dir set"""
        with pytest.raises(NameError):
            sami2py.archive_dir = ''
            sami2py.utils.generate_path(tag='test', lon=0, year=2012, day=277)

    def test_generate_path_w_numeric_tag(self):
        """Tests generation of a path with a numeric tag"""
        with pytest.raises(TypeError):
            sami2py.utils.generate_path(tag=7, lon=0, year=2012, day=277)

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

        with open(sami2py._archive_path, 'r') as archive_file:
            archive_dir = archive_file.readline()
        assert archive_dir == test_data_dir

        if os.path.isdir(tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=tmp_archive_dir)
        else:
            with open(sami2py._archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

    def test_set_archive_dir_exception(self):
        """if the provided path is invalid a value error should be produced"""
        with pytest.raises(ValueError):
            sami2py.utils.set_archive_dir('dummy_invalid_path')


class TestGetUnformattedData():
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
        glat = np.loadtxt(os.path.join(self.model_pathF, 'glatf.dat'))
        assert ret_data.size == glat.size

    def test_get_with_reshape_true(self):
        """Test a successful get of unformatted data with the reshape flag
        set to True
        """
        dim0 = 98 * 101 * 7 + 2  # nf*nz*ni + 2
        dim1 = 6             # nt
        udata = sami2py.utils.get_unformatted_data(self.model_pathU, 'deni',
                                                   dim=(dim0, dim1),
                                                   reshape=True)
        fdata = np.loadtxt(os.path.join(self.model_pathF, 'denif.dat'))
        # unformatted test data has 6 time steps, formatted has 2
        assert udata.size == 3 * fdata.size

    def test_reshape_exception(self):
        """Reshape should raise an error if invalid dimensions are provided"""
        with pytest.raises(ValueError):
            dim0 = 2
            dim1 = 2
            sami2py.utils.get_unformatted_data(self.model_pathU, 'deni',
                                               dim=(dim0, dim1),
                                               reshape=True)

    def file_open_error(self):
        """File open should raise an error if invalid file path provided"""
        with pytest.raises(IOError):
            sami2py.utils.get_unformatted_data(self.model_pathU, 'glat')


class TestFourierFunction():
    """Test basic functionality of the return_fourier function"""
    def setup(self):
        """Setup the x and coefficient variables with general x values and
        coefficients"""
        self.x = np.array([0.11, 0.36, 0.61, 0.86, 1.12, 1.37, 1.62, 1.88,
                           2.13, 2.38, 2.64, 2.89, 3.14, 3.4, 3.65, 3.9,
                           4.16, 4.41, 4.66, 4.92, 5.17, 5.42, 5.68, 5.93,
                           6.18, 6.44, 6.69, 6.94, 7.2, 7.45, 7.7, 7.95,
                           8.2, 8.46, 8.71, 8.96, 9.21, 9.46, 9.72, 9.97,
                          10.23, 10.49, 10.74, 10.99, 11.24, 11.49, 11.74,
                          11.99, 12.24, 12.49, 12.74, 13., 13.26, 13.51, 13.76,
                          14.02, 14.27, 14.53, 14.79, 15.04, 15.29, 15.54,
                          15.79, 16.04, 16.29, 16.54, 16.79, 17.04, 17.29,
                          17.54, 17.8, 18.05, 18.31, 18.56, 18.81, 19.07,
                          19.32, 19.57, 19.83, 20.08, 20.33, 20.59, 20.84,
                          21.09, 21.35, 21.6, 21.85, 22.11, 22.36, 22.61,
                          22.86, 23.12, 23.37, 23.624, 23.87])

        self.coeffs = np.zeros((10, 2))

    def teardown(self):
        """Resets both parameters"""
        del self.x, self.coeffs

    def test_cos(self):
        """Test generation of a simple cosine
        """
        self.coeffs[0, 0] = 1.0

        y = sami2py.utils.return_fourier(self.x, self.coeffs)
        target = np.cos(np.pi * self.x / 12.)
        assert (y == target).all()

    def test_sin(self):
        """Test generation of a simple sine
        """
        self.coeffs[0, 1] = 1.0

        y = sami2py.utils.return_fourier(self.x, self.coeffs)
        target = np.sin(np.pi * self.x / 12.)
        assert (y == target).all()


class TestFourierFit():
    """Test basic functionality of the fourier fitting function"""
    def setup(self):
        """Setup the sin wave to be fit by the fitting routine"""
        self.lt = np.linspace(0, 24, 49)
        self.coeffs = np.zeros((10, 2))
        self.coeffs[0, 0] = 1
        self.v = sami2py.utils.return_fourier(self.lt, self.coeffs)

    def teardown(self):
        """Delete the setup parameters"""
        del self.lt, self.coeffs, self.v

    def test_fit(self):
        """Test the goodness of fit"""
        v0, fit_coefs, cov = sami2py.utils.fourier_fit(self.lt, self.v, 10)
        max_diff = np.max(np.abs(self.coeffs.flatten() - fit_coefs.flatten()))
        assert max_diff < .0000001
        assert v0 < .00000001

    def test_warning(self):
        """Test that the warning is generated properly"""
        nan_drifts = np.array([np.nan])
        with pytest.warns(Warning):
            v0, fit_coefs, cov = sami2py.utils.fourier_fit(self.lt, nan_drifts,
                                                           10)
            assert v0 == 0
            assert (fit_coefs == np.zeros((10, 2))).all()
            assert (cov == np.zeros((10, 2))).all()
