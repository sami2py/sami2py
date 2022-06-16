#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Wrapper for running sami2 model.

Functions
---------
generate_path(tag, lon, year, day)
    Generates path to archived model runs based on input paramters.

set_archive_dir(path=None, store=None)
    Allows user to specify the location where the model outputs will be stored

return_fourier(x, coeffs)
    Returns Fourier Series up to NumF Coefficients

get_unformatted_data(dat_dir, var_name, nz, nf, ni, nt, reshape=False)
    routine to interpret unformatted binary files created by the SAMI2 model


"""

import numpy as np
import os
import warnings

from scipy.optimize import curve_fit


def generate_path(tag, lon, year, day, test=False):
    """Create a path based on run tag, date, and longitude.

    Parameters
    ----------
    tag : str
        specifies name of model run
    lon : int or float
        longitude of model run
    year : int
        year of model run
    day : int
        day of year of model run
    test : bool
        If True, use directory for test data.  If False, use archive_dir
        (default = False)

    Returns
    -------
    archive_path : str
        Complete path pointing to model archive for a given run

    Note
    ----
    The longitude value will be rounded to an integer for creating the path,
    but the simulation will store and use this as a float so that intersections
    with specific ground-based stations can be performed.

    Examples
    --------
        import sami2py
        sami2py.utils.set_archive_dir(path='path_name_here')
        path = sami2py.utils.generate_path(tag='run_name', lon=0, year=2012,
                                           day=210)
    Will return 'path_name_here/run_name/lon000/2012_210'

    """

    if not isinstance(tag, str):
        raise TypeError

    if test:
        from sami2py import test_data_dir
        top_directory = test_data_dir
    else:
        from sami2py import archive_dir
        top_directory = archive_dir

    # Check if top_directory is empty string, ie, user has not specified
    # a directory through set_archive_dir
    if top_directory:
        str_fmt1 = 'lon{lon:03d}'
        str_fmt2 = '{year:4d}_{day:03d}'
        archive_path = os.path.join(top_directory, tag,
                                    str_fmt1.format(lon=int(lon)),
                                    str_fmt2.format(year=year,
                                                    day=day))
    else:
        raise NameError(''.join(('Archive Directory Not Specified: ',
                                 'Run sami2py.utils.set_archive_dir')))

    return archive_path


def set_archive_dir(path=None, store=True):
    # type: (str, bool) -> None
    """Set the top level directory sami2py uses to look for data and reload.

    Parameters
    ----------
    path : string
        valid path to directory sami2py uses to look for data
    store : bool
        if True, store data directory for future runs

    Examples
    --------
    Should be run upon first installation.  Will prompt users if not run.
        import sami2py
        sami2py.utils.set_archive_dir(path='path_name_here')
    """
    import sami2py

    path = os.path.expanduser(path)
    if os.path.isdir(path):
        if store:
            with open(os.path.join(sami2py.sami2py_dir, 'archive_path.txt'),
                      'w') as archive_file:
                archive_file.write(path)
        sami2py.archive_dir = path
    else:
        raise ValueError('Path does not lead to a valid directory.')


def return_fourier(x, coeffs):
    """Return a Fourier series up to NumF coefficients.

    Parameters
    ----------
    x : 1d ndarray
        solar local time in hours (slt)
    coeffs : array
        10x2 array of fourier coefficients

    Returns
    --------
    y : array
        result of the fourier series

    """

    def cos_a(x, n):
        """Calculate a simple cosine.

        Parameters
        ----------
        x : np.ndarray
            Values of Solar Local Time
        n : int
            Number of oscillations

        Returns
        -------
        y : np.ndarray
            Cosine wave over x.

        """

        return np.cos(n * np.pi * x / 12.0)

    def sin_a(x, n):
        """Calculate a simple sine.

        Parameters
        ----------
        x : np.ndarray
            Values of Solar Local Time
        n : int
            Number of oscillations

        Returns
        -------
        y : np.ndarray
            Sine wave over x.

        """

        return np.sin(n * np.pi * x / 12.0)

    shape = coeffs.shape

    y = 0.0 * x
    for i in range(0, shape[0]):
        y += coeffs[i, 0] * cos_a(x, i + 1) + coeffs[i, 1] * sin_a(x, i + 1)

    return y


def get_unformatted_data(dat_dir, var_name, reshape=False, dim=(0, 0)):
    """Interpret unformatted binary files created by the SAMI2 model.

    Parameters
    -----------
    data_dir : str
        directory where the SAMI2 data is stored
    var_name : str
        name of unformatted data variable to be loaded
    nz : int
        number of mesh points along the geomagnetic field line
    nf : int
        number of mesh points transverse to the geomagnetic field line i.e.
        number of field lines
    ni : int
        number of ion species
    nt : int
        number of time steps
    reshape : bool
        if true the data is reshaped by the mesh geometry

    Returns
    -----------
    float_data : numpy.ndarray
        unformatted data organized into a numpy array for handling in python

    """

    binary_file = open(os.path.join(dat_dir, var_name + 'u.dat'), 'rb')
    float_data = np.fromfile(binary_file, dtype='float32')
    binary_file.close()

    if reshape:
        float_data = np.reshape(float_data, dim, order='F')
        return float_data[1:-1, :]
    else:
        return float_data[1:-1]


def _make_fourier(na, nb):
    """Make a fourier series to use in the curve fits.

    Parameters
    ----------
    na : int
        number of cosine terms/coefficients
    nb : int
        number of sin terms/coefficients

    """

    def fourier(x, *a):
        ret = a[0]
        for deg in range(0, na):
            ret += a[deg + 1] * np.cos((deg + 1) * np.pi * x / 12)
        for deg in range(na, na + nb):
            ret += a[deg + 1] * np.sin((deg - na + 1) * np.pi * x / 12)
        return ret

    return fourier


def fourier_fit(local_times, drifts, num_co):
    """Determine the terms in a fourier fit to data.

    Parameters
    ----------
    local_times : array-like
        xdim for fit; local time values
    drifts : array-like
        ydim for fit; median drift values from data
    num_co : int
        'number of coefficients) how many sin/cosine pairs for the fit

    Returns
    -------
    ve01 : float
        linear offset of the fourier fit
    coefficients : num_co by 2 array like
        coefficients to describe the fourier function that fits the drifts
    covariance : num_co by 2 array like
        covariance of the coefficients

    """

    coefficients = np.zeros((num_co, 2))
    covariance = np.zeros((num_co, 2))
    ind, = np.where(~np.isnan(drifts))
    if ind.size < num_co * 2 + 1:
        warnings.warn('not enough viable drift data, '
                      'returning zero value \"flat fit\"', Warning)
        return 0, coefficients, covariance
    # popt contains the coeficients. First ten are cosines, second ten are sins
    popt, pcov = curve_fit(_make_fourier(num_co, num_co), local_times[ind],
                           drifts[ind], [0.0] * (num_co * 2 + 1))
    # format the coefficients for input ito the SAMI2 model
    # the shape is np.zeroes((10,2))
    ve01 = popt[0]
    for n in range(1, num_co * 2):
        i = (n - 1) % num_co
        j = int((n - 1) / num_co)
        coefficients[i, j] = popt[n]
        covariance[i, j] = pcov[n, n]

    return ve01, coefficients, covariance
