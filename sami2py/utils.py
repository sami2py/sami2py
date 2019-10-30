#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
""" Wrapper for running sami2 model

Functions
-------------------------------------------------------------------------------
generate_path(tag, lon, year, day)
    Generates path to archived model runs based on input paramters.

set_archive_dir(path=None, store=None)
    Allows user to specify the location where the model outputs will be stored

get_unformatted_data(dat_dir, var_name, nz, nf, ni, nt, reshape=False)
    routine to interpret unformatted binary files created by the SAMI2 model
-------------------------------------------------------------------------------

Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 1 Dec 2017, Goddard Space Flight Center (GSFC)
-------------------------------------------------------------------------------
"""

import os


def generate_path(tag, lon, year, day, test=False):
    """Creates a path based on run tag, date, and longitude

    Parameters
    ----------
    tag : (string)
        specifies name of model run
    lon : (int)
        longitude of model run
    year : (int)
        year of model run
    day : (int)
        day of year of model run
    test : (bool)
        If True, use directory for test data.  If False, use archive_dir
        (default = False)

    Returns
    -------
    archive_path : (string)
        Complete path pointing to model archive for a given run
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
        str_fmt = 'lon{lon:03d}/{year:4d}_{day:03d}/'
        archive_path = os.path.join(top_directory, tag,
                                    (str_fmt.format(lon=lon,
                                                    year=year,
                                                    day=day)))
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


def get_unformatted_data(dat_dir, var_name, reshape=False, dim0=0, dim1=0):
    """Routine to interpret unformatted binary files created by the SAMI2 model

    Parameters
    -----------
    data_dir : (str)
        directory where the SAMI2 data is stored
    var_name : (str)
        name of unformatted data variable to be loaded
    nz : (int)
        number of mesh points along the geomagnetic field line
    nf : (int)
        number of mesh points transverse to the geomagnetic field line i.e.
        number of field lines
    ni : (int)
        number of ion species
    nt : (int)
        number of time steps
    reshape : (bool)
        if true the data is reshaped by the mesh geometry

    Returns
    -----------
    float_data : (numpy.ndarray)
        unformatted data organized into a numpy array for handling in python
    """
    import numpy as np

    binary_file = open(os.path.join(dat_dir, var_name + 'u.dat'), 'rb')
    float_data = np.fromfile(binary_file, dtype='float32')
    binary_file.close()

    if reshape:
        float_data = np.reshape(float_data, (dim0, dim1), order='F')
        return float_data[1:-1, :]
    else:
        return float_data[1:-1]
