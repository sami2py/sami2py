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
-------------------------------------------------------------------------------

Classes
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------

Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 1 Dec 2017, Goddard Space Flight Center (GSFC)
-------------------------------------------------------------------------------

References
-------------------------------------------------------------------------------


"""


def generate_path(tag, lon, year, day, test=False):
    """
    Creates a path based on run tag, date, and longitude

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
    path : (string)
        Complete path pointing to model archive for a given run
    """
    from os import path

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
        archive_path = path.join(top_directory, tag,
                                 (str_fmt.format(lon=lon,
                                                 year=year,
                                                 day=day)))
    else:
        raise NameError(''.join(('Archive Directory Not Specified: ',
                                 'Run sami2py.utils.set_archive_dir')))

    return archive_path


def set_archive_dir(path=None, store=True):
    """
    Set the top level directory pysat uses to look for data and reload.

    Parameters
    ----------
    path : string
        valid path to directory pysat uses to look for data
    store : bool
        if True, store data directory for future runs
    """
    import os
    import sami2py

    if os.path.isdir(path):
        if store:
            with open(os.path.join(os.path.expanduser('~'),
                                   '.sami2py', 'archive_path.txt'),
                      'w') as archive_file:
                archive_file.write(path)
        sami2py.archive_dir = path
    else:
        raise ValueError('Path does not lead to a valid directory.')


def load_msis_scalars(year, case='JTE2010'):
    """
    Loads a set of predetermined msis scalars based on previous studies

    Parameters
    ----------
    year : (int)
        The year of interest
    case : (str)
        A label for he case study of interest
        (default = 'JTE2010')

    """

    if case is 'JTE2010':
        snn = [0.9697, 0.9697, 0.9697, 0.9697, 0.7979, 0.9697, 0.9697, 0.9697]
        Tinf_scl = 0.9603
    else:
        snn = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        Tinf_scl = 1.0
