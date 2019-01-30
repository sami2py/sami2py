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

    Returns
    -------
    path : (string)
        Complete path pointing to model archive for a given run
    """
    from os import path

    if test:
        from sami2py import test_data_dir
        top_directory = test_data_dir
    else:
        from sami2py import archive_dir
        top_directory = archive_dir

    if top_directory:
        archive_path = path.join(top_directory, tag,
                                 ('lon%03d/%4d_%03d/' % (lon, year, day)))
    else:
        raise NameError(''.join(('Archive Directory Not Specified: ',
                                 'Run sami2py.utils.set_archive_dir')))

    return archive_path


def set_archive_dir(path=None, store=None):
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

    if store is None:
        store = True

    if os.path.isdir(path):
        if store:
            with open(os.path.join(os.path.expanduser('~'),
                                   '.sami2py', 'archive_path.txt'),
                      'w') as archive_file:
                archive_file.write(path)
        sami2py.archive_dir = path
    else:
        raise ValueError('Path does not lead to a valid directory.')
