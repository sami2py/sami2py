#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
#-----------------------------------------------------------------------------
""" Wrapper for running sami2 model

Functions
-------------------------------------------------------------------------------

generate_path(tag, lon, year, day)
    Generates path to archived model runs based on input paramters.


Classes
-------------------------------------------------------------------------------


Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 1 Dec 2017, Goddard Space Flight Center (GSFC)

References
-------------------------------------------------------------------------------


"""

import platform

def generate_path(tag, lon, year, day):
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

    if 'gs674-jklenmbp' in platform.node():
        basedir = '/Users/jklenzin/data/sami2/'
    else:
        basedir = '/Volumes/drive/models/sami2/'

    path = basedir + tag + ('/lon%03d/%4d_%03d/' % (lon, year, day))
    return path
