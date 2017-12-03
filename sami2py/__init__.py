#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
#-----------------------------------------------------------------------------
"""
sami2py
-----------


Functions
---------------------------------------------------------------------------
---------------------------------------------------------------------------

Classes
---------------------------------------------------------------------------

---------------------------------------------------------------------------

Modules
---------------------------------------------------------------------------
---------------------------------------------------------------------------
"""
import logging
import os

__version__ = str('0.1a1')

# Imports
#---------------------------------------------------------------------

# get home directory
home_dir = os.path.expanduser('~')
# set sami2py directory path in home directory
sami2py_dir = os.path.join(home_dir, '.sami2py')
# make sure a sami2py directory exists
if not os.path.isdir(sami2py_dir):
    # create directory
    os.mkdir(sami2py_dir)
    # create file
    with open(os.path.join(sami2py_dir, 'model_path.txt'),'w') as f:
        f.write('')
    print('Created .sami2py directory in user home directory to store settings.')
    model_dir=''
else:
    # load up stored data path
    with open(os.path.join(sami2py_dir, 'model_path.txt'),'r') as f:
        model_dir = f.readline()

if model_dir == '':
    print(''.join(('Run sami2py.utils.set_model_dir to set the path',
        ' to top-level directory that will/does contain science data.')))

# import main functions
try:
    from sami2py import (run_model, utils, model)
    from sami2py.run_model import (run_model)
    from sami2py.model import (model)
except ImportError as e:
    logging.exception('problem importing sami2py: ' + str(e))
