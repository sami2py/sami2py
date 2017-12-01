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
import platform
import os
import shutil
import numpy as np


__version__ = str('0.1a1')

# Imports
#---------------------------------------------------------------------

try:
    from sami2py import (models, run_model, _generate_path)
    from sami2py.models import (Sami2Model)
except ImportError as e:
    logging.exception('problem importing sami2py: ' + str(e))
