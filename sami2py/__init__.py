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

__version__ = str('0.1a1')

# Imports
#---------------------------------------------------------------------

try:
    from sami2py import (run_model, common, model)
    from sami2py.run_model import (run_model)
    #from sami2py.common import (generate_path)
    from sami2py.model import (model)
except ImportError as e:
    logging.exception('problem importing sami2py: ' + str(e))
