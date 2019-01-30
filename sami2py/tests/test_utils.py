#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & AGB
# Full license can be found in License.md
# -----------------------------------------------------------------------------
""" Tests the utilities functions
"""

from __future__ import (print_function)
import sami2py
from nose.tools import assert_raises, raises
import nose.tools
import numpy as np


class TestUtils():

    @raises(TypeError)
    def test_generate_path_w_numeric_tag(self):
        """Tests generation of a path with a numeric tag"""

        sami2py.utils.generate_path(tag=7, lon=0, year=2012, day=277)

    @raises(TypeError)
    def test_generate_path_w_nonnumeric_lon(self):
        """Tests generation of a path with a nonnumeric longitude"""

        sami2py.utils.generate_path(tag='test', lon='0', year=2012, day=277)

    @raises(TypeError)
    def test_generate_path_w_nonnumeric_year(self):
        """Tests generation of a path with a nonnumeric year"""

        sami2py.utils.generate_path(tag='test', lon=0, year='2012', day=277)

    @raises(TypeError)
    def test_generate_path_w_nonnumeric_day(self):
        """Tests generation of a path with a nonnumeric longitude"""

        sami2py.utils.generate_path(tag='test', lon=0, year=2012, day='277')
