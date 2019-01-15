'''
    Unit tests for run_model.py
'''
import sami2py
import os
import numpy as np
from nose.tools import raises


class test_basic_model_run():

    def setup(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat', 'sami2low-1.00.namelist']
        filelist_ufmt = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                         'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                         'time.dat', 'sami2low-1.00.namelist']
        for filename in filelist_fmt:
            open(filename, 'w').close()
        for filename in filelist_ufmt:
            open(filename, 'w').close()

    def teardown(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat', 'sami2low-1.00.namelist']
        filelist_ufmt = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                         'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                         'time.dat', 'sami2low-1.00.namelist']
        for filename in filelist_fmt:
            os.remove(filename)
        for filename in filelist_ufmt:
            os.remove(filename)

    def test_run_model_namelist(self):
        sami2py.run_model(year=2012, day=211)
        namelist_file = open(fortran_dir+'sami2low-1.00.namelist', 'rt')
        ref_namelist = open('tests/test_data/reference_sami2low-1.00.namelist',
                            'rt')
        assert namelist_file == ref_namelist
        assert os.stat('tests/test_data/glonf.dat')

    def test_run_model_dat_files(self):
        sami2py.run_model(year=2012, day=211)
        assert os.stat('tests/test_data/glonf.dat')
