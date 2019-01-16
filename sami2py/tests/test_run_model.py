'''
    Unit tests for run_model.py
'''
import sami2py
import os
import filecmp
import numpy as np
from nose.tools import raises
from sami2py import fortran_dir, test_data_dir
from sami2py.utils import generate_path


class test_basic_model_run():

    def setup(self):
        self.path = generate_path('test', 0, 2012, 211, True)
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat']
        for filename in filelist_fmt:
            open(fortran_dir+'/'+filename, 'w').close()

    def teardown(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat']
        for filename in filelist_fmt:
            os.remove(self.path+filename)
            os.remove(fortran_dir+'/'+filename)

    def test_run_model_namelist(self):
        sami2py.run_model(year=2012, day=211, test=True)
        namelist_file = self.path+'sami2low-1.00.namelist'
        ref_namelist = test_data_dir+'/reference_sami2low-1.00.namelist'
        assert filecmp.cmp(namelist_file, ref_namelist)
        assert os.stat(self.path+'glonf.dat')

    def test_run_model_dat_files(self):
        sami2py.run_model(year=2012, day=211, test=True)
        assert os.stat(self.path+'glonf.dat')
