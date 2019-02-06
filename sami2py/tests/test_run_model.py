'''
    Unit tests for run_model.py
'''
import sami2py
import os
import filecmp
import shutil
import numpy as np
from nose.tools import raises
from sami2py import fortran_dir, test_data_dir
from sami2py.utils import generate_path


class TestBasicModelRun():
    '''basic tests of the run_model script'''
    def setup(self):
        '''setup function run before each test method to setup files needed
           to run the test effectively
        '''
        self.path = generate_path('test', 0, 2012, 211, True)
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat']
        for filename in filelist_fmt:
            open(fortran_dir + '/' + filename, 'w').close()

    def teardown(self):
        '''teardown function run before each test method to remove files needed
           to run the test effectively
        '''
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat']
        for filename in filelist_fmt:
            os.remove(fortran_dir + '/' + filename)
        if os.path.exists(self.path):
            path_to_remove = os.path.split(self.path)[0]
            path_to_remove = os.path.split(path_to_remove)[0]
            shutil.rmtree(path_to_remove)

    def test_run_model_namelist(self):
        '''the test to ensure that the namelist file is generated properly
        '''
        sami2py.run_model(year=2012, day=211, test=True)
        namelist_file = self.path + 'sami2py-1.00.namelist'
        ref_namelist = test_data_dir + '/reference_sami2py-1.00.namelist'
        assert filecmp.cmp(namelist_file, ref_namelist)

    def test_run_model_dat_files(self):
        '''test to ensure that the dat files are copied properly
        '''
        sami2py.run_model(year=2012, day=211, test=True)
        assert os.stat(self.path + 'glonf.dat')

    @raises(ValueError)
    def test_input_format(self):
        '''test for error output upon incorrect input format
           file.write should throw the error when using string formatting to
           create the file name. Will happen for any variable in the namelist
           set with the wrong type
        '''
        sami2py.run_model(year='2012', day='211', test=True)

    def test_fortran_executable(self):
        '''Short run of fortran executable to ensure the code compiles
           and runs
        '''
        tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=test_data_dir)
        sami2py.run_model(year=2012, day=211,
                          dthr=0.05, hrinit=0.0, hrpr=0.0, hrmax=.11)
        if os.path.isdir(tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=tmp_archive_dir)
        else:
            with open(sami2py.archive_path, 'w') as f:
                f.write('')
                sami2py.archive_dir = ''
