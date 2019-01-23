'''
   Unit tests for utils.py
'''

import os
import sami2py
from nose.tools import raises
from sami2py.utils import generate_path, set_archive_dir


class test_generate_path():
    def test_generate_path(self):
        '''test the basic functionality of generating a path
        '''
        from sami2py import test_data_dir
        out = generate_path(tag='test', lon=0, year=1900, day=0, test=True)
        assert out == test_data_dir + 'test/lon000/1900_000/'

    @raises(NameError)
    def test_generate_path_exception(self):
        '''if no archive directory has been specified a name error should be
           produced
        '''
        tmp_archive_dir = sami2py.archive_dir
        sami2py.archive_dir = ''
        out = generate_path(tag='test', lon=0, year=0, day=0)
        sami2py.archive_dir = tmp_archive_dir


class test_archive_dir():
    def test_set_archive_dir(self):
        '''test that set_archive_dir has set and stored the archive directory
        '''
        from sami2py import test_data_dir
        tmp_archive_dir = sami2py.archive_dir
        set_archive_dir(path=test_data_dir)
        home_dir = os.path.expanduser('~')
        sami2py_dir = os.path.join(home_dir, '.sami2py')
        archive_path = os.path.join(sami2py_dir, 'archive_path.txt')
        with open(archive_path, 'r') as f:
            archive_dir = f.readline()
        assert archive_dir == test_data_dir
        # return the archive dir to its previous value
        set_archive_dir(path=tmp_archive_dir)

    @raises(ValueError)
    def test_set_archive_dir_exception(self):
        '''if the provided path is invalid a value error should be produced
        '''
        set_archive_dir('dummy_invalid_path')
