'''
   Unit tests for utils.py
'''

from nose.tools import raises
from sami2py.utils import generate_path, set_archive_dir


class test_generate_path():
    def test_generate_path(self):
        from sami2py import test_data_dir
        out = generate_path(tag='test', lon=0, year=0, day=0, test=True)
        assert out == test_data_dir + 'test/lon000/0000_000/'

    @raises(NameError)
    def test_generate_path_exception(self):
        out = generate_path(tag='test', lon=0, year=0, day=0)


class test_archive_dir():
    def test_set_archive_dir(self):
        return

    def test_set_archive_dir_exception(self):
        return
