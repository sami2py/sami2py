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
        from sami2py import test_data_dir
        set_archive_path(path=test_data_dir)
        home_dir = os.path.expanduser('~')
        sami2py_dir = os.path.join(home_dir, '.sami2py')
        archive_path = os.path.join(sami2py_dir, 'archive_path.txt')
        with open(archive_path, 'r') as f:
            archive_dir = f.readline()
        assert archive_dir == test_data_dir
        return

    @raises(ValueError)
    def test_set_archive_dir_exception(self):
        set_archive_dir('dummy_invalid_path')
