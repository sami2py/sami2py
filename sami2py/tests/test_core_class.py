"""Unit tests for `Model` object in _core_class.py."""

import os

import pytest
import xarray as xr

import sami2py


class TestModelObject(object):
    """Test basic model object functionality."""

    def setup(self):
        """Create a clean testing setup before each method."""

        self.tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=sami2py.test_data_dir)
        self.lon = 256.1
        self.year = 1999
        self.day = 256

        return

    def teardown(self):
        """Clean up the test env after each method."""

        if os.path.isdir(self.tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=self.tmp_archive_dir)
        else:
            archive_path = os.path.join(sami2py.sami2py_dir, 'archive_path.txt')
            with open(archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

        return

    def test_model_input_exception(self):
        """Test for error if the file does not exist."""

        with pytest.raises(IOError):
            sami2py.Model(tag='none', lon=428, day=428, year=1969)

        return

    def test_model_instantiation(self):
        """Test that model object is instantiated as a sami2py_model."""

        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True, outn=True)
        assert isinstance(model, sami2py.Model)

        return

    def test_check_standard_model(self):
        """Test the standard model output."""

        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True)
        keys = model.check_standard_model()

        if self.day == 256:
            # This day should be the standard output
            assert keys == list()
        else:
            # This day uses a modified output
            assert 'EUV Multiplier' in keys

        return

    def test_model_repr(self):
        """Test that __repr__ returns a string of information."""

        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True)
        repr_str = model.__repr__()
        assert type(repr_str) is str

        return

    def test_to_netcdf(self):
        """Test that output file is correctly generated."""

        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True, outn=True)
        model.to_netcdf()
        path = 'sami2py_output.nc'
        savedat = xr.load_dataset(path)
        os.remove(path)

        assert model.data == savedat

        return

    def test_to_netcdf_w_path(self):
        """Test that output file is correctly generated."""

        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True, outn=True)
        path = 'custom_filename.nc'
        model.to_netcdf(path=path)
        savedat = xr.load_dataset('custom_filename.nc')
        os.remove(path)

        assert model.data == savedat

        return


class TestModelObjectUnformatted(TestModelObject):
    """Test basic model object functionality."""

    def setup(self):
        """Create a clean testing setup before each method."""

        self.tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=sami2py.test_data_dir)
        self.lon = 256.1
        self.year = 1999
        self.day = 257

        return

    def teardown(self):
        """Clean up the test env after each method."""

        if os.path.isdir(self.tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=self.tmp_archive_dir)
        else:
            archive_path = os.path.join(sami2py.sami2py_dir, 'archive_path.txt')
            with open(archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

        return
