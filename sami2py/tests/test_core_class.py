"""Unit tests for model.py
"""
import os
import sami2py
import pytest


class TestModelObject():
    """Test basic model object functionality
    """
    def setup(self):
        """Set up .dat files in properly named director
           for model object to load model
        """
        self.tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=sami2py.test_data_dir)
        self.lon = 256
        self.year = 1999
        self.day = 256

    def teardown(self):
        """Undo any changes made to the archive directory
        """
        if os.path.isdir(self.tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=self.tmp_archive_dir)
        else:
            archive_path = os.path.join(sami2py.sami2py_dir, 'archive_path.txt')
            with open(archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

    def test_model_input_exception(self):
        """File not found error should be produced if the file does not exist
        """
        with pytest.raises(IOError):
            sami2py.Model(tag='none', lon=428, day=428, year=1969)

    def test_model_instantiation(self):
        """Test that model object is instantiated as a sami2py_model
        """
        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True, outn=True)
        assert isinstance(model, sami2py.Model)

    def test_model_plot(self):
        """Basic test that a reasonable plot has been created by testing the
           resulting axis limits
        """
        import matplotlib

        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True)
        fig = model.plot_lat_alt()
        assert isinstance(fig, matplotlib.figure.Figure)

        tol = 1.e-4
        xlims = fig.axes[0].get_xlim()
        ylims = fig.axes[0].get_ylim()
        assert abs(xlims[0] + 36.3329) < tol
        assert abs(xlims[1] - 19.37387) < tol
        assert abs(ylims[0] - 84.98926) < tol
        assert abs(ylims[1] - 1999.998) < tol
        matplotlib.pyplot.close(fig)

    def test_check_standard_model(self):
        """Test that standard model outputs nothing if there are no changes to
           the standard model / changes to EUV in unformatted version
        """
        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True)
        keys = model.check_standard_model()

        if self.day == 256:
            assert keys == list()
        else:
            assert 'EUV Multiplier' in keys

    def test_model_repr(self):
        """Test that __repr__ returns a string of information."""
        model = sami2py.Model(tag='test', lon=self.lon, year=self.year,
                              day=self.day, test=True)
        repr_str = model.__repr__()
        assert type(repr_str) is str


class TestModelObjectUnformatted(TestModelObject):
    """Test basic model object functionality
    """
    def setup(self):
        """Set up .dat files in properly named director
           for model object to load model
        """
        self.tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=sami2py.test_data_dir)
        self.lon = 256
        self.year = 1999
        self.day = 257

    def teardown(self):
        """Undo any changes made to the archive directory
        """
        if os.path.isdir(self.tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=self.tmp_archive_dir)
        else:
            archive_path = os.path.join(sami2py.sami2py_dir, 'archive_path.txt')
            with open(archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''
