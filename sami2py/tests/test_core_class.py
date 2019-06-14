"""Unit tests for model.py
"""
import os
import sys
import sami2py
from sami2py.utils import generate_path
import numpy as np
from nose.tools import raises


class TestModelObject():
    """Test basic model object functionality
    """
    def setup(self):
        """Set up .dat files in properly named director
           for model object to load model
        """
        self.tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=sami2py.test_data_dir)

    def teardown(self):
        """Undo any changes made to the archive directory
        """
        if os.path.isdir(self.tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=self.tmp_archive_dir)
        else:
            env_dir = sys.prefix
            sami2py_dir = os.path.join(env_dir, '.sami2py')
            archive_path = os.path.join(sami2py_dir, 'archive_path.txt')
            with open(archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

    @raises(IOError)
    def test_model_input_exception(self):
        """File not found error should be produced if the file does not exist
        """
        sami2py.Model(tag='none', lon=428, day=428, year=1969)

    def test_model_instantiation(self):
        """Test that model object is instantiated as a sami2py_model
        """
        model = sami2py.Model(tag='test', lon=256, year=1999, day=256,
                              test=True)
        assert isinstance(model, sami2py.Model)

    def test_model_instantiation_with_unformatted_files(self):
        """Test that model object is instantiated as a sami2py_model
        """
        model = sami2py.Model(tag='test', lon=256, year=1999, day=256,
                              test=True, format=False)
        assert isinstance(model, sami2py.Model)

    def test_model_plot(self):
        """Basic test that a reasonable plot has been created by testing the
           resulting axis limits
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        model = sami2py.Model(tag='test', lon=256, year=1999, day=256,
                              test=True)
        model.plot_lat_alt()
        fig = plt.gcf()
        xlims = fig.axes[0].get_xlim()
        ylims = fig.axes[0].get_ylim()
        assert xlims == (-36.3329, 19.37387)
        assert ylims == (84.98926, 1999.998)
        plt.close()

    def test_check_standard_model(self):
        """Test that standard model outputs nothing if there are no changes to
           the standard model
        """
        model = sami2py.Model(tag='test', lon=256, year=1999, day=256,
                              test=True)
        keys = model.check_standard_model()
        assert keys == list()

    def test_model_repr(self):
        """Test that __repr__ returns a string of information."""
        model = sami2py.Model(tag='test', lon=256, year=1999, day=256,
                              test=True)
        repr_str = model.__repr__()
        assert type(repr_str) is str
