'''
    Unit tests for model.py
'''
import os
import shutil
import sami2py
import numpy as np
import matplotlib.pyplot as plt
from nose.tools import raises
from sami2py.utils import generate_path
from sami2py.model import get_unformatted_data


class TestModelObject():
    '''test basic model object functionality
    '''
    def setup(self):
        '''set up .dat files in properly named director
           for model object to load model
        '''
        self.path = generate_path('test', 256, 1999, 256, test=True)

    @raises(FileNotFoundError)
    def test_model_input_exception(self):
        '''file not found error should be produced if the file does not exist
        '''
        sami2py.model(tag='none', lon=420, day=420, year=1969)

    def test_model_instantiation(self):
        '''test that model object is instantiated as a sami2py_model
        '''
        S = sami2py.model(tag='test', lon=256, year=1999, day=256, test=True)
        assert isinstance(S, sami2py.model)

    def test_model_plot(self):
        '''basic test that a reasonable plot has been created by testing the
           resulting axis limits
        '''
        S = sami2py.model(tag='test', lon=256, year=1999, day=256, test=True)
        S.plot_lat_alt()
        fig = plt.gcf()
        xlims = fig.axes[0].get_xlim()
        ylims = fig.axes[0].get_ylim()
        assert xlims == (-36.3329, 19.37387)
        assert ylims == (84.98926, 1999.998)
        plt.close()

    def test_check_standard_model(self):
        '''test that standard model outputs nothing if there are no changes to
           the standard model
        '''
        S = sami2py.model(tag='test', lon=256, year=1999, day=256, test=True)
        keys = S.check_standard_model()
        assert keys == list()

    def test_successful_get(self):
        '''test a successful get of unformatted data
        '''
        nf = 98
        nz = 101
        ni = 7
        nt = 0
        ret_data = get_unformatted_data(self.path, 'glat', nz, nf, ni, nt)
        glat = np.loadtxt(self.path+'glatf.dat')
        assert ret_data.size == glat.size

    @raises(ValueError)
    def test_reshape_exception(self):
        '''reshape should raise an error if invalid dimensions are provided
        '''
        nf = 500
        nz = 500
        ni = 10
        nt = 0
        get_unformatted_data(self.path, 'glat', nz, nf, ni, nt, reshape=True)

    @raises(FileNotFoundError)
    def file_open_error(self):
        '''file open should raise an error if invalid file path provided
        '''
        nf = 98
        nz = 101
        ni = 7
        nt = 0
        get_unformatted_data(self.path, 'glat', nz, nf, ni, nt)
