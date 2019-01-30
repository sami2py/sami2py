'''
    Unit tests for model.py
'''
import os
import shutil
import sami2py
from sami2py.utils import generate_path
from matplotlib.testing.decorators import image_comparison


class TestModelObject():
    '''test basic model object functionality
    '''
    def setup(self):
        '''set up .dat files in properly named director
           for model object to load model
        '''
        self.path = generate_path('test', 256, 1999, 256, test=True)
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat']
        for filename in filelist_fmt:
            open(self.path+'/'+filename, 'w').close()

    def teardown(self):
        '''remove directory tree with the empty .dat files
        '''
        if os.path.isdir(self.path):
            shutil.rmtree(self.path)

    @raises(NameError)
    def test_model_input_exception(self):
        '''If _load_model can't find a file a name error should be produced by
           generate_path in utils.py
        '''
        sami2py.model(tag='none', lon=420, day=420, year=1969)

    def test_model_instantiation(self):
        '''test that model object is instantiated as a sami2py_model
        '''
        S = sami2py.model(tag='test', lon=256, year=1999, day=256)
        assert isinstasnce(S, 'sami2py_model')

    @image_comparison(ionosphere_images=['blank_model_plot'])
    def test_model_plot(self):
        S = sami2py.model(tag='test', lon=256, year=1999, day=256)
        S.plot_lat_alt()

    def test_check_standard_model(self):
        S = sami2py.model(tag='test', lon=256, year=1999, day=256)
        keys = S.check_standard_model()
        assert keys == list()


class TestGetUnformattedData():
    '''Test basic functionality of get_unformatted_data function
    '''
    def test_successful_get(self):
        '''test a successful get of unformatted data
        '''
        return

    @raises()
    def test_reshape_exception(self):
        '''reshape should raise an error if invalid dimensions are provided
        '''
        return

    @raises()
    def file_open_error(self):
        '''file open should raise an error if invalid file path provided
        '''
        return
