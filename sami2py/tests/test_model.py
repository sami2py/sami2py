'''
    Unit tests for model.py
'''
import os
import shutil


class test_model_object():
    def setup():
        '''set up .dat files in properly named director
           for model object to load model
        '''
        self.path = generate_path('test', 256, 1999, 256, True)
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat']
        for filename in filelist_fmt:
            open(self.path+'/'+filename, 'w').close()

    def teardown():
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
        return

    def test_model_plot(self):
        return

    def test_check_standard_model(self):
        return


class test_get_unformatted_data():
    def test_successful_get(self):
        return

    @raises()
    def test_reshape_exception(self):
        return

    @raises()
    def file_open_error(self):
        return
