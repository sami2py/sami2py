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
