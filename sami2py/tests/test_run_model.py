'''
    Unit tests for run_model.py
'''
import sami2py
import os
import numpy as np


class test_basic_model_run():

    def setup(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat', 'sami2low-1.00.namelist']
        filelist_ufmt = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                         'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                         'time.dat', 'sami2low-1.00.namelist']
        for filename in filelist_fmt:
            open(filename, 'w').close()
        for filename in filelist_ufmt:
            open(filename, 'w').close()

    def teardown(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat', 'sami2low-1.00.namelist']
        filelist_ufmt = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                         'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                         'time.dat', 'sami2low-1.00.namelist']
        for filename in filelist_fmt:
            os.remove(filename)
        for filename in filelist_ufmt:
            os.remove(filename)

    def test_run_model_namelist(self):
        sami2py.run_model(year=2012, day=211)
        namelist_file = open(fortran_dir+'sami2low-1.00.namelist', 'rt')
        ref_namelist = open('tests/test_data/reference_sami2low-1.00.namelist',
                            'rt')
        assert namelist_file == ref_namelist
        assert os.stat('tests/test_data/glonf.dat')

    def test_run_model_dat_files(self):
        sami2py.run_model(year=2012, day=211)
        assert os.stat('tests/test_data/glonf.dat')


class test_generate_drift_info():
    def test_exb_file_creation(self):
        exb_drifts = np.zeros(10, 2)
        out = sami2py._generate_drift_info(False, exb_drifts)
        assert os.stat(exb.inp)

    def test_exb_info_str(self):
        exb_drifts = np.zeros(10, 2)
        out = sami2py._generate_drift_info(False, exb_drifts)
        assert out is '.false.'

    @raises(Exception)
    def test_exb_exception(self):
        exb_drifts_wrong = np.zeros(2, 5)
        out = sami2py._generate_drift_info(False, exb_drifts)

    def test_fejer(self):
        out = sami2py._generate_drift_info(True)
        assert out is '.true.'


class test_generate_format_info():
    def test_fmtout_bool(self):
        out = sami2py._generate_format_info(True)
        assert out is '.true.'
        out = sami2py._generate_format_info(False)
        assert out is '.false.'
        return


class test_generate_namelist():
    def setup(self):
        self.info = {'year': 2012, 'day': 211, 'lat': 0, 'lon': 0,
                     'alt': 300, 'f107': 120, 'f107a': 120, 'ap': 0,
                     'rmin': 100, 'rmax': 2000, 'gams': 3, 'gamp': 3,
                     'altmin': 85., 'dthr': 0.25, 'hrinit': 0., 'hrpr': 24.,
                     'hrmax': 48., 'dt0': 30., 'maxstep': 100000000,
                     'denmin': 1.e-6, 'nion1': 1, 'nion2': 7,
                     'mmass': 48, 'h_scale': 1, 'o_scale': 1,
                     'no_scale': 1, 'o2_scale': 1,
                     'he_scale': 1, 'n2_scale': 1,
                     'n_scale': 1, 'exb_scale': 1, 've01': 0,
                     'alt_crit': 150., 'cqe': 7.e-14, 'euv_scale': 1,
                     'Tinf_scale': 1, 'Tn_scale': 1,
                     'wind_scale': 1, 'hwm_model': 14}

    def teardown(self):
        del self.info

    def test_hwm_check(self):
        '''test that the hwm model check works and doesnt
        '''
        os.chdir(test_data)
        self.info['hwm_model'] = 45
        sami2py._generate_namelist(self.info)
        namelist_file = open('sami2low-1.00.namelist', 'rt')
        ref_namelist = open('reference_sami2low-1.00.namelist', 'rt')
        assert namelist_file == ref_namelist

    def test_generate_namelist(self):
        '''assert that the namelist stored in the fortran_dir is the same as
           the nameList in the test_dir
        '''
        os.chdir(test_data)
        sami2py._generate_namelist(self.info)
        namelist_file = open('sami2low-1.00.namelist', 'rt')
        ref_namelist = open('reference_sami2low-1.00.namelist', 'rt')
        assert namelist_file == ref_namelist


class test_basic_archive_model():
    def setup(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat', 'sami2low-1.00.namelist']
        filelist_ufmt = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                         'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                         'time.dat', 'sami2low-1.00.namelist']
        for filename in filelist_fmt:
            open(filename, 'w').close()
        for filename in filelist_ufmt:
            open(filename, 'w').close()

    def teardown(self):
        filelist_fmt = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                        'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                        'time.dat', 'sami2low-1.00.namelist']
        filelist_ufmt = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                         'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                         'time.dat', 'sami2low-1.00.namelist']
        for filename in filelist_fmt:
            os.remove(filename)
        for filename in filelist_ufmt:
            os.remove(filename)

    def test_path_behavior(self):
        '''test behavior if path exists and if it doesnt
           need to make test.dat files for archive_model to work.
        '''
        path = 'test_data/empty_test_dir'
        sami2py._archive_model('test_data/empty_test_dir', False, True, True)
        assert os.stat('test_data/empty_test_dir')
        os.remove('test_data/empty_test_dir')

    @raises(FileNotFoundError)
    def test_cleanup(self):
        '''test that clean does in fact remove files, what could go wrong
        '''
        sami2py._archive_model('test_data/empty_test_dir', True, True, True)
        os.stat('glonf.dat')


class test_archive_model_exceptions():
    @raises(FileNotFoundError)
    def test_copy_no_file(self):
        '''test proper copying and relocation of dat files, what could go wrong
        '''
        sami2py.archive_model('test_data/empty_test_dir', False, True, True)

    @raises(FileNotFoundError)
    def test_non_fejer(self):
        '''test that exb.inp is properly moved or not depending on flag
        '''
        sami2py.archive_model('test_data/empty_test_dir', False, False, True)
