"""Unit tests for `run_model` in _core.py."""

import numpy as np
import os
import shutil
import warnings

import pytest

import sami2py
from sami2py import fortran_dir
from sami2py import test_data_dir
from sami2py.utils import generate_path


def cmp_lines(path_1, path_2):
    """Compare content of two files."""

    l1 = l2 = True
    with open(path_1, 'r') as f1, open(path_2, 'r') as f2:
        while l1 and l2:
            l1 = f1.readline()
            l2 = f2.readline()
            if l1 != l2:
                return False
    return True


class TestBasicModelRun(object):
    """Basic tests of the run_model script."""

    def setup_method(self):
        """Create a clean testing setup before each method."""

        self.format = True
        self.ref_file = 'ref_f_sami2py-1.00.namelist'
        self.model_path = generate_path(tag='test', lon=0, year=2012, day=211,
                                        test=True)
        if not os.path.exists(self.model_path):
            os.makedirs(self.model_path)
        self.filelist = ['glonf.dat', 'glatf.dat', 'zaltf.dat', 'denif.dat',
                         'dennf.dat', 'u4f.dat', 'vsif.dat', 'tif.dat',
                         'tef.dat', 'time.dat']
        for filename in self.filelist:
            open(os.path.join(fortran_dir, filename), 'w').close()

        return

    def teardown_method(self):
        """Clean up the test env after each method."""

        for filename in self.filelist:
            os.remove(os.path.join(fortran_dir, filename))
        if os.path.exists(self.model_path):
            shutil.rmtree(self.model_path)
        del self.format, self.ref_file, self.model_path, self.filelist

        return

    def test_run_model_namelist(self):
        """Test that the namelist file is generated properly."""

        sami2py.run_model(tag='test', lon=0, year=2012, day=211, test=True,
                          fmtout=self.format)
        namelist_file = os.path.join(self.model_path, 'sami2py-1.00.namelist')
        ref_namelist = os.path.join(test_data_dir, self.ref_file)
        assert cmp_lines(namelist_file, ref_namelist)

        return

    def test_run_model_namelist_w_invalid_hwm(self):
        """Test that an invalid HWM version reverts to version 14."""

        sami2py.run_model(tag='test', lon=0, year=2012, day=211, test=True,
                          fmtout=self.format, hwm_model=15)
        namelist_file = os.path.join(self.model_path, 'sami2py-1.00.namelist')
        ref_namelist = os.path.join(test_data_dir, self.ref_file)
        assert cmp_lines(namelist_file, ref_namelist)

        return

    def test_run_model_dat_files(self):
        """Test that the dat files are copied properly."""

        sami2py.run_model(tag='test', lon=0, year=2012, day=211, test=True,
                          fmtout=self.format, outn=True)
        if self.format:
            fname = 'glonf.dat'
        else:
            fname = 'glonu.dat'
        assert os.stat(os.path.join(self.model_path, fname))

        return

    def test_run_model_exb_files(self):
        """Test that the ExB files are copied properly."""

        sami2py.run_model(tag='test', lon=0, year=2012, day=211, test=True,
                          fmtout=self.format,
                          fejer=False, exb_drifts=np.zeros((10, 2)))
        assert os.stat(os.path.join(self.model_path, 'exb.inp'))

        return

    def test_run_model_exb_wrong_size(self):
        """Test that the ExB has proper shape."""

        with pytest.raises(Exception):
            sami2py.run_model(year=2012, day=211, test=True,
                              fmtout=self.format, fejer=False,
                              exb_drifts=np.zeros((1, 2)))

        return

    def test_input_format(self):
        """Test for error output upon incorrect input format.

        Note
        ----
        file.write should throw the error when using string formatting to
        create the file name. Will happen for any variable in the namelist
        set with the wrong type

        """
        with pytest.raises(ValueError):
            sami2py.run_model(tag='test', year='2012', day='211', test=True,
                              fmtout=self.format)

        return

    def test_fortran_executable(self):
        """Test that fortran executable will compile and run."""

        tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=test_data_dir)
        sami2py.run_model(tag='test', year=2012, day=211, fmtout=self.format,
                          dthr=0.05, hrinit=0.0, hrpr=0.0, hrmax=.11)
        if os.path.isdir(tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=tmp_archive_dir)
        else:
            with open(sami2py._archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

        return


class TestBasicModelRunUnformatted(TestBasicModelRun):
    """Basic tests of the run_model script w/ unformatted output."""

    def setup_method(self):
        """Create a clean testing setup before each method."""

        self.format = False
        self.ref_file = 'ref_u_sami2py-1.00.namelist'
        self.model_path = generate_path(tag='test', lon=0, year=2012, day=211,
                                        test=True)
        if not os.path.exists(self.model_path):
            os.makedirs(self.model_path)
        self.filelist = ['glonu.dat', 'glatu.dat', 'zaltu.dat', 'deniu.dat',
                         'dennu.dat', 'u4u.dat', 'vsiu.dat', 'tiu.dat',
                         'teu.dat', 'time.dat']
        for filename in self.filelist:
            open(os.path.join(fortran_dir, filename), 'w').close()

        return

    def teardown_method(self):
        """Clean up the test env after each method."""

        for filename in self.filelist:
            os.remove(os.path.join(fortran_dir, filename))
        if os.path.exists(self.model_path):
            shutil.rmtree(self.model_path)
        del self.format, self.ref_file, self.model_path, self.filelist


class TestDriftGeneration(object):
    """Tests for the _core function _generate_drift_info."""

    def setup_method(self):
        """Create a clean testing setup before each method."""

        f = open('exb.inp', 'w')
        f.close()

        return

    def teardown_method(self):
        """Clean up the test env after each method."""

        os.remove('exb.inp')

        return

    def test_none_drifts(self):
        """Test usage of None in drift generation."""

        empty_drift = np.zeros((10, 2))
        sami2py._core._generate_drift_info(False, None)
        drifts = np.loadtxt('exb.inp')
        assert np.array_equal(drifts, empty_drift)

        return

    def test_default_drifts(self):
        """Test usage of default in drift generation."""

        default = np.zeros((10, 2))
        default[0, 0] = -30
        sami2py._core._generate_drift_info(False, 'default')
        drifts = np.loadtxt('exb.inp')
        assert np.array_equal(drifts, default)

        return

    def test_bad_string(self):
        """Test that invalid keyward will error."""

        with pytest.raises(Exception):
            sami2py._core._generate_drift_info(False, 'really_cool_drifts')

        return


class TestDeprecation(object):
    """Unit test for deprecation warnings."""

    def setup_method(self):
        """Set up the unit test environment for each method."""

        warnings.simplefilter("always", DeprecationWarning)

        self.model_path = generate_path(tag='test', lon=0, year=2012, day=211,
                                        test=True)
        if not os.path.exists(self.model_path):
            os.makedirs(self.model_path)

        # Set tmp directory
        self.tmp_archive_dir = sami2py.archive_dir
        sami2py.utils.set_archive_dir(path=test_data_dir)

        return

    def teardown_method(self):
        """Clean up the unit test environment after each method."""

        if os.path.isdir(self.tmp_archive_dir):
            sami2py.utils.set_archive_dir(path=self.tmp_archive_dir)
        else:
            with open(sami2py._archive_path, 'w') as archive_file:
                archive_file.write('')
                sami2py.archive_dir = ''

        return

    @pytest.mark.parametrize("key", ['ExB_drifts', 'Tinf_scale', 'Tn_scale'])
    def test_key_deprecation(self, key):
        """Check that camel case variables are deprecated.

        parameters
        ----------
        key : str
            Name of deprecated key

        """

        dep_keys = {'ExB_drifts': np.ones((10, 2)),
                    'Tinf_scale': 0.75,
                    'Tn_scale': 0.75}
        kwargs = {key: dep_keys[key]}

        with warnings.catch_warnings(record=True) as war:
            # Using minimum runtime since the check has occurred before
            # sami2 executable is run
            sami2py.run_model(tag='test', hrmax=0.0, **kwargs)

        warn_msg = "keyword `{:}` is deprecated".format(key)
        msg_found = []
        for item in war:
            msg_found.append(warn_msg in str(item.message))

        # Check that desired warning appears in list of warnings
        assert np.any(msg_found)

        return

    def test_key_error(self):
        """Test that invalid keys error for the deprecated keyword handler."""

        with pytest.raises(KeyError):
            # Using minimum runtime since the check has occurred before
            # sami2 executable is run
            sami2py.run_model(tag='test', hrmax=0.0, dinosaur=True)
        return
