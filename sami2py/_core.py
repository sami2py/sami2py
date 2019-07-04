#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Wrapper for running sami2 model

Functions
-------------------------------------------------------------------------------

run_model(tag='model_run', lat=0, lon=0, alt=300, year=2018, day=1,
          f107=120, f107a=120, ap=0,
          rmin=100, rmax=2000, gams=3, gamp=3, altmin=85.,
          dthr=0.25, hrinit=0., hrpr=24., hrmax=48.,
          dt0=30., maxstep=100000000, denmin=1.e-6,
          nion1=1, nion2=7, mmass=48, h_scale=1, o_scale=1,
          no_scale=1, o2_scale=1, he_scale=1, n2_scale=1, n_scale=1,
          Tinf_scale=1, Tn_scale=1., euv_scale=1,
          wind_scale=1, hwm_model=14,
          fejer=True, ExB_drifts=np.zeros((10, 2)), ve01=0., exb_scale=1,
          alt_crit=150., cqe=7.e-14,
          clean=False, test=False, fmtout=True, outn=False)

    Initializes a run of the SAMI2 model and archives the data.
-------------------------------------------------------------------------------

Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 1 Dec 2017, Goddard Space Flight Center (GSFC)
-------------------------------------------------------------------------------
"""
import os
import subprocess
import numpy as np
from sami2py import fortran_dir, __version__
from .utils import generate_path


def run_model(tag='model_run', lat=0, lon=0, alt=300, year=2018, day=1,
              f107=120, f107a=120, ap=0,
              rmin=100, rmax=2000, gams=3, gamp=3, altmin=85.,
              dthr=0.25, hrinit=0., hrpr=24., hrmax=48.,
              dt0=30., maxstep=100000000, denmin=1.e-6,
              nion1=1, nion2=7, mmass=48, h_scale=1, o_scale=1,
              no_scale=1, o2_scale=1, he_scale=1, n2_scale=1, n_scale=1,
              Tinf_scale=1, Tn_scale=1., euv_scale=1,
              wind_scale=1, hwm_model=14,
              fejer=True, ExB_drifts=np.zeros((10, 2)), ve01=0., exb_scale=1,
              alt_crit=150., cqe=7.e-14,
              clean=False, test=False, fmtout=True, outn=False):
    """Runs SAMI2 and archives the data in path

    Parameters
    ----------
    tag : (string)
        Name of run for data archive.  First-level directory under save
        directory
        Note that this is not passed through to sami2 executable
        (default = 'model_run')
    lat : (float)
        latitude intercept of sami2 plane
        (default = 0)
    lon : (float)
        longitude intercept of sami2 plane
        (default = 0)
    alt : (float)
        The input altitude in km.
        (default = 300)
    year : (int)
        year of desired run, integer
    day : (int)
        day of year from Jan 1, acceptable range is [1,366]

    f107 : (float)
        Daily F10.7 solar flux value in SFU
        (default = 120)
    f107a : (float)
        81-day average of F10.7 in SFU
        (default = 120)
    ap : (float)
        quasi-logarithmic geomagnetic index of 3-hour range relative to an
        assumed quiet-day curve.  Integer version of Kp index.
        (default = 0)

    rmin : (float)
        Maximum altitude of the lowest field line in km
        (default = 100)
    rmax : (float)
        Maximum altitude of the highest field line in km
        This has to be less than 20,000 km.
        (default = 2000)
    gams : (int)
        Determines grid spacing along the geomagnetic field. As this
        parameter is increased, the spacing between grid points along
        the field line increases at high altitudes. As it is decreased,
        the spacing becomes more uniform.
        (default=3)
    gamp : (int)
        Determines grid spacing orthogonal to the geomagnetic field.
        As this parameter is increased, the spacing between field lines
        increases at high altitudes. As it is decreased, the spacing
        becomes more uniform.
        (default=3)
    altmin : (float)
        Altitude of the base of a field line (km).
        (default=85)

    dthr : (float)
        Defines how often the data is output (hr).
        (default = 0.25)
    hrinit : (float)
        Universal time at the start of the run (hr).
        (default=0)
    hrpr : (float)
        The time period that elapses before the data is output (hr).
        (default = 24)
    hrmax : (float)
        The number of hours for the run (hr). The first 24 hrs
        allows transients to clear the system.
        (default = 48)
    dt0 : (float)
         The maximum time step allowed (sec). This shouldn't be changed.
        (default=30)
    maxstep : (int)
        The maximum number of time steps allowed.
        (default = 100000000)
    denmin : (float)
        Miniumum ion density allowed.
        (default=1.e-6)

    nion1 : (int)
        Minimum ion specie index.
        1: H+, 2: O+, 3: NO+, 4: O2+, 5: He+, 6: N2+, 7: N+
        (default=1)
    nion2 : (int)
        Maximum ion specie index (see above). One can use 4 and consider
        only the dominant ions in the ionosphere (H, O, NO, O2). This will
        speed up the run time of the code by about 30%.
        (default=7)
    mmass : (int)
        Average neutral mass density.
        (default = 48)

    h_scale : (float)
        Multiplier to scale MSIS neutral H densities
        (default = 1)
    o_scale : (float)
        Multiplier to scale MSIS neutral O densities
        (default = 1)
    no_scale : (float)
        Multiplier to scale MSIS neutral NO densities
        (default = 1)
    o2_scale : (float)
        Multiplier to scale MSIS neutral O2 densities
        (default = 1)
    he_scale : (float)
        Multiplier to scale MSIS neutral He densities
        (default = 1)
    n2_scale : (float)
        Multiplier to scale MSIS neutral N2 densities
        (default = 1)
    n_scale : (float)
        Multiplier to scale MSIS neutral all other densities
        (default = 1)

    Tinf_scale : (float)
        Multiplier to scale Exospheric temperature in MSIS
        (default = 1)
    Tn_scale : (float)
        Multiplier to scale Neutral temperature in MSIS
        (default = 1)
    euv_scale : (float)
        Multiplier to scale total ionization in EUVAC
        (default = 1)
    wind_scale : (float)
        Multiplier to scale Neutral Winds from HWM
        (default = 1)
    hwm_model : (int)
        Specifies which version of HWM to use.
        Allowable values are 93, 7, 14
        (default = 14)

    fejer : (boolean)
        A True value will use the Fejer-Scherliess model of ExB drifts
        A False value will use a user-specified Fourier series for ExB drifts
        (default = True)
    ExB_drifts : (10x2 ndarray of floats)
        Matrix of Fourier series coefficients dependent on solar local time
        (SLT) in hours where
        ExB_total = ExB_drifts[i,0]*cos((i+1)*pi*SLT/12)
                  + ExB_drifts[i,1]*sin((i+1)*pi*SLT/12)
        (default = np.zeros((10,2)))
    ve01 : (float)
        Constant offset for Fourier ExB drifts (m/s)
        (default=0)
    exb_scale : (float)
        Multiplier for ExB model to scale vertical drifts
        (default=1)
    alt_crit : (float)
        The E x B drift is exponentially decreased below this
        altitude with a scale length 20 km.  [This is done to
        allow rmin to be less than 150 km without using an
        extremely small time step.]
        (default=150)
    cqe : (float)
        Constant used in the subroutine 'etemp' associated
        with photoelectron heating. The typical range is
        3e-14 -- 8e-14. The higher this value, the lower
        the electron temperature above 300 km.
        (default=7e-14)

    clean : (boolean)
        A True value will delete the local files after archiving
        A False value will not delete local save files
        Note that this is not passed through to sami2 executable
        (default = False)
    test : (boolean)
        A True value will not run the sami2 executable.  Used for debugging
        the framework.
        A False value will run the sami2 executable
        Note that this is not passed through to sami2 executable
        (default = False)
    fmtout : (boolean)
        If true, sami2 will output as text files.
        If false, sami2 will output as binary.
    outn : (boolean)
        if true, sami2 will include neutral density and wind in output
        if false, sami2 will not include neutral density and wind

    Methods
    ----------
    _generate_namelist(info)

    archive_model(path,clean,fejer)
    """

    current_dir = os.getcwd()
    os.chdir(fortran_dir)

    info = {'year': year, 'day': day, 'lat': lat, 'lon': lon, 'alt': alt,
            'f107': f107, 'f107a': f107a, 'ap': ap,
            'rmin': rmin, 'rmax': rmax, 'gams': gams, 'gamp': gamp,
            'altmin': altmin,
            'dthr': dthr, 'hrinit': hrinit, 'hrpr': hrpr, 'hrmax': hrmax,
            'dt0': dt0, 'maxstep': maxstep, 'denmin': denmin,
            'nion1': nion1, 'nion2': nion2, 'mmass': mmass, 'h_scale': h_scale,
            'o_scale': o_scale, 'no_scale': no_scale, 'o2_scale': o2_scale,
            'he_scale': he_scale, 'n2_scale': n2_scale, 'n_scale': n_scale,
            'exb_scale': exb_scale, 've01': ve01, 'alt_crit': alt_crit,
            'cqe': cqe, 'euv_scale': euv_scale,
            'Tinf_scale': Tinf_scale, 'Tn_scale': Tn_scale,
            'wind_scale': wind_scale, 'hwm_model': hwm_model}

    info['fejer'] = _generate_drift_info(fejer, ExB_drifts)
    info['fmtout'] = _generate_fortran_bool(fmtout)
    info['outn'] = _generate_fortran_bool(outn)
    _generate_namelist(info)
    archive_path = generate_path(tag, lon, year, day, test)
    if not test:
        _ = subprocess.check_call('./sami2py.x')

    _archive_model(archive_path, clean, fejer, fmtout, outn)

    os.chdir(current_dir)


def _generate_drift_info(fejer, ExB_drifts=None):
    """Generates the information regarding the ExB drifts used by the model.
    This information is later stored in the namelist file for SAMI2
    """
    drift_info = _generate_fortran_bool(fejer)
    if not fejer:
        if ExB_drifts.shape != (10, 2):
            raise Exception('Invalid ExB drift shape!  Must be 10x2 ndarray.')
        np.savetxt('exb.inp', ExB_drifts)
    return drift_info


def _generate_fortran_bool(pybool):
    """Generates a fortran bool as a string for the namelist generator
    """
    if pybool:
        fbool = '.true.'
    else:
        fbool = '.false.'
    return fbool


def _generate_namelist(info):
    """Generates namelist file for sami2

    Parameters
    ----------
    info : (dict)
        Contains variables for each line of the namelist file
    """

    # Check HWM model parameters
    if info['hwm_model'] not in [93, 7, 14]:
        print('Invalid HWM Model.  Defaulting to HWM14')
        info['hwm_model'] = 14

    # Print out namelist file
    file = open('sami2py-1.00.namelist', 'w')

    file.write('&go\n')
    file.write(('  fmtout   = {:s},\n').format(info['fmtout']))  # 1
    file.write(('  maxstep  =  {:d},\n').format(info['maxstep']))  # 2
    file.write(('  hrmax    =  {:f},\n').format(info['hrmax']))  # 3
    file.write(('  dt0      =  {:f},\n').format(info['dt0']))  # 4
    file.write(('  dthr     =  {:f},\n').format(info['dthr']))  # 5
    file.write(('  hrpr     =  {:f},\n').format(info['hrpr']))  # 6
    file.write(('  grad_in  =  {:f},\n').format(info['alt']))  # 7
    file.write(('  glat_in  =  {:f},\n').format(info['lat']))  # 8
    file.write(('  glon_in  =  {:f},\n').format(info['lon']))  # 9
    file.write(('  fejer    =  {:s},\n').format(info['fejer']))  # 10
    file.write(('  rmin     =  {:f},\n').format(info['rmin']))  # 11
    file.write(('  rmax     =  {:f},\n').format(info['rmax']))  # 12
    file.write(('  altmin   =  {:f},\n').format(info['altmin']))  # 13
    file.write(('  fbar     =  {:f},\n').format(info['f107a']))  # 14
    file.write(('  f10p7    =  {:f},\n').format(info['f107']))  # 15
    file.write(('  ap       =  {:d},\n').format(info['ap']))  # 16
    file.write(('  year     =  {:d},\n').format(info['year']))  # 17
    file.write(('  day      =  {:d},\n').format(info['day']))  # 18
    file.write(('  mmass    =  {:d},\n').format(info['mmass']))  # 19
    file.write(('  nion1    =  {:d},\n').format(info['nion1']))  # 20
    file.write(('  nion2    =  {:d},\n').format(info['nion2']))  # 21
    file.write(('  hrinit   =  {:f},\n').format(info['hrinit']))  # 22
    file.write(('  tvn0     =  {:f},\n').format(info['wind_scale']))  # 23
    file.write(('  tvexb0   =  {:f},\n').format(info['exb_scale']))  # 24
    file.write(('  ve01     =  {:f},\n').format(info['ve01']))  # 25
    file.write(('  gams     =  {:d},\n').format(info['gams']))  # 26
    file.write(('  gamp     =  {:d},\n').format(info['gamp']))  # 27
    temp_str = '  snn      =  {h:f},{o:f},{no:f},{o2:f},{he:f},{n2:f},{n:f},\n'
    file.write(temp_str.format(h=info['h_scale'],
                               o=info['o_scale'],
                               no=info['no_scale'],
                               o2=info['o2_scale'],
                               he=info['he_scale'],
                               n2=info['n2_scale'],
                               n=info['n_scale']))  # 28
    file.write(('  stn      =  {:f},\n').format(info['Tn_scale']))  # 29
    file.write(('  denmin   =  {:e},\n').format(info['denmin']))  # 30
    file.write(('  alt_crit =  {:f},\n').format(info['alt_crit']))  # 31
    file.write(('  cqe      =  {:e},\n').format(info['cqe']))  # 32
    file.write(('  Tinf_scl =  {:f},\n').format(info['Tinf_scale']))  # 33
    file.write(('  euv_scl  =  {:f},\n').format(info['euv_scale']))  # 34
    file.write(('  hwm_mod  =  {:d},\n').format(info['hwm_model']))  # 35
    file.write(('  outn     = {:s}\n').format(info['outn']))  # 36
    file.write('&end\n')

    file.close()


def _archive_model(path, clean, fejer, fmtout, outn):
    """Moves the model output files to a common archive

    Parameters
    ----------
    path : (string)
        full path of file destination
    clean : (boolean)
        If True, then delete dat files locally
    fejer : (boolean)
        Specifies whether Fejer-Scherliess model is used
        If False, then 'exb.inp' is also archived
    """
    import shutil
    import subprocess

    if fmtout:
        filelist = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                    'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                    'time.dat', 'sami2py-1.00.namelist']
        if outn:
            filelist.append('dennf.dat')
            filelist.append('u4f.dat')
    else:
        filelist = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                    'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                    'time.dat', 'sami2py-1.00.namelist']
        if outn:
            filelist.append('dennu.dat')
            filelist.append('u4u.dat')

    if os.path.isfile(filelist[0]):
        try:
            os.stat(path)
        except FileNotFoundError:
            os.makedirs(path)

        hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
        with open(os.path.join(path, 'version.txt'), 'w+') as f:
            f.write('sami2py v' + __version__ + '\n')
            f.write('short hash ' + hash.decode("utf-8"))

        for list_file in filelist:
            shutil.copyfile(list_file, os.path.join(path, list_file))
        if clean:
            for list_file in filelist[:-1]:
                os.remove(list_file)
        if not fejer:
            shutil.copyfile('exb.inp', os.path.join(path, 'exb.inp'))
    else:
        print('No files to move!')
