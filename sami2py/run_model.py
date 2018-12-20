#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
#-----------------------------------------------------------------------------
""" Wrapper for running sami2 model

Functions
-------------------------------------------------------------------------------

run_model(year, day, lat=0, lon=0, alt=300,
              f107=120, f107a=120, ap=0,
              rmin=100, rmax=2000, gams=3, gamp=3, altmin=85.,
              dthr=0.25, hrinit=0., hrpr=24., hrmax=48.,
              dt0=30., maxstep=100000000, denmin=1.e-6,
              nion1=1, nion2=7, mmass=48, h_scale=1, o_scale=1,
              no_scale=1, o2_scale=1, he_scale=1, n2_scale=1, n_scale=1,
              Tinf_scale=1, Tn_scale=1., euv_scale=1,
              wind_scale=1, hwm_model=14,
              fejer=True, ExB_drifts=np.zeros((10,2)), ve01=0., exb_scale=1,
              alt_crit=150., cqe=7.e-14,
              tag='test', clean=False, test=False)

    Initializes a run of the SAMI2 model and archives the data.
-------------------------------------------------------------------------------

Classes
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 1 Dec 2017, Goddard Space Flight Center (GSFC)
-------------------------------------------------------------------------------

References
-------------------------------------------------------------------------------


"""
import os
import numpy as np
from sami2py import fortran_dir
from .utils import generate_path


def run_model(year, day, lat=0, lon=0, alt=300,
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
              tag='test', clean=False, test=False, fmtout=True):
    """
    Runs SAMI2 and archives the data in path

    Parameters
    ----------
    year : (int)
        year of desired run, integer
    day : (int)
        day of year from Jan 1, acceptable range is [1,366]
    lat : (float)
        latitude intercept of sami2 plane
        (default = 0)
    lon : (float)
        longitude intercept of sami2 plane
        (default = 0)
    alt : (float)
        The input altitude in km.
        (default = 300)

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

    tag : (string)
        Name of run for data archive.  First-level directory under save directory
        (default = 'test')
    clean : (boolean)
        A True value will delete the local files after archiving
        A False value will not delete local save files
        (default = False)
    test : (boolean)
        A True value will not run the sami2 executable.  Used for debugging the framework.
        A False value will run the sami2 executable
        (default = False)
    fmtout : (boolean)
        If true, sami2 will output as text files.
        If false, sami2 will output as binary.


    Methods
    ----------
    _generate_namelist(info)

    archive_model(path,clean,fejer)
    """

    current_dir = os.getcwd()
    os.chdir(fortran_dir)

    info = {'year':year, 'day':day, 'lat':lat, 'lon':lon, 'alt':alt,
            'f107':f107, 'f107a':f107a, 'ap':ap,
            'rmin':rmin, 'rmax':rmax, 'gams':gams, 'gamp':gamp, 'altmin':altmin,
            'dthr':dthr, 'hrinit':hrinit, 'hrpr':hrpr, 'hrmax':hrmax,
            'dt0':dt0, 'maxstep':maxstep, 'denmin':denmin,
            'nion1':nion1, 'nion2':nion2, 'mmass':mmass, 'h_scale':h_scale,
            'o_scale':o_scale, 'no_scale':no_scale, 'o2_scale':o2_scale,
            'he_scale':he_scale, 'n2_scale':n2_scale, 'n_scale':n_scale,
            'exb_scale':exb_scale, 've01':ve01, 'alt_crit':alt_crit, 'cqe':cqe,
            'Tinf_scale':Tinf_scale, 'Tn_scale':Tn_scale, 'euv_scale':euv_scale,
            'wind_scale':wind_scale, 'hwm_model':hwm_model}
    if fejer:
        info['fejer'] = '.true.'
    else:
        info['fejer'] = '.false.'
        if ExB_drifts.shape != (10, 2):
            print('Invalid ExB drift shape!  Must be 10x2 ndarray.')
        np.savetxt('exb.inp', ExB_drifts)

    if fmtout:
        info['fmtout'] = '.true.'
    else:
        info['fmtout'] = '.false.'


    _generate_namelist(info)
    path = generate_path(tag, lon, year, day)
    if not test:
        os.system('./sami2low.x')
    _archive_model(path, clean, fejer, fmtout)

    os.chdir(current_dir)


def _generate_namelist(info):
    """
    Generates namelist file for sami2

    Parameters
    ----------
    info : (dict)
        Contains variables for each line of the namelist file
    """

    # Check HWM model parameters
    if not (info['hwm_model'] in [93, 7, 14]):
        print('Invalid HWM Model.  Defaulting to HWM14')
        info['hwm_model'] = 14

    # Print out namelist file
    file = open('sami2low-1.00.namelist', 'w')

    file.write('&go\n')
    file.write('  fmtout   = %s,\n' % info['fmtout']) #1
    file.write('  maxstep  =  %d,\n' % info['maxstep']) #2
    file.write('  hrmax    =  %f,\n' % info['hrmax']) #3
    file.write('  dt0      =  %f,\n' % info['dt0']) #4
    file.write('  dthr     =  %f,\n' % info['dthr']) #5
    file.write('  hrpr     =  %f,\n' % info['hrpr']) #6
    file.write('  grad_in  =  %f,\n' % info['alt']) #7
    file.write('  glat_in  =  %f,\n' % info['lat']) #8
    file.write('  glon_in  =  %f,\n' % info['lon']) #9
    file.write('  fejer    =  %s,\n' % info['fejer']) #10
    file.write('  rmin     =  %f,\n' % info['rmin']) #11
    file.write('  rmax     =  %f,\n' % info['rmax']) #12
    file.write('  altmin   =  %f,\n' % info['altmin']) #13
    file.write('  fbar     =  %f,\n' % info['f107a']) #14
    file.write('  f10p7    =  %f,\n' % info['f107']) #15
    file.write('  ap       =  %d,\n' % info['ap']) #16
    file.write('  year     =  %d,\n' % info['year']) #17
    file.write('  day      =  %d,\n' % info['day']) #18
    file.write('  mmass    =  %d,\n' % info['mmass']) #19
    file.write('  nion1    =  %d,\n' % info['nion1']) #20
    file.write('  nion2    =  %d,\n' % info['nion2'])  #21
    file.write('  hrinit   =  %f,\n' % info['hrinit']) #22
    file.write('  tvn0     =  %f,\n' % info['wind_scale']) #23
    file.write('  tvexb0   =  %f,\n' % info['exb_scale']) #24
    file.write('  ve01     =  %f,\n' % info['ve01']) #25
    file.write('  gams     =  %d,\n' % info['gams']) #26
    file.write('  gamp     =  %d,\n' % info['gamp']) #27
    file.write('  snn      =  %f,%f,%f,%f,%f,%f,%f,\n'
               % (info['h_scale'],
                  info['o_scale'],
                  info['no_scale'],
                  info['o2_scale'],
                  info['he_scale'],
                  info['n2_scale'],
                  info['n_scale'])) #28
    file.write('  stn      =  %f,\n' % info['Tn_scale']) #29
    file.write('  denmin   =  %e,\n' % info['denmin']) #30
    file.write('  alt_crit =  %f,\n' % info['alt_crit']) #31
    file.write('  cqe      =  %e,\n' % info['cqe']) #32
    file.write('  Tinf_scl =  %f,\n' % info['Tinf_scale']) #33
    file.write('  euv_scl  =  %f,\n' % info['euv_scale']) #34
    #file.write('  hwm_scl  =  %f,\n' % info['wind_scale']) # Duplicate!
    file.write('  hwm_mod  =  %d\n' % info['hwm_model']) #35
    file.write('&end\n')

    file.close()


def _archive_model(path, clean, fejer, fmtout):
    """ Moves the model output files to a common archive

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

    if fmtout:
        filelist = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                    'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                    'time.dat', 'sami2low-1.00.namelist']
    else:
        filelist = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                    'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                    'time.dat', 'sami2low-1.00.namelist']

    try:
        os.stat(path)
    except FileNotFoundError:
        os.makedirs(path)

    for list_file in filelist:
        shutil.copyfile(list_file, path+list_file)
    if clean:
        for list_file in filelist[:-1]:
            os.remove(list_file)
    if not fejer:
        shutil.copyfile('exb.inp', path+'exb.inp')
