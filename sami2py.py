#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
#-----------------------------------------------------------------------------
""" Tools for calculating integrated solar indices.

Functions
-------------------------------------------------------------------------------

def run_model(year, day, lat=0, lon=0, alt=300,
              f107=120, f107a=120, ap=0,
              rmin=100, rmax=2000, gams=3, gamp=3, altmin=85.,
              dthr=0.25, hrinit=0., hrpr=24., hrmax=48.,
              dt0=30., maxstep=100000000, denmin=1.e-6,
              nion1=1, nion2=7, mmass=48.0, h_scale=1, o_scale=1,
              no_scale=1, o2_scale=1, he_scale=1, n2_scale=1, n_scale=1,
              Tinf_scale=1, Tn_scale=1., euv_scale=1,
              wind_scale=1, hwm_model=14,
              fejer=True, ExB_drifts=np.zeros((10,2)), ve01=0., exb_scale=1,
              alt_crit=150., cqe=7.e-14,
              tag='test', clean=False, test=False)

    Initializes a run of the SAMI2 model and archives the data.


_generate_path(tag, lon, year, day)
    Generates path to archived model runs based on input paramters.

Classes
-------------------------------------------------------------------------------

model   Loads, reshapes, and holds SAMI2 output for a given model run
        specified by the user.

Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 22 Nov 2017, Goddard Space Flight Center (GSFC)

References
-------------------------------------------------------------------------------


"""

import platform
import os
import shutil
import numpy as np


class model:


    def __init__(self, tag, lon, year, day):
        """ Loads a previously run sami2 model and sorts into
            appropriate array shapes

        Parameters
        ----------
        tag : (string)
            name of run (top-level directory)
        lon : (int)
            longitude reference
        year : (int)
            year
        day : (int)
            day of year from Jan 1

        Returns
        ---------
        self : model class object containing OCB file data

        Attributes
        ----------
        ut : (1D ndarray)
            Universal Time (hrs)
        slt : (1D ndarray)
            Solar Local Time (hr)

        glat : (2D ndarray)
            Geographic Latitude (deg)
        glon : (2D ndarray)
            Geographic Longitude (deg)
        zalt : (2D ndarray)
            Altitude (km)

        deni : (4D ndarray)
            Ion density by species (cm^-3)
        vsi : (4D ndarray)
            Ion Velocity by species (m/s)
        ti : (4D ndarray)
            Ion Temperature by species (K)
        te : (3D ndarray)
            Electron Temperature (K)
        """

        self.tag = tag
        self.lon0 = lon
        self.year = year
        self.day = day

        self._load_model()

    # End __init__ method

    def __repr__(self):

        out = ['']
        out.append('Model Run = %s\n' % self.tag)
        out.append('Day %03d, %4d\n' % (self.day,self.year))
        out.append('Longitude = %d deg\n' % self.lon0)
        out.append('%d time steps from %4.1f to %4.1f UT\n\n'
                % (len(self.ut), min(self.ut), max(self.ut)))

        out.append('Solar Activity\n')
        out.append('--------------\n')
        out.append('F10.7: %5.1f sfu\n' % self.MetaData['F10.7'])
        out.append('F10.7A: %5.1f sfu\n' % self.MetaData['F10.7A'])
        out.append('ap: %d \n\n' % self.MetaData['ap'])

        out.append('Component Models Used\n')
        out.append('---------------------\n')
        out.append('Neutral Atmosphere: %s\n' % self.MetaData['Neutral Atmosphere Model'])
        out.append('Winds: %s\n' % self.MetaData['Wind Model'])
        out.append('Photoproduction: %s\n' % self.MetaData['EUV Model'])
        out.append('ExB Drifts: %s\n\n' % self.MetaData['ExB model'])

        return ''.join(out)

    def _calculate_slt(self):
        """ Calculates Solar Local Time for reference point of model

        Parameters
        ----------
        None

        Returns
        -------
        self.slt : (float)
            Solar Local Time in hours

        """

        slt = np.mod((self.ut*60 + self.lon0*4),1440)/60.0
        m = 2*np.pi*self.day/365.242
        dt = -7.657*np.sin(m) + 9.862*np.sin(2*m + 3.599)
        self.slt = slt - dt/60.0

    def _load_model(self):
        """ Loads model results

        Parameters
        ----------
        None

        Returns
        -------

        """

        nf = 98
        nz = 101
        ni = 7

        path = _generate_path(self.tag,self.lon0,self.year,self.day)

        # Get NameList
        file = open(path + 'sami2low-1.00.namelist')
        self.namelist = file.readlines()
        file.close()

        self.MetaData = dict()
        self._generate_metadata(self.namelist)

        # Get times
        time = np.loadtxt(path+'time.dat')
        self.ut = time[:,1] + time[:,2]/60 + time[:,3]/3600;
        self._calculate_slt()
        nt = len(self.ut)

        # Get Location
        glat = np.loadtxt(path+'glatf.dat')
        self.glat = np.reshape(glat,(nz,nf),order="F")
        glon = np.loadtxt(path+'glonf.dat')
        self.glon = np.reshape(glon,(nz,nf),order="F")
        zalt = np.loadtxt(path+'zaltf.dat')
        self.zalt = np.reshape(zalt,(nz,nf),order="F")
        del glat, glon, zalt

        # Get plasma values
        deni = np.loadtxt(path+'denif.dat')
        self.deni = np.reshape(deni,(nz,nf,ni,nt),order="F")
        vsi = np.loadtxt(path+'vsif.dat')
        self.vsi = np.reshape(vsi,(nz,nf,ni,nt),order="F")
        ti = np.loadtxt(path+'tif.dat')
        self.ti = np.reshape(ti,(nz,nf,ni,nt),order="F")
        te = np.loadtxt(path+'tef.dat')
        self.te = np.reshape(te,(nz,nf,nt),order="F")
        del deni, vsi, ti, te

    def _generate_metadata(self,namelist):
        """ Reads the namelist and generates MetaData based on Parameters
        """

        import re

        self.MetaData['F10.7A'] = float(re.findall(r"\d*\.\d+|\d+", namelist[14])[0])
        self.MetaData['F10.7'] = float(re.findall(r"\d*\.\d+|\d+", namelist[15])[2])
        self.MetaData['ap'] = int(re.findall(r"\d+", namelist[16])[0])

        self.MetaData['Neutral Atmosphere Model'] = 'NRLMSISe-2000'
        self.MetaData['EUV Model'] = 'EUVAC'
        neutral_scalars = re.findall(r"\d*\.\d+|\d+", namelist[28])
        self.MetaData['H Multiplier'] = float(neutral_scalars[0])
        self.MetaData['O Multiplier'] = float(neutral_scalars[1])
        self.MetaData['NO Multiplier'] = float(neutral_scalars[2])
        self.MetaData['O2 Multiplier'] = float(neutral_scalars[3])
        self.MetaData['He Multiplier'] = float(neutral_scalars[4])
        self.MetaData['N2 Multiplier'] = float(neutral_scalars[5])
        self.MetaData['N Multiplier'] = float(neutral_scalars[6])

        if '.true.' in namelist[10]:
            self.MetaData['ExB model'] = 'Fejer-Scherliess'
        else:
            self.MetaData['ExB model'] = 'Fourier Series'
            self.MetaData['Fourier Coeffs'] = np.loadtxt(path+'exb.inp')

        wind_model = int(re.findall(r"\d+",namelist[36])[0])
        self.MetaData['Wind Model'] = ('HWM-%02d' % wind_model)

def _generate_path(tag, lon, year, day):
    """
    Creates a path based on run tag, date, and longitude

    Parameters
    ----------
    tag : (string)
        specifies name of model run
    lon : (int)
        longitude of model run
    year : (int)
        year of model run
    day : (int)
        day of year of model run

    Returns
    -------
    path : (string)
        Complete path pointing to model archive for a given run
    """

    if platform.node()=='gs674-jklenmbp.home':
        basedir = '/Users/jklenzin/data/sami2/'
    else:
        basedir = '/Volumes/drive/models/sami2/'

    path = basedir + tag + ('/lon%03d/%4d_%03d/' % (lon, year, day))
    return path


def run_model(year, day, lat=0, lon=0, alt=300,
              f107=120, f107a=120, ap=0,
              rmin=100, rmax=2000, gams=3, gamp=3, altmin=85.,
              dthr=0.25, hrinit=0., hrpr=24., hrmax=48.,
              dt0=30., maxstep=100000000, denmin=1.e-6,
              nion1=1, nion2=7, mmass=48.0, h_scale=1, o_scale=1,
              no_scale=1, o2_scale=1, he_scale=1, n2_scale=1, n_scale=1,
              Tinf_scale=1, Tn_scale=1., euv_scale=1,
              wind_scale=1, hwm_model=14,
              fejer=True, ExB_drifts=np.zeros((10,2)), ve01=0., exb_scale=1,
              alt_crit=150., cqe=7.e-14,
              tag='test', clean=False, test=False):
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
        Local time at the start of the run (hr).
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
    mmass : (float)
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


    Methods
    ----------
    _generate_namelist(info)

    archive_model(path,clean,fejer)

    """

    def _generate_namelist(info):
        """
        Generates namelist file for sami2

        Parameters
        ----------
        info : (dict)
            Contains variables for each line of the namelist file
        """

        # Check HWM model parameters
        if (info['hwm_model'] in [93, 7, 14])==False:
            print('Invalid HWM Model.  Defaulting to HWM14')
            info['hwm_model']=14

        # Print out namelist file

        file = open('sami2low-1.00.namelist','w')

        file.write('&go\n')
        file.write('  fmtout   = .true.,\n')
        file.write('  maxstep  =  %d,\n' % info['maxstep'])
        file.write('  hrmax    =  %f,\n' % info['hrmax'])
        file.write('  dt0      =  %f,\n' % info['dt0'])
        file.write('  dthr     =  %f,\n' % info['dthr'])
        file.write('  hrpr     =  %f,\n' % info['hrpr'])
        file.write('  grad_in  =  %f,\n' % info['alt'])
        file.write('  glat_in  =  %f,\n' % info['lat'])
        file.write('  glon_in  =  %f,\n' % info['lon'])
        file.write('  fejer    =  %s,\n' % info['fejer'])
        file.write('  rmin     =  %f,\n' % info['rmin'])
        file.write('  rmax     =  %f,\n' % info['rmax'])
        file.write('  altmin   =  %f,\n' % info['altmin'])
        file.write('  fbar     =  %f,\n' % info['f107a'])
        file.write('  f10p7    =  %f,\n' % info['f107'])
        file.write('  ap       =  %d,\n' % info['ap'])
        file.write('  year     =  %d,\n' % info['year'])
        file.write('  day      =  %d,\n' % info['day'])
        file.write('  mmass    =  %f,\n' % info['mmass'])
        file.write('  nion1    =  %d,\n' % info['nion1'])
        file.write('  nion2    =  %d,\n' % info['nion2'])
        file.write('  hrinit   =  %f,\n' % info['hrinit'])
        file.write('  tvn0     =  %f,\n' % info['wind_scale'])
        file.write('  tvexb0   =  %f,\n' % info['exb_scale'])
        file.write('  ve01     =  %f,\n' % info['ve01'])
        file.write('  gams     =  %d,\n' % info['gams'])
        file.write('  gamp     =  %d,\n' % info['gamp'])
        file.write('  snn      =  %f,%f,%f,%f,%f,%f,%f,\n'
                   % (info['h_scale'],
                      info['o_scale'],
                      info['no_scale'],
                      info['o2_scale'],
                      info['he_scale'],
                      info['n2_scale'],
                      info['n_scale']))
        file.write('  stn      =  %f,\n' % info['Tn_scale'])
        file.write('  denmin   =  %e,\n' % info['denmin'])
        file.write('  alt_crit =  %f,\n' % info['alt_crit'])
        file.write('  cqe      =  %e,\n' % info['cqe'])
        file.write('  Tinf_scl =  %f,\n' % info['Tinf_scale'])
        file.write('  euv_scl  =  %f,\n' % info['euv_scale'])
        file.write('  hwm_scl  =  %f,\n' % info['wind_scale']) # Duplicate!
        file.write('  hwm_mod  =  %d\n' % info['hwm_model'])
        file.write('&end\n')

        file.close()


    def _generate_sami2_path():
        """
        Creates a path based on platform

        parameters
        ----------
        None

        Returns
        -------
        path : (string)
            complete path pointing to sami2 fortran and executables

        """

        if platform.node()=='gs674-jklenmbp.home':
            path = '/Users/jklenzin/code/sami2py/fortran/'
        else:
            path = '/Users/jklenzing/code/sami2low-1.00/fortran/'

        return path


    def archive_model(path,clean,fejer):
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

        filelist = ['glonf.dat','glatf.dat','zaltf.dat',
                    'vsif.dat','time.dat','tif.dat','tef.dat',
                    'denif.dat','sami2low-1.00.namelist']
        try:
            os.stat(path)
        except:
            os.makedirs(path)

        for i in range(0,len(filelist)):
            shutil.copyfile(filelist[i], path+filelist[i])
        if clean:
            for i in range(0,len(filelist)-1):
                os.remove(filelist[i])
        if fejer==False:
            shutil.copyfile('exb.inp', path+'exb.inp')


    current_dir = os.getcwd()
    model_path = _generate_sami2_path()
    os.chdir(model_path)

    info = {'year':year, 'day':day, 'lat':lat, 'lon':lon, 'alt':alt,
            'f107':f107, 'f107a':f107a, 'ap':ap,
            'rmin':rmin, 'rmax':rmax, 'gams':gams, 'gamp':gamp, 'altmin':altmin,
            'dthr':dthr, 'hrinit':hrinit, 'hrpr':hrpr, 'hrmax':hrmax,
            'dt0':dt0, 'maxstep':maxstep, 'denmin':denmin,
            'nion1':nion1, 'nion2':nion2, 'mmass':mmass,'h_scale':h_scale,
            'o_scale':o_scale, 'no_scale':no_scale, 'o2_scale':o2_scale,
            'he_scale':he_scale, 'n2_scale':n2_scale, 'n_scale':n_scale,
            'exb_scale':exb_scale, 've01':ve01, 'alt_crit':alt_crit, 'cqe':cqe,
            'Tinf_scale':Tinf_scale, 'Tn_scale':Tn_scale, 'euv_scale':euv_scale,
            'wind_scale':wind_scale, 'hwm_model':hwm_model}
    if fejer:
        info['fejer'] = '.true.'
    else:
        info['fejer'] = '.false.'
        if ExB_drifts.shape!=(10,2):
            print('Invalid ExB drift shape!  Must be 10x2 ndarray.')
        np.savetxt('exb.inp',ExB_drifts)


    _generate_namelist(info)
    path = _generate_path(tag,lon,year,day)
    if test==False:
        os.system('./sami2low.x')
    archive_model(path,clean,fejer)

    os.chdir(current_dir)
