#!/usr/bin/env python
#---------------------------------------------------------------------------
# sami2py
#
# Author: Jeff Klenzing, NASA/GSFC,  2017
#
#
# Comments: Tools for specifying, running, and archiving sami2low-1.00
#
# Classes: model
#
# Methods: _generate_path, run_model
#
#---------------------------------------------------------------------------

import platform
import os
import shutil
import numpy as np


class model:


    def __init__(self, tag, lon, year, day):
        """
        Loads a previously run sami2 model and sorts into
           appropriate array shapes

        Input:
            tag  = name of run (top-level directory)
            lon  = longitude reference
            year = year
            day  = day of year

        Properties:
            ut   = Universal Time (hrs), 1D ndarray
            slt  = Solar Local Time (hr), 1D ndarray

            glat = Geographic Latitude (deg), 2D ndarray
            glon = Geographic Longitude (deg), 2D ndarray
            zalt = Altitude (km), 2D ndarray

            deni = Ion density by species, 4D ndarray
            vsi  = Ion Velocity by species, 4D ndarray
            ti   = Ion Temperature by species, 4D ndarray
            te   = Electron Temperature, 3D ndarray
        """

        def _calculate_slt(ut, glon, day):

            slt = np.mod((ut*60 + glon*4),1440)/60.0
            m = 2*np.pi*day/365.242
            dt = -7.657*np.sin(m) + 9.862*np.sin(2*m + 3.599)
            slt = slt - dt/60.0

            return slt

        nf = 98
        nz = 101
        ni = 7

        self.tag = tag
        self.lon0 = lon
        self.year = year
        self.day = day

        path = _generate_path(tag,lon,year,day)

        # Get NameList
        file = open(path + 'sami2low-1.00.namelist')
        self.namelist = file.readlines()
        file.close()

        if self.namelist[10][-6:-2]=='true':
            self.fejer=True
        else:
            self.fejer=False
            self.ExB_drifts = np.loadtxt(path+'exb.inp')
        self.hwm_mod = 'HWM' + self.namelist[36][-3:-1]

        # Get times
        time = np.loadtxt(path+'time.dat')
        self.ut = time[:,1] + time[:,2]/60 + time[:,3]/3600;
        self.slt = _calculate_slt(self.ut,lon,day)
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

    # End __init__ method

    def __repr__(self):

        out = ['']
        out.append('Model Run = %s\n' % self.tag)
        out.append('Day %03d, %4d\n' % (self.day,self.year))
        out.append('Longitude = %d deg\n' % self.lon0)
        out.append('%d time steps from %4.1f to %4.1f UT\n\n'
                % (len(self.ut), min(self.ut), max(self.ut)))

        if self.fejer:
            out.append('Fejer ExB model used\n')
        else:
            out.append('Fourier ExB model used\n')
        out.append('Wind Model used: %s' % self.hwm_mod)

        return ''.join(out)
# End model class


def _generate_sami2_path():
    """
    Creates a path based on platform

    tag: string specifying name of model run
    info: structure to generate namelist
    """

    if platform.node()=='gs674-jklenmbp.home':
        path = '/Users/jklenzin/code/sami2py/fortran/'
    else:
        path = '/Users/jklenzing/code/sami2low-1.00/fortran/'

    return path

# End _generate_path method

def _generate_path(tag, lon, year, day):
    """
    Creates a path based on run tag, date, and longitude

    tag  = string specifying name of model run
    lon  = longitude of model run
    year = year of model run
    day  = day of year of model run
    """

    if platform.node()=='gs674-jklenmbp.home':
        basedir = '/Users/jklenzin/data/sami2/'
    else:
        basedir = '/Volumes/drive/models/sami2/'

    return basedir + tag + ('/lon%03d/%4d_%03d/' % (lon, year, day))

# End _generate_path method


def run_model(year, day, lat=0, lon=0,
              rmin=100, rmax=2000, hrmax=24.5,
              f107=120, f107a=120, ap=0,
              nx=1, ox=1, exb_scale=1,
              fejer=True, ExB_drifts=np.zeros((10,2)),
              Tinf_scale=1, euv_scale=1, hwm_scale=1, hwm_mod=14,
              tag='test', clean=False, test=False):
    """
    Runs SAMI2 and archives the data in path

    Methods: generate_namelist, archive_model

    Inputs

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

        rmin : (float)
            apex altitude of bottom field line in km
            (default = 100)
        rmax : (float)
            apex altitude of top field line in km
            (default = 2000)
        hrmax : (float)
            total time to run model in hours
            (note that output occurs after 24 hours)
            (default = 24.5)

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

        ox : (float)
            Scaled input to modify MSIS neutral monatomic oxygen densities
            (default = 1)
        nx : (float)
            Scaled input to modify MSIS neutral all other densities
            (default = 1)

        exb_scale : (float)
            Multiplier for ExB model to scale vertical drifts
            (default = 1)
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

        Tinf_scale : (float)
            Multiplier to scale Exospheric temperature in MSIS
            (default = 1)
        euv_scale : (float)
            Multiplier to scale total ionization in EUVAC
            (default = 1)
        hwm_scale : (float)
            Multiplier to scale Neutral Winds from HWM
            (default = 1)
        hwm_mod : (int)
            Specifies which version of HWM to use.
            Allowable values are 93, 7, 14
            (default = 14)

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

    """

    def generate_namelist(info):
        """
        Generates namelist file for sami2
        """

        # Check HWM model parameters
        if (info['hwm_mod'] in [93, 7, 14])==False:
            print('Invalid HWM Model.  Defaulting to HWM14')
            info['hwm_mod']=14

        # Print out namelist file

        file = open('sami2low-1.00.namelist','w')

        file.write('&go\n')
        file.write('  fmtout   = .true.,\n')
        file.write('  maxstep  =  100000000,\n')
        file.write('  hrmax    =  %4.1f,\n' % info['hrmx'])
        file.write('  dt0      =  30.,\n')
        file.write('  dthr     =  .2,\n')
        file.write('  hrpr     =  24.,\n')
        file.write('  grad_in  =  300.,\n')
        file.write('  glat_in  =  %6.2f,\n' % info['lat'])
        file.write('  glon_in  =  %6.2f,\n' % info['lon'])
        file.write('  fejer    = ' + info['fejer'] + '\n')
        file.write('  rmin     =  %7.1f,\n' % info['rmin'])
        file.write('  rmax     =  %7.1f,\n' % info['rmax'])
        file.write('  altmin   =   85.,\n')
        file.write('  fbar     =  %5.1f,\n' % info['f107a'])
        file.write('  f10p7    =  %5.1f,\n' % info['f107'])
        file.write('  ap       =  %d,\n' % info['ap'])
        file.write('  year     = %4d,\n' % info['year'])
        file.write('  day      =   %3d,\n' % info['day'])
        file.write('  mmass    =   48 ,\n')
        file.write('  nion1    =    1,\n')
        file.write('  nion2    =    7,\n')
        file.write('  hrinit   =    0.,\n')
        file.write('  tvn0     =    1,\n')
        file.write('  tvexb0   =    %5.2f,\n' % info['exb_scale'])
        file.write('  ve01     =    0.,\n')
        file.write('  gams     =    3,\n')
        file.write('  gamp     =    3,\n')
        file.write('  snn      =    %4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,\n'
                   % (info['nx'],info['ox'],info['nx'],info['nx'],
                   info['nx'],info['nx']))
        file.write('  stn      =    1.,\n')
        file.write('  denmin   =    1.e-6,\n')
        file.write('  alt_crit =    150.,\n')
        file.write('  cqe      =   7.e-14,\n')
        file.write('  Tinf_scl =  %4.2f,\n' % info['Tinf_scale'])
        file.write('  euv_scl  =  %6.2f,\n' % info['euv_scale'])
        file.write('  hwm_scl  =  %6.2f,\n' % info['hwm_scale'])
        file.write('  hwm_mod  = %d\n' % info['hwm_mod'])
        file.write('&end\n')

        file.close()

    # End generate_namelist method

    def archive_model(path,clean,fejer):

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

    # End archive_model method

    current_dir = os.getcwd()
    model_path = _generate_sami2_path()
    os.chdir(model_path)

    info = {'year':year, 'day':day, 'lat':lat, 'lon':lon,
            'hrmx':hrmax, 'rmin':rmin, 'rmax':rmax,
            'f107':f107, 'f107a':f107a, 'ap':ap,
            'nx':nx, 'ox':ox, 'exb_scale':exb_scale,
            'Tinf_scale':Tinf_scale,'euv_scale':euv_scale,
            'hwm_scale':hwm_scale, 'hwm_mod':hwm_mod}
    if fejer:
        info['fejer'] = '.true.'
    else:
        info['fejer'] = '.false.'
        if ExB_drifts.shape!=(10,2):
            print('Invalid ExB drift shape!  Must be 10x2 ndarray.')
        np.savetxt('exb.inp',ExB_drifts)

    generate_namelist(info)
    path = _generate_path(tag,lon,year,day)
    if test==False:
        os.system('./sami2low.x')
    archive_model(path,clean,fejer)

    os.chdir(current_dir)

# End run_model method
