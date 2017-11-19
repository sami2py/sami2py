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
# Methods: run_model
#
#---------------------------------------------------------------------------

import numpy as np
import os, shutil

#class model:
#    def __init__(self, path='/Volumes/drive/sami2/temp/'):

def generate_path(tag, info):
    '''
    Creates a path based on run tag, date, and longitude

    tag: string specifying name of model run
    info: structure to generate namelist
    '''
    import platform

    if platform.node()=='gs674-jklenmbp.home':
        basedir = '/Users/jklenzin/data/sami2/'
    else:
        basedir = '/Volumes/drive/models/sami2/'

    return basedir + tag + '/'


def run_model(year, day, lat=0, lon=0, f107=120, Tinf_scl=1,
              euv_scl=1, hwm_scl=1, hwm_mod=14, tag='test', clean=False):
    '''
    Runs SAMI2 and archives the data in path

    Methods: generate_namelist, archive_model

    '''
    def generate_namelist(info):
        '''
        Generates namelist file for sami2
        '''

        # Check HWM model parameters
        if ~(info['hwm_mod'] in [93, 7, 14]):
            #disp('Invalid HWM Model.  Defaulting to HWM14')
            info['hwm_mod']=14

        # Print out namelist file

        file = open('sami2low-1.00.namelist','w')

        file.write('&go\n')
        file.write('  fmtout   = .true.,\n')
        file.write('  maxstep  =  100000000,\n')
        file.write('  hrmax    =  24.5,\n')
        file.write('  dt0      =  30.,\n')
        file.write('  dthr     =  .2,\n')
        file.write('  hrpr     =  24.,\n')
        file.write('  grad_in  =  300.,\n')
        file.write('  glat_in  =  %6.2f,\n' % info['lat'])
        file.write('  glon_in  =  %6.2f,\n' % info['lon'])
        file.write('  fejer    = .true.\n')
        file.write('  rmin     =  100,\n')
        file.write('  rmax     =  1000,\n')
        file.write('  altmin   =   85.,\n')
        file.write('  fbar     =  117.4778,\n')
        file.write('  f10p7    =  %5.1f,\n' % info['f107'])
        file.write('  ap       =  0,\n')
        file.write('  year     = %4d,\n' % info['year'])
        file.write('  day      =   %3d,\n' % info['day'])
        file.write('  mmass    =   48 ,\n')
        file.write('  nion1    =    1,\n')
        file.write('  nion2    =    7,\n')
        file.write('  hrinit   =    0.,\n')
        file.write('  tvn0     =    1,\n')
        file.write('  tvexb0   =    1,\n')
        file.write('  ve01     =    0.,\n')
        file.write('  gams     =    3,\n')
        file.write('  gamp     =    3,\n')
        file.write('  snn      =    1,1,1,1,1,1,\n')
        file.write('  stn      =    1.,\n')
        file.write('  denmin   =    1.e-6,\n')
        file.write('  alt_crit =    150.,\n')
        file.write('  cqe      =   7.e-14,\n')
        file.write('  Tinf_scl =  %6.2f,\n' % info['Tinf_scl'])
        file.write('  euv_scl  =  %6.2f,\n' % info['euv_scl'])
        file.write('  hwm_scl  =  %6.2f,\n' % info['hwm_scl'])
<<<<<<< HEAD
        file.write('  hwm_mod  = %d\n' % info['hwm_mod'])
=======
        file.write('  hwm_mod    = %d\n' % info['hwm_mod'])
>>>>>>> ab6679b817cc7b4487630f401c9755d8ed92da84
        file.write('&end')

        file.close()

    def archive_model(path='',clean=False):
        filelist = ['glonf.dat','glatf.dat','zaltf.dat',
                    'vsif.dat','time.dat','tif.dat','tef.dat',
                    'denif.dat','sami2low-1.00.namelist']
        try:
            os.stat(path)
        except:
            os.mkdir(path)
        for i in range(0,len(filelist)):
            shutil.copyfile(filelist[i], path+filelist[i])
        if clean:
            for i in range(0,len(filelist)-1):
                os.remove(filelist[i])

    info = {'year':year, 'day':day, 'lat':lat, 'lon':lon, 'f107':f107,
            'Tinf_scl':Tinf_scl,'euv_scl':euv_scl,'hwm_scl':hwm_scl,
            'hwm_mod':hwm_mod}
    generate_namelist(info)
    path = generate_path(tag,info)
    #os.system('./sami2low.x')
    archive_model(path, clean)
