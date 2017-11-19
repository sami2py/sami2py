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


def run_model(info, tag='test', clean=False):
    '''
    Runs SAMI2 and archives the data in path

    Methods: generate_namelist, archive_model

    '''
    def generate_namelist(info=info):
        '''
        Generates namelist file for sami2
        '''
        file = open('sami2low-1.00.namelist','w')

        file.write('&go')
        file.write('  fmtout   = .true.,')
        file.write('  maxstep  =  100000000,')
        file.write('  hrmax    =  24.5,')
        file.write('  dt0      =  30.,')
        file.write('  dthr     =  .2,')
        file.write('  hrpr     =  24.,')
        file.write('  grad_in  =  300.,')
        file.write('  glat_in  =  -19.07,')
        file.write('  glon_in  =  190.07,')
        file.write('  fejer    = .true.')
        file.write('  rmin     =  4120,')
        file.write('  rmax     =  18540,')
        file.write('  altmin   =   85.,')
        file.write('  fbar     =  117.4778,')
        file.write('  f10p7    =  110.9,')
        file.write('  ap       =  9,')
        file.write('  year     = 2012,')
        file.write('  day      =   356,')
        file.write('  mmass    =   48 ,')
        file.write('  nion1    =    1,')
        file.write('  nion2    =    7,')
        file.write('  hrinit   =    0.,')
        file.write('  tvn0     =    1,')
        file.write('  tvexb0   =    1,')
        file.write('  ve01     =    0.,')
        file.write('  gams     =    3,')
        file.write('  gamp     =    3,')
        file.write('  snn      =    1,1,1,1,1,1,')
        file.write('  stn      =    1.,')
        file.write('  denmin   =    1.e-6,')
        file.write('  alt_crit =    150.,')
        file.write('  cqe      =   7.e-14,')
        file.write('  Tinf_scl =  1,')
        file.write('  euv_scl  =  1,')
        file.write('  hwm_scl  =  1,')
        file.write('  hwm_mod    = 93')
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

    #generate_namelist()
    path = generate_path(tag,info)
    os.system('./sami2low.x')
    archive_model(path, clean)
