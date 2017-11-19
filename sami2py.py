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
    #def generate_namelist(info=info):

    def archive_model(path='',clean):
        filelist = ['glonf.dat','glatf.dat','zaltf.dat',
                    'vsif.dat','time.dat','tif.dat','tef.dat',
                    'denif.dat','sami2low-1.00.namelist']
        if ~os.direxists(path):
            os.mkdir(path)
        for i in range(0,len(filelist)):
            shutil.copyfile(filelist[i], path+filelist[i])
        if clean:
            for i in range(0,len(filelist)-1):
                os.remove(filelist[i])

    #generate_namelist()
    path = generate_path(tag,info)
    os.system(''./sami2low.x')
    archive_model(path, clean)
