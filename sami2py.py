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


def run_model(info, path='/Volumes/drive/sami2/temp/'):
    '''
    Runs SAMI2 and archives the data in path

    Methods: generate_namelist, archive_model

    '''
    #def generate_namelist(info=info):

    def archive_model(path=path):
        filelist = ['glonf.dat','glatf.dat','zaltf.dat',
                    'vsif.dat','time.dat','tif.dat','tef.dat',
                    'denif.dat','sami2low-1.00.namelist']
        for i in range(0,len(filelist)):
            shutil.copyfile(filelist[i], path+filelist[i])

    #generate_namelist()
    os.system("./sami2low.x")
    archive_model()
