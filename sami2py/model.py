#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
#-----------------------------------------------------------------------------
""" Wrapper for running sami2 model

Functions
-------------------------------------------------------------------------------


Classes
-------------------------------------------------------------------------------

model
    Loads, reshapes, and holds SAMI2 output for a given model run
    specified by the user.


Moduleauthor
-------------------------------------------------------------------------------
Jeff Klenzing (JK), 1 Dec 2017, Goddard Space Flight Center (GSFC)

References
-------------------------------------------------------------------------------


"""
from .utils import generate_path
import numpy as np

class model(object):


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

        path = generate_path(self.tag,self.lon0,self.year,self.day)

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

        if self.MetaData['fmtout']:
            # Get Location
            glat = np.loadtxt(path+'glatf.dat')
            glon = np.loadtxt(path+'glonf.dat')
            zalt = np.loadtxt(path+'zaltf.dat')

            # Get plasma values
            deni = np.loadtxt(path+'denif.dat')
            vsi = np.loadtxt(path+'vsif.dat')
            ti = np.loadtxt(path+'tif.dat')
            te = np.loadtxt(path+'tef.dat')
        else:
            # Get Location
            f = open(path+'glatu.dat','rb')
            glat = np.fromfile(f, dtype='float32')[1:-1]
            f.close()

            f = open(path+'glonu.dat','rb')
            glon = np.fromfile(f, dtype='float32')[1:-1]
            f.close()

            f = open(path+'zaltu.dat','rb')
            zalt = np.fromfile(f, dtype='float32')[1:-1]
            f.close()

            # Get plasma values
            f = open(path+'deniu.dat','rb')
            temp = np.fromfile(f, dtype='float32')
            f.close()
            deni = temp.reshape((nz*nf*ni+2),nt,order='F')[1:-1,:]

            f = open(path+'vsiu.dat','rb')
            temp = np.fromfile(f, dtype='float32')
            f.close()
            vsi = temp.reshape((nz*nf*ni+2),nt,order='F')[1:-1,:]

            f = open(path+'tiu.dat','rb')
            temp = np.fromfile(f, dtype='float32')
            f.close()
            ti = temp.reshape((nz*nf*ni+2),nt,order='F')[1:-1,:]

            f = open(path+'teu.dat','rb')
            temp = np.fromfile(f, dtype='float32')
            f.close()
            te = temp.reshape((nz*nf+2),nt,order='F')[1:-1,:]

        self.glat = np.reshape(glat,(nz,nf),order="F")
        self.glon = np.reshape(glon,(nz,nf),order="F")
        self.zalt = np.reshape(zalt,(nz,nf),order="F")
        self.deni = np.reshape(deni,(nz,nf,ni,nt),order="F")
        self.vsi = np.reshape(vsi,(nz,nf,ni,nt),order="F")
        self.ti = np.reshape(ti,(nz,nf,ni,nt),order="F")
        self.te = np.reshape(te,(nz,nf,nt),order="F")
        del glat, glon, zalt, deni, vsi, ti, te

    def _generate_metadata(self,namelist):
        """ Reads the namelist and generates MetaData based on Parameters
        """

        import re

        self.MetaData['fmtout'] = ('.true.' in namelist[1])

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

    def check_standard_model(self, model_type="all"):
        """ Checks for standard atmospheric inputs

        parameters
        -----------
        model_type : (str)
            Limit check to certain models (default='all')

        Returns
        --------
        mod_keys : (list)
            List of modified keyword for self.MetaData, empty
            if no modifications were made
        """
        mod_keys = list()
        meta_keys = self.MetaData.keys()

        for mkey in meta_keys:
            if mkey.find('Multiplier') > 0:
                if self.MetaData[mkey] != 1:
                    mod_keys.append(mkey)

        return mod_keys
