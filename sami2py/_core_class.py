#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK & JH
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Wrapper for running sami2 model.

Classes
-------
Model
    Loads, reshapes, and holds SAMI2 output for a given model run
    specified by the user.

"""

import numpy as np
from os import path
import re

import xarray as xr

import sami2py


class Model(object):
    """Python object to handle SAMI2 model output data."""

    def __init__(self, tag, lon, year, day, outn=False, test=False):
        """Load a previously run sami2 model and sorts into array shapes.

        Parameters
        ----------
        tag : str
            name of run (top-level directory)
        lon : float
            longitude reference
        year : int
            year
        day : int
            day of year from Jan 1
        outn : bool
            if true : look for neutral density and wind files
            if false :  only look for default sami2 output
        test : bool
            if true : use test model output
            if false : look for user made model output

        Returns
        -------
        self : model class object containing OCB file data

        Attributes
        ----------
        ut : 1D ndarray
            Universal Time (hrs)
        slt : 1D ndarray
            Solar Local Time (hr)
        glat : 2D ndarray
            Geographic Latitude (deg)
        glon : 2D ndarray
            Geographic Longitude (deg)
        zalt : 2D ndarray
            Altitude (km)
        deni : 4D ndarray
            Ion density for each species (cm^-3)
        vsi : 4D ndarray
            Ion Velocity for each species (cm/s)
        ti : 4D ndarray
            Ion Temperature for each species (K)
        te : 3D ndarray
            Electron Temperature (K)

        Examples
        --------
        To load a previous model run:
            ModelRun = sami2py.Model(tag='run_name', lon=0, year=2012, day=210)

        """

        self.tag = tag
        self.lon0 = lon
        self.year = year
        self.day = day
        self.outn = outn
        self.test = test

        self._load_model()

    def __repr__(self):
        """Make a printable representation of a Model object.

        Returns
        -------
        out : str
            string containing a printable representation of a Model object

        Examples
        --------
        Load the model
        ::

            ModelRun = sami2py.Model(tag='run_name', lon=0, year=2012, day=210)

        Check the model information
        ::

            ModelRun

        """

        out = ['']
        out.append('Model Run Name = ' + self.tag)
        out.append(('Day {day:03d}, {year:4d}').format(day=self.day,
                                                       year=self.year))
        out.append(('Longitude = {:5.1f} deg').format(self.lon0))
        temp_str = '{N:d} time steps from {t0:4.1f} to {tf:4.1f} UT'
        out.append(temp_str.format(N=len(self.ut),
                                   t0=min(self.ut),
                                   tf=max(self.ut)))
        out.append('Ions Used: ' + self.MetaData['Ions Used'])

        out.append('\nSolar Activity')
        out.append('--------------')
        temp_str = 'F10.7: {f107:5.1f} sfu'
        out.append(temp_str.format(f107=self.MetaData['F10.7']))
        temp_str = 'F10.7A: {f107a:5.1f} sfu'
        out.append(temp_str.format(f107a=self.MetaData['F10.7A']))
        out.append(('ap: {ap:d}').format(ap=self.MetaData['ap']))

        out.append('\nComponent Models Used')
        out.append('---------------------')
        out.append(' '.join(('Neutral Atmosphere:',
                             self.MetaData['Neutral Atmosphere Model'])))
        out.append('Winds: ' + self.MetaData['Wind Model'])
        out.append('Photoproduction: ' + self.MetaData['EUV Model'])
        out.append('ExB Drifts: ' + self.MetaData['ExB model'])

        mod_keys = self.check_standard_model()
        if mod_keys:
            out.append('\nMultipliers used')
            out.append('----------------')
            for mkey in mod_keys:
                out.append(('{s}: {f}').format(s=mkey,
                                               f=self.MetaData[mkey]))
        else:
            out.append('\nNo modifications to empirical models')

        return '\n'.join(out)

    def _calculate_slt(self):
        """Calculate Solar Local Time for reference point of model."""

        local_time = np.mod((self.ut * 60 + self.lon0 * 4), 1440) / 60.0
        mean_anomaly = 2 * np.pi * self.day / 365.242
        delta_t = (-7.657 * np.sin(mean_anomaly)
                   + 9.862 * np.sin(2 * mean_anomaly + 3.599))
        self.slt = local_time - delta_t / 60.0

    def _load_model(self):
        """Load model results."""

        nf = 98
        nz = 101
        ni = 7

        model_path = sami2py.utils.generate_path(self.tag, self.lon0, self.year,
                                                 self.day, self.test)

        # Get NameList
        namelist_file = open(path.join(model_path, 'sami2py-1.00.namelist'))
        self.namelist = namelist_file.readlines()
        namelist_file.close()

        self.MetaData = dict()
        self._generate_metadata(self.namelist, model_path)

        # Get times
        time = np.loadtxt(path.join(model_path, 'time.dat'))
        self.ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600

        self._calculate_slt()
        nt = len(self.ut)

        if self.MetaData['fmtout']:
            # Get Location
            glat = np.loadtxt(path.join(model_path, 'glatf.dat'))
            glon = np.loadtxt(path.join(model_path, 'glonf.dat'))
            zalt = np.loadtxt(path.join(model_path, 'zaltf.dat'))

            # Get plasma values
            deni = np.loadtxt(path.join(model_path, 'denif.dat'))
            vsi = np.loadtxt(path.join(model_path, 'vsif.dat'))
            ti = np.loadtxt(path.join(model_path, 'tif.dat'))
            te = np.loadtxt(path.join(model_path, 'tef.dat'))

            # get neutral values
            if self.outn:
                denn = np.loadtxt(path.join(model_path, 'dennf.dat'))
                u4 = np.loadtxt(path.join(model_path, 'u4f.dat'))
        else:
            # Get Location
            glat = sami2py.utils.get_unformatted_data(model_path, 'glat')
            glon = sami2py.utils.get_unformatted_data(model_path, 'glon')
            zalt = sami2py.utils.get_unformatted_data(model_path, 'zalt')

            # Get plasma values
            dim0 = nz * nf * ni + 2
            dim1 = nt
            deni = sami2py.utils.get_unformatted_data(model_path, 'deni',
                                                      dim=(dim0, dim1),
                                                      reshape=True)
            vsi = sami2py.utils.get_unformatted_data(model_path, 'vsi',
                                                     dim=(dim0, dim1),
                                                     reshape=True)
            ti = sami2py.utils.get_unformatted_data(model_path, 'ti',
                                                    dim=(dim0, dim1),
                                                    reshape=True)

            # Electron Temperatures have only one species
            dim0 = nz * nf + 2
            te = sami2py.utils.get_unformatted_data(model_path, 'te',
                                                    dim=(dim0, dim1),
                                                    reshape=True)
            if self.outn:
                # Multiple neutral species
                dim0 = nz * nf * ni + 2
                denn = sami2py.utils.get_unformatted_data(model_path, 'denn',
                                                          dim=(dim0, dim1),
                                                          reshape=True)
                # Only one wind
                dim0 = nz * nf + 2
                u4 = sami2py.utils.get_unformatted_data(model_path, 'u4',
                                                        dim=(dim0, dim1),
                                                        reshape=True)

        glat = np.reshape(glat, (nz, nf), order="F")
        glon = np.reshape(glon, (nz, nf), order="F")
        zalt = np.reshape(zalt, (nz, nf), order="F")
        deni = np.reshape(deni, (nz, nf, ni, nt), order="F")
        vsi = np.reshape(vsi, (nz, nf, ni, nt), order="F")
        ti = np.reshape(ti, (nz, nf, ni, nt), order="F")
        te = np.reshape(te, (nz, nf, nt), order="F")
        self.data = xr.Dataset({'deni': (['z', 'f', 'ion', 'ut'], deni,
                                         {'units': 'N/cc',
                                          'long_name': 'ion density'}),
                                'vsi': (['z', 'f', 'ion', 'ut'], vsi,
                                        {'units': 'cm/s',
                                         'long_name': 'ion velocity'}),
                                'ti': (['z', 'f', 'ion', 'ut'], ti,
                                       {'units': 'K',
                                        'long_name': 'ion temperature'}),
                                'te': (['z', 'f', 'ut'], te,
                                       {'units': 'K',
                                        'long_name': 'electron temperature'}),
                                'slt': (['ut'], self.slt,
                                        {'units': 'hrs',
                                         'long_name': 'solar local time'})},
                               coords={'glat': (['z', 'f'], glat,
                                                {'units': 'deg',
                                                 'long_name': 'Geo Latitude',
                                                 'value_min': -90.0,
                                                 'value_max': 90.0}),
                                       'glon': (['z', 'f'], glon,
                                                {'units': 'deg',
                                                 'long_name': 'Geo Longitude',
                                                 'value_min': 0.0,
                                                 'value_max': 360.0}),
                                       'zalt': (['z', 'f'], zalt,
                                                {'units': 'km',
                                                 'long_name': 'Altitude'}),
                                       'ut': (['ut'], self.ut,
                                              {'units': 'hrs',
                                               'long_name': 'Universal Time'})})
        if self.outn:
            denn = np.reshape(denn, (nz, nf, ni, nt), order="F")
            self.data['denn'] = (('z', 'f', 'ion', 'ut'), denn,
                                 {'units': 'N/cc',
                                  'long_name': 'neutral density'})
            u4 = np.reshape(u4, (nz, nf, nt), order="F")
            self.data['u4'] = (('z', 'f', 'ut'), u4,
                               {'units': 'm/s',
                                'long_name': 'neutral wind velocity'})

        if self.MetaData['ExB model'] == 'Fourier Series':
            exb = sami2py.utils.return_fourier(self.data['slt'],
                                               self.MetaData['Fourier Coeffs'])
            self.data['exb'] = (('ut'), exb.data,
                                {'units': 'm/s',
                                 'long_name': 'ExB Foureir Coefficients'})

    def _generate_metadata(self, namelist, model_path):
        """Read the namelist and generates MetaData based on Parameters.

        Parameters
        -----------
        namelist : list
            variable namelist from SAMI2 model

        """

        def find_float(name, ind):
            """Search for regular expression float vals."""

            return float(re.findall(r"\d*\.\d+|\d+", name)[ind])

        def find_int(name, ind):
            """Search for regular expression int vals."""

            return int(re.findall(r"\d+", name)[ind])

        self.MetaData['Model Run Name'] = self.tag
        self.MetaData['Day'] = '{day:03d}, {year:4d}'.format(day=self.day,
                                                              year=self.year)
        self.MetaData['Longitude'] = '{:.1f}'.format(self.lon0)

        self.MetaData['fmtout'] = ('.true.' in namelist[1])

        self.MetaData['F10.7A'] = find_float(namelist[14], 0)
        self.MetaData['F10.7'] = find_float(namelist[15], 2)
        self.MetaData['ap'] = find_int(namelist[16], 0)

        self.MetaData['Neutral Atmosphere Model'] = 'NRLMSISe-2000'
        self.MetaData['EUV Model'] = 'EUVAC'

        # Ions Used
        nion1 = wind_model = find_int(namelist[20], 1) - 1
        nion2 = wind_model = find_int(namelist[21], 1) - 1
        ions = ['H+', 'O+', 'NO+', 'O2+', 'He+', 'N2+', 'N+']
        self.MetaData['Ions Used'] = ', '.join(ions[nion1:nion2])

        # Multipliers
        neutral_scalars = re.findall(r"\d*\.\d+|\d+", namelist[28])
        self.MetaData['H Multiplier'] = float(neutral_scalars[0])
        self.MetaData['O Multiplier'] = float(neutral_scalars[1])
        self.MetaData['NO Multiplier'] = float(neutral_scalars[2])
        self.MetaData['O2 Multiplier'] = float(neutral_scalars[3])
        self.MetaData['He Multiplier'] = float(neutral_scalars[4])
        self.MetaData['N2 Multiplier'] = float(neutral_scalars[5])
        self.MetaData['N Multiplier'] = float(neutral_scalars[6])
        self.MetaData['T_exo Multiplier'] = find_float(namelist[33], 0)
        self.MetaData['T_n Multiplier'] = find_float(namelist[29], 0)
        self.MetaData['EUV Multiplier'] = find_float(namelist[34], 0)
        self.MetaData['ExB Drift Multiplier'] = find_float(namelist[24], 1)
        self.MetaData['Wind Multiplier'] = find_float(namelist[23], 1)

        if '.true.' in namelist[10]:
            self.MetaData['ExB model'] = 'Fejer-Scherliess'
            self.MetaData['Fourier Coeffs'] = np.loadtxt(path.join(model_path,
                                                                   'exb.inp'))
        else:
            self.MetaData['ExB model'] = 'Fourier Series'
            self.MetaData['Fourier Coeffs'] = np.loadtxt(path.join(model_path,
                                                                   'exb.inp'))

        wind_model = find_int(namelist[35], 0)
        self.MetaData['Wind Model'] = ('HWM-{:02d}').format(wind_model)

        # Model Geometry
        self.MetaData['rmin'] = find_float(namelist[11], 0)
        self.MetaData['rmax'] = find_float(namelist[12], 0)
        self.MetaData['gams'] = find_int(namelist[26], 0)
        self.MetaData['gamp'] = find_int(namelist[27], 0)
        self.MetaData['altmin'] = find_float(namelist[13], 0)

        # Model runtime
        self.MetaData['dthr'] = find_float(namelist[5], 0)
        self.MetaData['hrinit'] = find_float(namelist[22], 0)
        self.MetaData['hrpr'] = find_float(namelist[6], 0)
        self.MetaData['hrmax'] = find_float(namelist[3], 0)
        self.MetaData['dt0'] = find_float(namelist[4], 0)
        self.MetaData['maxstep'] = find_int(namelist[2], 0)
        self.MetaData['denmin'] = find_float(namelist[30], 0)

    def check_standard_model(self, model_type="all"):
        """Check for standard atmospheric inputs.

        Parameters
        ----------
        model_type : str
            Limit check to certain models (default='all')
            Not currently implemented

        Returns
        -------
        mod_keys : list
            List of modified keyword for self.MetaData, empty
            if no modifications were made

        Examples
        --------
        Load the model
        ::

            ModelRun = sami2py.Model(tag='run_name', lon=0, year=2012, day=210)

        Check the model information for changes to the standard inputs
        ::

            ModelRun.check_standard_model()

        """

        mod_keys = list()
        meta_keys = list(self.MetaData.keys())

        # See if Fourier Coefficients are used
        mkey = 'Fourier Coeffs'
        if mkey in meta_keys:
            mod_keys.append(mkey)
            # Since this item is an array, remove for multiplier check
            meta_keys.remove(mkey)

        # Check for scalar multipliers
        for mkey in meta_keys:
            if ((mkey.find('Multiplier') > 0) & (self.MetaData[mkey] != 1)):
                mod_keys.append(mkey)

        return mod_keys

    def to_netcdf(self, path=''):
        """Save core data as a netcdf file."""

        if path == '':
            path = 'sami2py_output.nc'
        attrs = {}
        keys = self.MetaData.keys()
        for key in keys:
            new_key = key.replace(' ', '_').replace('.', '_')
            if key == 'Fourier Coeffs':
                terms = [", ".join(item) for item
                         in self.MetaData[key].astype(str)]
                coeffs = '; '.join(terms)
                attrs[new_key] = coeffs
            else:
                attrs[new_key] = self.MetaData[key]

        attrs['fmtout'] = str(attrs['fmtout'])

        self.data.attrs = attrs
        self.data.to_netcdf(path=path, format='NETCDF4')
