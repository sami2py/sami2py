[![Build Status](https://travis-ci.com/jklenzing/sami2py.svg?branch=develop)](https://travis-ci.com/jklenzing/sami2py)

# Overview

Sami2py is a python module that runs the SAMI2 model, as well as archives, loads and plots the resulting modeled values. SAMI2 is a model developed by the Naval Research Laboratory to simulate the motions of plasma in a 2D ionospheric environment along a dipole magnetic field [Huba et al, 2000].  SAMI2 solves for the chemical and dynamical evolution of seven ion species in this environment (H<sup>+</sup>, He^+, N^+, O^+, N_2^+, NO^+, and O_2^+).

The implementation used here includes the ability to scale the neutral atmosphere in which the ions form through direct modification of the exospheric neutral temperature for extreme solar minimum conditions, as discussed by Emmert et al [20??].  This implementation is based on the version used in Klenzing et al [2013].

# Installation

First, checkout the repository:

```
    $ git clone git@gitlab.com:jklenzing/sami2py.git;
```

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

```
    $ cd sami2py/
    $ python setup.py install
```

Additionally, you must make and install the fortran executables.

```
  $ make -C sami2py/fortran compile
```

Now you can run the sami2 executable (sami2low.x) from anywhere.

# Example

In iPython, run:

```
  $ import sami2py
  $ sami2py.run_model(day=210, year=2012, lon=0)
```
Note that the sami2 model runs for 24 hours to clear transients, then begins to output data.

Now load the resultant data:

```
  $ ModelRun = sami2py.model(tag='test', lon=lon, year=year, day=day)

```

Plotting features coming soon...

# How to Cite
When referring to this software package, please cite the original paper by Huba et al [2000] https://doi.org/10.1029/2000JA000035 as well as the package by Klenzing et al [2018]. ***add doi here***.

Additionally, please include the following text in the acknowledgements: "This
work uses the SAMI2 ionosphere model written and developed by the Naval Research Laboratory."
