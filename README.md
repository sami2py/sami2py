# sami2py: sami2py is another model of the ionosphere python style

[![Build Status](https://travis-ci.com/jklenzing/sami2py.svg?token=2cVxNXGyxT253xN5hYXp&branch=master)](https://travis-ci.org/jklenzing/sami2py)

# Overview

Sami2py is a python module that runs the sami2 model, as well as archives, loads and plots the resulting modeled values.

# Installation

First, checkout the repository:

```
    $ git clone https://github.com/jklenzing/sami2py.git;
```

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

```
    $ cd sami2py/
    $ python setup.py install
```

Additionally, you must make and install the fortran executables.

```
  $ cd sami2py/fortran
  $ make clean
  $ make compile
```

Now you can run the sami2 executable (sami2low.x) from anywhere.

# Example

In iPython, run:

```
  $ import sami2py

  $ day = 210
  $ year = 2012
  $ lon = 0

  $ sami2py.run_model(day=day,year=year,lon=lon)
```
Note that the sami2 model runs for 24 hours to clear transients, then begins to output data.

Now load the resultant data:

```
  $ S = sami2py.model(tag='test',lon=lon,year=year,day=day)

```

Plotting features coming soon...
