Installation
============

First, checkout the repository::

  git clone git@gitlab.com:jklenzing/sami2py.git

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

  cd sami2py/
  python setup.py install

Additionally, you must make and install the fortran executables.

  make -C sami2py/fortran compile

Now you can run the sami2 executable (sami2py.x) from anywhere.
