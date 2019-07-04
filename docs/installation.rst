Installation
============

First, checkout the repository:

::

  git clone https://github.com/jklenzing/sami2py.git

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

::

  cd sami2py/
  python setup.py install

If something has gone wrong, you may be prompted to manually install the fortran executables.

::

  make -C sami2py/fortran compile
