Installation
============

First, checkout the repository:

.. code-block:: console

  git clone https://github.com/sami2py/sami2py.git

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

.. code-block:: console

  cd sami2py/
  python setup.py install

If something has gone wrong, you may be prompted to manually install the fortran executables.

.. code-block:: console

  make -C sami2py/fortran compile


or, on windows,

.. code-block:: console

  make -C sami2py\fortran compile


Fortran Compilers
-----------------

By default, sami2py uses gfortran and make to compile the fortran executables.
If you don't have a fortran compiler, gfortran is included as part of the latest
gcc package.  You can get this from several locations.

For Mac OS X, you can install gcc through package managers such as `brew <https://brew.sh/>`_.

For windows, make sure that mingw-64 is installed.  This may need to be reinstalled to make sure links work properly.  See discussion at https://www.scivision.dev/cmake-install-windows

Additionally, make is required to compile the code.  You can get make through pip.

.. code-block:: console

  pip install make

If you have a different compiler, you can modify the first line of the fortran
Makefile accordingly by using "gf" to point to your compiler of choice.  Note
that several options are included to ensure the compile is successful.

.. code-block::

  gf = gfortran -fno-range-check -fno-automatic -ffixed-line-length-none
