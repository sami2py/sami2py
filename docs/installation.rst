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
