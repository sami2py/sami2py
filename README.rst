sami2py: sami2py is another model of the ionosphere python style
================================================================

.. list-table::
    :stub-columns: 1

    * - docs
      - | |rtd| |doi|
    * - tests
      - | |pytest|
        | |coveralls| |codecov|
        | |codeclimate|

.. |rtd| image:: https://readthedocs.org/projects/sami2py/badge/?version=latest
    :target: http://sami2py.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |pytest| image:: https://github.com/sami2py/sami2py/actions/workflows/main.yml/badge.svg
    :target: https://github.com/sami2py/sami2py/actions/workflows/main.yml
    :alt: Pytest with Flake8

.. |coveralls| image:: https://coveralls.io/repos/github/sami2py/sami2py/badge.svg?branch=main
    :target: https://coveralls.io/github/sami2py/sami2py?branch=main
    :alt: Coverage Status

.. |codecov| image:: https://codecov.io/gh/sami2py/sami2py/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/sami2py/sami2py
    :alt: Coverage Status

.. |codeclimate| image:: https://api.codeclimate.com/v1/badges/306cb2d5c709707f7b64/maintainability
   :target: https://codeclimate.com/github/sami2py/sami2py
   :alt: CodeClimate Quality Status

.. |doi| image:: https://zenodo.org/badge/167871330.svg
  :target: https://zenodo.org/badge/latestdoi/167871330


Overview
--------

Sami2py is a python module that runs the SAMI2 model, as well as archives and
loads the resulting modeled values. SAMI2 is a model developed by the Naval
Research Laboratory to simulate the motions of plasma in a 2D ionospheric
environment along a dipole magnetic field [Huba et al, 2000].  SAMI2 solves for
the chemical and dynamical evolution of seven ion species in this environment
(H\ :sup:`+`\, He\ :sup:`+`\, N\ :sup:`+`\, O\ :sup:`+`\, N\ :sub:`2`\ :sup:`+`\,
NO\ :sup:`+`\, and O\ :sub:`2`\ :sup:`+`\).

The implementation used here includes several added options to the original
release of SAMI2.  A full list is included in
https://sami2py.readthedocs.io/en/latest/modifications.html, but several of
these include:
 - The ability to scale the neutral atmosphere in which the ions form through direct modification of the exospheric neutral temperature for extreme solar minimum conditions, as discussed by Emmert et al [2010].
 - The ability to input custom ExB drifts as a Fourier series.

 This implementation is based on the matlab version used in Klenzing et al [2013].


Installation
------------

First, sami2py depends on fortran. Information on installing the GFortran compiler can be found here:

https://gcc.gnu.org/wiki/GFortranBinaries

Next, checkout the repository:

.. code-block:: console

  git clone https://github.com/sami2py/sami2py.git

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

.. code-block:: console

  cd sami2py
  python setup.py install

If something has gone wrong, you may be prompted to manually install the fortran executables.

.. code-block:: console

  make -C sami2py/fortran compile

or, on windows,

.. code-block:: console

  make -C sami2py\fortran compile

Note that you will need a fortran compiler (gfortran is the default setup) and make installed on your system.  For more information, please refer to the
`documentation <https://sami2py.readthedocs.io/en/latest/installation.html#fortran-compilers>`_.


Example
-------

In iPython, run:

.. code-block:: python

  import sami2py

If this is your first import of sami2py, it will remind you to set the top level directory that will hold the model output.  This should be a string containing the path to the directory you want to store the data in, such as ``path='/Users/me/data/sami2py'`` or ``path='C:\home\data'``.  This should be outside the main code directory, so model output files are not confused with model inputs or source code.  If you are using Git, it will also ensure that Git does not try to store your local code runs within the repository.

.. code-block:: python

  sami2py.utils.set_archive_dir(path=path)

sami2py will raise an error if this is not done before trying to run the model.

.. code-block:: python

  sami2py.run_model(tag='run_name', lon=0, year=2012, day=210)

Note that the sami2 model runs for 24 hours to clear transients, then begins to output data.

Now load the resultant data:

.. code-block:: python

  ModelRun = sami2py.Model(tag='run_name', lon=0, year=2012, day=210)

How to Cite
-----------
When referring to this software package, please cite the original paper by Huba et al [2000] https://doi.org/10.1029/2000JA000035 as well as the package by Klenzing et al [2019] https://doi.org/10.5281/zenodo.2875799. Note that this doi will always point to the latest version of the code.  The specific version doi can be found at the top of this page.

Additionally, please include the following text in the acknowledgements: "This
work uses the SAMI2 ionosphere model written and developed by the Naval Research Laboratory."

References
----------
- Huba, J.D., G. Joyce, and J.A. Fedder, Sami2 is Another Model of the Ionosphere (SAMI2): A new low‐latitude ionosphere model, *J. Geophys. Res.*, 105, Pages 23035-23053, https://doi.org/10.1029/2000JA000035, 2000.
- Emmert, J.T., J.L. Lean, and J.M. Picone, Record‐low thermospheric density during the 2008 solar minimum, *Geophys. Res. Lett.*, 37, https://doi.org/10.1029/2010GL043671, 2010.
- Klenzing, J., A. G. Burrell, R. A. Heelis, J. D. Huba, R. Pfaff, and F. Simões, Exploring the role of ionospheric drivers during the extreme solar minimum of 2008, *Ann. Geophys.*, 31, 2147-2156, https://doi.org/10.5194/angeo-31-2147-2013, 2013.
