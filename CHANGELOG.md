# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.3.0] - 2022-08-31
* Maintenance
  * Implement flake8-docstring and hacking packages to lint code
  * Update docstring standards
  * Updated NEP 29 compliance in CI tests
  * Add CI tests for Mac Os X
  * Fixed a bug in Windows CI environment in usage of mingw-64
  * Remove deprecated plotting functions (moved to `sami2py_vis`)

## [0.2.5] - 2021-10-19
* Add basic metadata for ExB Fourier Coefficients
* Fixed a bug with xarray 0.19.0
* Update NEP-29 minimum version to numpy 1.18

## [0.2.4] - 2021-06-23
* Updated the unit attributes to consistently use 'hrs' instead of 'hours'.

## [0.2.3] - 2021-06-16
* Updated Variable and dataset attributes for netcdf export
  * Updated netcdf file in test_data
* Added default drift fourier coefficient array accessible from run model
* Using minimum test version of numpy in accordance with NEP 29
* Removed the sami2py-1.00.namelist and version.txt files from run_name
* Bug Fix
  * Pull version info from a single location
  * fixed a bug in compiling readthedocs
  * added default exb file: setup.py generates a exb.inp file
  * Use integer for longitude in directory structure
  * Improved windows compatibility
* Added ability to input custom ExB Drifts as a Fourier Series
  * return_fourier function in utils.py
  * plot_exb function in _core_class.py
  * Testing return_fourier function in test_utils.py
* Added ability to generate Fourier Coefficients from time series ExB drifts
  * fourier_fit function in utils.py
  * private __make_fourier function in utils.py
  * Testing fourier_fit function in test_utils.py
  * scipy dependency added
* Added default exb file: setup.py generates a exb.inp file
* Added deprecation warnings to plotting functions
* Adjusted output step size in model to more closely match desired cadence
* Migrated CI testing to Github Actions
* Documentation
  * Added zenodo integration
  * Improved maintainability of conf.py

## [0.2.2] - 2020-07-17
* Added simple port of core data to netcdf file
* Increased unformatted test data to 6 time steps
* Add namelist and exb.inp in fortran dir to gitignore
* Documentation Changes
  * Primary branch now `main`
  * Improved discussion of install / usage for first time users

## [0.2.1] - 2020-04-13
* Documentation Changes
  * Updates to docstrings
* Non-breaking changes
  * Stylistic updates to conform to flake8
* Testing changes
  * Update Travis CI to use flake8 tests as part of CI
  * Update Travis CI to test for numpy versions as in NEP 029

## [0.2.0] - 2019-12-17
* API changes
  * Store loaded data in xarray object
  * Use consistent keyword order in run_model and Model
  * xarray and pandas are now required packages.  
  * Model.plot_lat_alt() now returns the figure object
  * Pass dimensions into `get_unformatted_data` as tuple
* Non-breaking changes
  * Output version / short hash for each model run
  * Move package metadata to setup.cfg
  * Auto-build fortran executables as part of setup.py
  * Added CHANGELOG.md
  * Switched to pytest for unit testing
  * Removes python 3.4 testing from Travis
  * Adds manual install of pandas / xarray to Travis workflow to fix setup
  * Streamline `_archive_model` to generate filelist and move/copy files
  * Reduce duplication in `_generate_metadata` with `find_int` and `find_float` functions
  * Streamline `Model.check_standard_model` and add check for Fourier Coefficients
  * Add deprecation warning to plot_alt_lat
  * Adds Appveyor support
* Bugs fixed
  * Fix for windows directories
  * New directory structure in .sami2py adds virtual environment flexibility
  * Fixed a bug in feedback when there are no files to move

## [0.1.2] - 2019-07-02
* Patch to fix loading of unformatted output files.

## [0.1.1] - 2019-05-20
* Patch to documentation.

## [0.1.0] - 2019-05-17
* Initial release
