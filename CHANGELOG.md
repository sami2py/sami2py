# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.1] - 2020-04-13
- Documentation Changes
  - Updates to docstrings
- Non-breaking changes
  - Stylistic updates to conform to flake8
- Testing changes
  - Update Travis CI to use flake8 tests as part of CI
  - Update Travis CI to test for numpy versions as in NEP 029

## [0.2.0] - 2019-12-17
- API changes
  - Store loaded data in xarray object
  - Use consistent keyword order in run_model and Model
  - xarray and pandas are now required packages.  
  - Model.plot_lat_alt() now returns the figure object
  - Pass dimensions into `get_unformatted_data` as tuple
- Non-breaking changes
  - Output version / short hash for each model run
  - Move package metadata to setup.cfg
  - Auto-build fortran executables as part of setup.py
  - Added CHANGELOG.md
  - Switched to pytest for unit testing
  - Removes python 3.4 testing from Travis
  - Adds manual install of pandas / xarray to Travis workflow to fix setup
  - Streamline `_archive_model` to generate filelist and move/copy files
  - Reduce duplication in `_generate_metadata` with `find_int` and `find_float` functions
  - Streamline `Model.check_standard_model` and add check for Fourier Coefficients
  - Add deprecation warning to plot_alt_lat
  - Adds Appveyor support
- Bugs fixed
  - Fix for windows directories
  - New directory structure in .sami2py adds virtual environment flexibility
  - Fixed a bug in feedback when there are no files to move

## [0.1.2] - 2019-07-02
- Patch to fix loading of unformatted output files.

## [0.1.1] - 2019-05-20
- Patch to documentation.

## [0.1.0] - 2019-05-17
- Initial release
