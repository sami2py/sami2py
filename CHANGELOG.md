# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.0] - TBD
- API changes
  - Store loaded data in xarray object
  - Use consistent keyword order in run_model and Model
  - xarray (<0.12) and pandas (<0.25) are now required packages.  Upper limits enforced for Travis CI testing.  Will remove in future once python 2.7 is deprecated and things settle out.
  - Model.plot_lat_alt() now returns the figure object
- Non-breaking changes
  - Output version / short hash for each model run
  - Move package metadata to setup.cfg
  - Auto-build fortran executables as part of setup.py
  - Added CHANGELOG.md
  - Switched to pytest for unit testing
  - Removes python 3.4 testing from Travis
  - Adds manual install of pandas / xarray to Travis workflow to fix setup
  - Add deprecation warning to plot_alt_lat

## [0.1.2] - 2019-07-02
- Patch to fix loading of unformatted output files.

## [0.1.1] - 2019-05-20
- Patch to documentation.

## [0.1.0] - 2019-05-17
- Initial release
