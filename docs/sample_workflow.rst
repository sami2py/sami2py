Sample Workflow
===============

In iPython, run:

.. code:: python

  import sami2py

sami2py will remind you to set the top level directory that will hold the model output.

.. code:: python

  sami2py.utils.set_archive_dir(path=path)

sami2py will raise an error if this is not done before trying to run the model.

.. code:: python

  sami2py.run_model(tag='run_name', lon=0, year=2012, day=210)

Note that the sami2 model runs for 24 hours to clear transients, then begins to output data.

Now load the resultant data:

.. code:: python

  ModelRun = sami2py.Model(tag='run_name', lon=0, year=2012, day=210)

The data is stored as `ModelRun.data`, which is an `xarray.Dataset`.  Information about the run is stored as 'ModelRun.MetaData', which is a human-readable dictionary of the namelist.

The MetaData can be accessed directly via the dictionary, or through the __repr__.  Typing

.. code:: python

  ModelRun

yields

.. code:: python

  Model Run Name = test
  Day 256, 1999
  Longitude = 256 deg
  2 time steps from  0.1 to  0.1 UT
  Ions Used: H+, O+, NO+, O2+, He+, N2+

  Solar Activity
  --------------
  F10.7: 120.0 sfu
  F10.7A: 120.0 sfu
  ap: 0

  Component Models Used
  ---------------------
  Neutral Atmosphere: NRLMSISe-2000
  Winds: HWM-14
  Photoproduction: EUVAC
  ExB Drifts: Fejer-Scherliess

  No modifications to empirical models

Full description coming soon
