Sample Workflow
===============

In iPython, run:

```
  import sami2py
  sami2py.run_model(year=2012, day=210, lon=0, tag='test')
```
Note that the sami2 model runs for 24 hours to clear transients, then begins to output data.

Now load the resultant data:

```
  ModelRun = sami2py.model(tag='test', lon=0, year=2012, day=210)
```
