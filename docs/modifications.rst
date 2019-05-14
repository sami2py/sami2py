.. _modifications:

Modifications from SAMI2-1.00
========================================

NRLMSISe-00
-----------
This version uses the official release of NRLMSISe-00 (with one modification, discussed below). The unmodified version can be found at https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSISE-00/

The exospheric neutral temperature can be directly scaled by the user for extreme solar minimum conditions, as discussed by Emmert et al [2010].  This is modified by the `Tinf_scale` keyword, and is passed through the namelist as Tinf_scl.

Photoionization
---------------
New namelist variable: euv_scl: Scales the resultant EUV spectra to change photoionization rates.

ExB Drifts
----------
Fejer-Scherliess remains the default model for drifts, but the user may now input a series of Fourier coefficients to describe the ExB drifts as a function of local time rather than use the Fejer-Scherliess model.

Horizontal Wind Model
---------------------
Sami2py uses the HWM-14 model by default.
