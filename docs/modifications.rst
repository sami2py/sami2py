.. _modifications:

Modifications from SAMI2-1.00
========================================

NRLMSISe-00
-----------
This version uses the official release of NRLMSISe-00 (with one modification, discussed below). The unmodified version can be found at https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSISE-00/

The exospheric neutral temperature can be directly scaled by the user for extreme solar minimum conditions, as discussed by Emmert et al [2010].  This is modified by the ``Tinf_scale`` keyword, and is passed through the namelist as Tinf_scl.

Photoionization
---------------
The photoionization rates can be modified by scaling the resultant EUV spectra.  Note that this occurs independently of any modifications to the neutral atmosphere through MSIS.  See Klenzing et al [2013] for examples. This is modified by the ``euv_scale`` keyword, and is passed through the namelist as euv_scl.

ExB Drifts
----------
Fejer-Scherliess remains the default model for drifts, but the user may now input a series of Fourier coefficients to describe the ExB drifts as a function of local time rather than use the Fejer-Scherliess model.  The coefficients are specified by the ``ExB_drifts`` keyword, which are passed into sami2 through the new exb.inp file.  The fourier coefficient routine replaces the sine wave triggered when .fejer. = false.

Horizontal Wind Model
---------------------
Sami2py uses the HWM-14 model by default.  Users may specify HWM93 or HWM07 for comparison through the ``hwm_model`` keyword, which is passed through the namelist as hwm_mod.
