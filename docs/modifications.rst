.. _modifications:

Modifications from SAMI2-1.00
========================================

- Update to the official release of NRLMSISe-00. https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSISE-00/
- New namelist variable: Tinf_scl: Scales the neutral atmosphere in which the ions form through direct modification of the exospheric neutral temperature for extreme solar minimum conditions, as discussed by Emmert et al [2010].
- New namelist variable: euv_scl: Scales the resultant EUV spectra to change photoionization rates.
- New input file: exb.inp: Inputs a series of Fourier coefficients to describe the ExB drifts as a function of local time rather than use the Fejer-Scherliess model. Fejer-Scherliess remains the default model for drifts.
- New namelist variable: hwm_mod: Selects version of horizontal wind model, including HWM93, HWM07, and HWM14.  HWM14 is set as the default option.
