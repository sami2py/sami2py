
.. |br| raw:: html

    <br>

Introduction
============

Sami2py is a python module that runs the SAMI2 model, as well as archives, loads and plots the resulting modeled values. SAMI2 is a model developed by the Naval Research Laboratory to simulate the motions of plasma in a 2D ionospheric environment along a dipole magnetic field [Huba et al, 2000].  SAMI2 solves for the chemical and dynamical evolution of seven ion species in this environment (H\ :sup:`+`\, He\ :sup:`+`\, N\ :sup:`+`\, O\ :sup:`+`\, N\ :sub:`2` :sup:`+`\, NO\ :sup:`+`\, and O\ :sub:`2` :sup:`+`\).

The implementation used here includes several added options to the original release of SAMI2.  A full list is included in :ref:`modifications`, but several of these include:
 - The ability to scale the neutral atmosphere in which the ions form through direct modification of the exospheric neutral temperature for extreme solar minimum conditions, as discussed by Emmert et al [2010].
 - The ability to switch between HWM93, HWM07, and HWM14 as a user option.
 This implementation is based on the matlab version used in Klenzing et al [2013].

The open-source fortran version of SAMI2 is found at https://www.nrl.navy.mil/ppd/branches/6790/sami2


How to Cite
-----------

When referring to this software package, please cite the original paper by Huba et al [2000] https://doi.org/10.1029/2000JA000035 as well as the package by Klenzing and Smith [2019] https://doi.org/10.5281/zenodo.2875800.

Additionally, please include the following text in the acknowledgements: "This work uses the SAMI2 ionosphere model written and developed by the Naval Research Laboratory."


References
----------

 - Huba, J.D., G. Joyce, and J.A. Fedder, Sami2 is Another Model of the Ionosphere (SAMI2): A new low‐latitude ionosphere model, *J. Geophys. Res.*, 105, Pages 23035-23053, https://doi.org/10.1029/2000JA000035, 2000.
 - Emmert, J.T., J.L. Lean, and J.M. Picone, Record‐low thermospheric density during the 2008 solar minimum, *Geophys. Res. Lett.*, 37, https://doi.org/10.1029/2010GL043671, 2010.
 - Klenzing, J., A. G. Burrell, R. A. Heelis, J. D. Huba, R. Pfaff, and F. Simões, Exploring the role of ionospheric drivers during the extreme solar minimum of 2008, *Ann. Geophys.*, 31, 2147-2156, https://doi.org/10.5194/angeo-31-2147-2013, 2013.
