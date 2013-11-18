stress_corr
===========

Computes the not-diagonal stress correlation function from output files of PIMAIM (xyxzyzstress.out and xxyyzzstress.out)


This program has been written by Mathieu Salanne, UPMC, PECSA, France.


stress_corr outputs 5 files:
- scorr.out_2zz-xx-yy
- scorr.out_xx-yy
- scorr.out_xy
- scorr.out_xz
- scorr.out_yz

In order to compute the viscosity, one has to make the average of the 5 (yes, five) files, and integrate this average.

All 5 stress correlation functions should tend toward 0 after a short time.
Be carefull to print the stress tensor frequently enough.

Conversion is:
convEta=((auEta*volA)/T)*6.260534902d7
