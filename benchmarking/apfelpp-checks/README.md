In this folder can be found the datafiles used to make comparisons
with Valerio Bertone (Feb 2023) against his implementation of the
strucutre functions in APPFEL++.

The code used to produce the datafiles can be found in
structure_functions_valerio_checks.f90.

We use the LHA toy pdf initialised at sqrt(2) GeV with αS(sqrt(2)) =
0.35 and with a variable flavour scheme.

In general we find that, up to NNLO, we have agreement at the or
better than permille level except close to zeroes. At N3LO the
discrepancies are slightly larger, but in general we do also agree at
the permille level. At this point we think it is due to numerical
limitations in the computaton of the rather wiggly functions. Further
tests are needed to see if we can get rid of these.

