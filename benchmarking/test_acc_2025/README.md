# README for the test_acc_2025/ directory

This directory contains the scripts used for the 2025 paper discussing
the v2 release of hoppet. The datafiles for the actual precision and
timing results are to be found at
https://github.com/hoppet-code/2025-prec-and-timing
For some plotting scripts, this is assumed to be in the
which is assumed to be in the `../../../2025-prec-and-timing/`
directory.

## Scripts for generating results

- `run-many.sh buildDir outputDir`: generates the reference runs and a 
  results for a range of grid spacings, and interpolation orders. 

- `run-many-lhapdf.sh buildDir outputDir`: generates lhapdf grids 
  as well as hoppet results for a range of grid spacing. Note that
  to use the lhapdf grids, your LHAPDF path should including directory
  that contains these files (or they should be moved into a standard
  LHAPDF location)

- `run-and-process-lhapdf-grids.sh buildDir refFile outputDir`:
  runs a LHAPDF to sample the different flavours at the same
  x,Q points as used in the hoppet precision tests, to allow
  comparisons of precision (compared to refFile). It also produces 
  timings.

These scripts use
[prec_and_timing.f90](../prec_and_timing.f90), 
[compare2files_v2.f90](../compare2files_v2.f90) and
[lhapdf_timings.cc](../lhapdf_timings.cc),
with the corresponding executables expected to be located in
`buildDir/benchmarking`.





