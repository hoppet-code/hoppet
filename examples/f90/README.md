# Fortran examples with HOPPET

If building with CMake, by default all the examples in this directory will be built.
Assuming CMake used a `build/` directory, the executables will be found
in the `build/example_f90/`.

## Normal QCD evolution

- [tabulation_example.f90](tabulation_example.f90): a basic example to create a table across x and Q, with the full modern Fortran interface

- [tabulation_example_streamlined.f90](tabulation_example_streamlined.f90): similar example, but using the simpler "streamlined" interface

- [tabulation_example_2.f90](tabulation_example_2.f90): a more structured alternative to `tabulation_example.f90`, illustrating how one might organise the different objects

- [tabulation_example_n3lo.f90](tabulation_example_n3lo.f90): similar to `tabulation_example.f90`, with N3LO evolution and internal options for variable-flavour vs. fixed-flavour number schemes

- [tabulation_example_up_and_down.f90](tabulation_example_up_and_down.f90): an example of carrying out evolution up in scale to create a pdf_table and then back down, to demonstrate the code's ability to carry out such operations

## QCD+QED evolution

- [tabulation_example_qed.f90](tabulation_example_qed.f90): an example with QCD+QED evolution in the full modern Fortran interface
- [tabulation_example_qed_streamlined.f90](tabulation_example_qed_streamlined.f90): an example with QCD+QED evolution in the streamlined interface


## Structure-function evaluation

- [structure_functions_example.f90](structure_functions_example.f90): example to illustrate determination of (massless) structure functions

- [structure_functions_example_flavour.f90](structure_functions_example_flavour.f90): example that shows how to access structure functions decomposed by flavour and then recombine them

- [structure_functions_example_n3lo_tests.f90](structure_functions_example_n3lo_tests.f90): writes a series of tables with different combinations of NNLO and N3LO coefficient functions and splitting functions, including scale variation in the coefficient functions.

## Loading a PDF set from LHAPDF

These examples are all in the [with-lhapdf/](with-lhapdf/) directory and are only
compiled if using the CMake build system and if LHAPDF has been found. The executables are placed in the `build/example_f90/` directory.

- [with-lhapdf/lhapdf_to_hoppet.f90](with-lhapdf/lhapdf_to_hoppet.f90):
  an example that transfers a single LHAPDF PDF into the HOPPET pdf_table 
  of the streamlined interface

- [with-lhapdf/lhapdf_to_hoppet_allmembers.f90](with-lhapdf/lhapdf_to_hoppet_allmembers.f90):
  an example that transfers all members of an LHAPDF set into a structure
  with an array of tables, allowing for faster evaluation of all members in one go

- [with-lhapdf/hoppet_lhapdf.f90](with-lhapdf/hoppet_lhapdf.f90):
  The module that is used to help with the above tasks. This module is not included
  in the HOPPET library. Instead, currently, a user should copy it to their
  own project
