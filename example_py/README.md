# Examples for the Python interface to HOPPET

- [tabulation_example.py](tabulation_example.py): the standard benchmark
  evolution example

- [tabulation_example_qed.py](tabulation_example_qed.py): illustrates
  evolution including QED effects

- [structure_function_example.py](structure_function_example.py):
  illustrates evaluation of structure functions

- [lhapdf_to_hoppet.py](lhapdf_to_hoppet.py): shows how to read in an
  LHAPDF grid and use it from HOPPET. This gives PDF evaluation that is
  about 5 times faster that LHAPDF's Python interface in version 6.5.5.
  Requires the LHAPDF Python interface.
