# NEWS for HOPPET

# Release 2.0.1, 18 October 2025

## Bug fixes

* resolved issues with Fortran-C interfaces for routines that take
  logical arguments, affecting structure function and QED
  components. Fixed a couple of structure function interface functions
  which had the wrong signature.
  
* other small fixes for issues picked up with link-time optimization
  (thanks to Alexander Puck Neuwirth for reporting both of the
  above). As a consequence the code can be compiled with LTO by adding
  the following flags.

```
cmake -DCMAKE_C_FLAGS="-flto" -DCMAKE_CXX_FLAGS="-flto" -DCMAKE_Fortran_FLAGS="-flto"
```


# Release 2.0.0, 10 October 2025

## Major new features:

* inclusion of approximate N3LO splitting functions, including code from 
  FHMPRUVV group (arXiv:2410.08089 and refs therein)

* inclusion of mass threshold matrices for VFNS N3LO evolution, including
  code from ABBFMSS group (arXiv:2510.02175 and refs therein)

* inclusion of structure function calculations, including code from (D)MVV
  group, (arXiv:1606.08907 and refs therein)

* QED evolution

* CMake build system (thanks to Andrii Verbitskyi for contributing the
  initial version). The v1 build system will be kept for the 2.0.x
  series but will not support all features.

* a Python interface (similar to the streamlined interface). The
  interface is built with SWIG (https://github.com/swig) and can be
  enabled with `-DHOPPET_BUILD_PYINTERFACE=ON`. Alternatively just `pip
  install hoppet`.
  
We are grateful to the authors of the N3LO coefficient functions,
splitting functions and mass threshold terms for providing their code. 

The updates are documented in https://arxiv.org/abs/2510.09310.

## Backwards incompatible changes

* hoppet_v1 now renamed hoppet (in library name, module name, C++
  include name, etc.)
* some QED defaults (lepton masses, effective light quark masses)  
  have been updated to 2024/25 values

## Authorship update: 

* Frederic Dreyer, Alexander Karlberg, Paolo Nason and Giulia Zanderighi
  are now authors

## Other changes:

* A compiler with support for Fortran 2008 is now required, and HOPPET relies
  on extensions for long-lines (part of the Fortran 2023 standard). Compilers
  also need to be able to run the preprocessor. All features are widely supported.

* access to a single flavour x,Q point of a pdf_table, with
  `EvalPdfTable_yQf` (and `_xQf` variant) and speed-up of
  `EvalPdfTable_yQ` by about a factor of three (300ns to 100 ns on an
  M2Pro). In the streamlined/C++ interface, use hoppetEvalIFlv(x,Q,iflv),
  with the predefined iflv_g, iflv_d, etc. symbols to get the right index. 
  
* Also addition of a variant of `EvalPdfTable_yQ` that can take
  an array of tables, e.g. for speed in evaluating the same x,Q point
  for error sets.

* Interpolation order for evaluation of PDF tables can now be overridden
  with the `PdfTableOverrideInterpOrders(...)` (F90 interface) and
  `hoppetSetYLnlnQInterpOrders(...)` calls (streamlined interface).
  Setting quadratic interpolation in y and Q, this reduces PDF
  evaluation time to about 60ns, and accuracy remains good with suitably
  fine grids (`dy=0.05`).

* functionality to print a PDF from a HOPPET table in LHAPDF6 format
  (`WriteLHAPDFFromPdfTable` and `hoppetWriteLHAPDFGrid`)

* ability to declare points where the integration of convolution
  functions should be split up (split optional arg to `InitGridConv`, etc.)

* the streamlined interface now has `hoppetDeleteAll()` function for cleaning up
  all memory it allocated.

* Implementation of continuous integration and a wider range of checks, including 
  a new `unit-tests/` directory

* reorganisation of the examples into `example_f90/`, `example_f77/`,
  `example_cpp/` and `example_py/`  

* various other small additions and bug fixes (see ChangeLog for details)


# Release 1.2.0, 30 October 2020

* hoppetEvalSplit now supports composite iloop value, e.g. iloop=12
  to get `PLO*PNLO*pdf`
* 4-loop running coupling is now available
* Frederic Dreyer and Alexander Karlberg joined as contributors
* small fixes to build system and examples, memory management


# Release 1.1.5  August 2012

* a PartonLuminosity(...) function has been added (see manual)
* modules are now installed by default with make install
* default module installation location is now PREFIX/include/hoppet/
* hoppet-config now accepts an --fflags option
* fixed NaN problems in tabulation at Q values below lower edge of table
* small fixes to build and check system, error messages


# Release 1.1.4  November 2011

* mass thresholds now available at MSbar masses as well as pole masses
* fixes for compatibility with gfortran 4.5 upwards


# Release 1.1.3  May 2010

* small bug fixes to build system
* addition of hoppet-config script


# Release 1.1.2  September 2009


* new build system (configure lookalike).

* added 1-loop time-like (fragmentation function) evolution following
  discussions with Wei-Tian Deng. Only partially tested.


# Release 1.1.1  August 2008


* added TruncatedMoment routine (e.g. for verifying sum-rules)

* added ability to extract a whole pdf(x,flav) at a given Q from a
  tabulation, with the EvalPdfTable_Q routine

* clarifications and other small changes in the documentation


# Release 1.1.0  April 2008


* Name changed to hoppet
* Many internal name changes to provide more consistent user-interface
* Some changes in setting of default grid and evolution parameters
* Documentation written


# Release 1.0 (named pdf-conv) (6 April 2006)


Imported from the PDFevln sources, with modifications of the directory
structure and addition of some f77 hooks.
