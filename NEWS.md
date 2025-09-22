# NEWS for HOPPET

# Release 2.0.0, XX September 2025

## Major new features:
* inclusion of QED evolution
* inclusion of structure function calculations, including code from (D)MVV
  group, (arXiv:1606.08907 and refs therein)
* inclusion of approximate N3LO splitting functions, including code from 
  FHMPRUVV group (arXiv:2410.08089 and refs therein)
* inclusion of mass threshold matrices for VFNS N3LO evolution, including
  code from ABBFMSS group (DESY 24-037 and refs therein)
* inclusion of a python interface (similar to the streamlined interface). The
  interface is built with SWIG (https://github.com/swig), and will not be
  compiled if SWIG is not present on the system

We are grateful to the authors of the N3LO coefficient functions,
splitting functions and mass threshold terms for providing their code. 

## Backwards incompatible changes
* hoppet_v1 now renamed hoppet (in library name, module name, C++
  include name, etc.)

## Authorship update: 
* Frederic Dreyer, Alexander Karlberg, Paolo Nason and Giulia Zanderighi
  are now authors

## Other changes:
* A compiler with support for Fortran 2008 is now required
* various small bug fixes (see ChangeLog for details)
* ability to declare points where the integration of convolution
  functions should be split up (split optional arg to InitGridConv, etc.)
* access to a single flavour x,Q point of a pdf_table, with EvalPdfTable_yQf
  (and _xQf variant) and speed-up of EvalPdfTable_yQ by about 20% for order=-6.
* streamlined interface now has hoppetDeleteAll() function for cleaning up
  all memory it allocated.
* Users have a choice between traditional makefiles and the CMake build
  system. CMake is now the preferred option. The two should not be
  used concurrently in the same directory.
* functionality to print a PDF from a hoppet table in the LHAPDF6 format is
  now provided, also in the streamlined interface. 


# Release 1.2.0, 30 October 2020

* hoppetEvalSplit now supports composite iloop value, e.g. iloop=12
  to get PLO*PNLO*pdf
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
