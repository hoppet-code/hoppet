HOPPET: Higher Order Perturbative Parton Evolution Toolkit
==========================================================

HOPPET is a Fortran 95 package for carrying out DGLAP evolution and
other common manipulations of parton distribution functions (PDFs).

Within HOPPET, PDFs are represented on a grid in x-space so as to
avoid limitations on the functional form of input
distributions. Decent speed and accuracy are obtained through the
representation of splitting functions in terms of their convolution
with a set of piecewise polynomial basis functions, and Runge-Kutta
techniques are used for the evolution in Q.

Unpolarised evolution is provided in the MSbar scheme to NNLO in the
coupling, including heavy-quark thresholds, and polarised evolution to
NLO. The code is structured so as to provide simple access to the
objects representing splitting function and PDFs, making it possible
for a user to extend the facilities already provided.

The latest version can always be obtained from

    git clone https://github.com/gavinsalam/hoppet

Details of changes are to be found in the file ChangeLog, while
summaries of changes between releases are in ReleaseNotes.


----------------------------------------------------------------------
F95 compilers
-------------

You will need a Fortran 95 compiler to compile this package. 

The gfortran compiler is widespread. Versions upwards of 4.3 are
OK. Older versions (which might be present on legacy systems) generate
incorrect code for HOPPET and should be avoided. The code is
acceptably fast, though in tests some years ago was not quite as fast
as with commercial compilers.

The code has also been tested with the intel (ifort, versions 8.1.033
upwards) and lahey (lf95) compilers.

Compilation
-----------
For details see the INSTALL file. To get moving quickly, just specify
an installation prefix and a fortran compiler (FC), and then do

    ./configure --prefix="..."  FC="..."
    make 
    make check
    make install      [if you're interested]
    make install-mod  [if you want the f90 module files installed too]

This is not autotools based: if you're used to more advanced usage of
autotools scripts, you'll be disappointed here...

Example programs
----------------

The main example program is examples_f90/tabulation_example.f90

An equivalent program with the startup, pdf initialisation and
evolution spread across different subroutines is given as
tabulation_example_2.

An equivalent program based on the streamlined interface is given as
tabulation_example_streamlined.

Some users may prefer a pure f77 interface. Corresponding examples are
to be found in the [examples_f77/](examples_f77) directory. Look inside
the suppplied Makefile and if need be edit it manually.

     cd ../example_f77
     # <edit the Makefile directly>
     # compile
     make tabulation_example
     # run the program. Should give output identical to that from
     # example_f90/tabulation_example
     ./tabulation_example

In the same directory there is a C++ example

     make cpp_tabulation_example
     ./cpp_tabulation_example

which again does the same things (though in this case it uses the
simpler of the two streamlined initialization calls).

Other programs provided in the example_f77/ directory illustrate the
use of the streamlined interface in conjunction with LHAPDF
(compare_lhapdf_hoppet.f), and show how to use the feature of getting
convolutions with splitting functions (convolution_example.f).

----------------------------------------------------------------------
Documentation
-------------

Detailed documentation is available as doc/HOPPET-v1-doc.tex .


----------------------------------------------------------------------
Benchmarking
------------

The benchmarking/ directory contains the programs used for the full
benchmarking, accuracy and precision testing. 


----------------------------------------------------------------------
Compiler warnings
-----------------

When the hoppet library is being built, on some systems the following
warnings may be reported 

    /usr/bin/ranlib: file: libhoppet_v1.a(hoppet_v1.o) has no symbols
    /usr/bin/ranlib: file: libhoppet_v1.a(types.o) has no symbols
    
    ranlib: file: libhoppet_v1.a(hoppet_v1.o) has no symbols
    ranlib: file: libhoppet_v1.a(types.o) has no symbols

These warnings do not indicate an actual problem. They are simply a
consequence of the fact that some of the source (.f90) files contain
only information about interfaces and constants, so that the resulting
object (.o) files are essentially empty.
