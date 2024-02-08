[![Build
Status](https://img.shields.io/github/actions/workflow/status/hoppet-code/hoppet/CI.yml?label=build&logo=github&style=flat-square)](https://github.com/hoppet-code/hoppet/actions/workflows/CI.yml)

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

    git clone https://github.com/hoppet-code/hoppet

Details of changes are to be found in the file ChangeLog, while
summaries of changes between releases are in ReleaseNotes.


----------------------------------------------------------------------
Note on F95 compilers
-------------

You will need a Fortran 95 compiler to compile this package. 

The gfortran compiler is widespread. Versions upwards of 4.3 are
OK. Older versions (which might be present on legacy systems) generate
incorrect code for HOPPET and should be avoided. The code is
acceptably fast, though in tests some years ago was not quite as fast
as with commercial compilers.

The code has also been tested with the intel (ifort, versions 8.1.033
upwards) and lahey (lf95) compilers.

Build with `configure` script
-----------------------------
For details see the INSTALL file. To get moving quickly, just specify
an installation prefix and a fortran compiler (FC), and then do

    ./configure --prefix="..."  FC="..."
    make 
    make check
    make install      [if you're interested]
    make install-mod  [if you want the f90 module files installed too]

This is not autotools based: if you're used to more advanced usage of
autotools scripts, you'll be disappointed here...

Build with `CMake`
------------------

   mkdir build
   cmake -S . -B build <extra flags>
   cmake --build  build -j 
   cmake --install build
   ctest --test-dir build

The extra flags might be:
- generic `CMake` flags, e.g. `-DCMAKE_INSTALL_PREFIX=/my/home/dir`, `-DCMAKE_Fortran_COMPILER=ifort`, `-DCMAKE_Fortran_FLAGS="-O2 -g"`, etc.
- flags pointing to the dependencies, `-DLHAPDF_DIR=/where/the/LHAPDF/is`,
- `HOPPET_USE_EXACT_COEF`    Use exact coefficient functions.
- `HOPPET_BUILD_EXAMPLES`    Build examples.
- `HOPPET_ENABLE_TESTING`    Enable testing. Requires building the examples.
- `HOPPET_BUILD_BENCHMARK`   Build benchmark.
- `HOPPET_ENABLE_FPES`       Enable trapping of the floating point exceptions. Recommended for usage with debug builds, 
   i.e. `-DCMAKE_BUILDTYPE=DEBUG`

Usage in `CMake`-based projects
-------------------------------
If hoppet is build with `CMake` it is possible to import its targets into other projects.
Namely, a minimal working project 

    cmake_minimum_required(VERSION 3.12)
    project(mytest LANGUAGES Fortran)
    find_package(hoppet)
    add_executable(tabulation_example tabulation_example.f90)
    target_link_libraries(tabulation_example PUBLIC hoppet::hoppet_static)

could be configured with 


    cmake -S . -B BUILD -Dhoppet_DIR=hoppet/installation/share/hoppet/cmake


Example programs
----------------

The main example program is
[example_f90/tabulation_example.f90](example_f90/tabulation_example.f90)

An equivalent program with the startup, pdf initialisation and
evolution spread across different subroutines is given as
tabulation_example_2.

An equivalent program based on the streamlined interface is given as
[tabulation_example_streamlined](example_f90/tabulation_example_streamlined.f90).

Some users may prefer a pure f77 interface. Corresponding examples are
to be found in the [example_f77/](example_f77) directory. 

For the `CMake` based build the examples are compiled  if 
`-DHOPPET_BUILD_EXAMPLES=ON` flag is set. To build the examples with 
`configure` look inside the suppplied Makefile and if need be edit it manually.

     cd ../example_f77
     # <edit the Makefile directly>
     # compile
     make tabulation_example
     # run the program. Should give output identical to that from
     # example_f90/tabulation_example
     ./tabulation_example

In the same directory there is a C++ [example](example_f77/cpp_tabulation_example.cc)

     make cpp_tabulation_example
     ./cpp_tabulation_example

which again does the same things (though in this case it uses the
simpler of the two streamlined initialization calls).

Other programs provided in the [example_f77/](example_f77) directory
illustrate the use of the streamlined interface in conjunction with
LHAPDF ([compare_lhapdf_hoppet.f](example_f77/compare_lhapdf_hoppet.f)),
and show how to use the feature of getting convolutions with splitting
functions (convolution_example.f).

There is also an example program illustrating the QED evolution
feature of hoppet in
[tabulation_example_qed_streamlined.f90](example_f90/tabulation_example_qed_streamlined.f90)
and one that uses the structure functions
[structure_functions_example.f90.](example_f90/structure_functions_example.f90).

----------------------------------------------------------------------
Documentation
-------------

Detailed documentation is available as
[doc/HOPPET-v1-doc.tex](doc/HOPPET-v1-doc.tex). If you use HOPPET in a
scientific publication, please refer to
http://arxiv.org/abs/0804.3755. 


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

----------------------------------------------------------------------
Branches
--------

The [master](https://github.com/hoppet-code/hoppet/tree/master) branch
provides the official HOPPET release. There are currently no other
branches under active development.

