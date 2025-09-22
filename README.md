[![build](https://github.com/hoppet-code/hoppet/actions/workflows/main.yml/badge.svg)](https://github.com/hoppet-code/hoppet/actions/workflows/main.yml)

# HOPPET: Higher Order Perturbative Parton Evolution Toolkit

HOPPET is a Fortran package for carrying out DGLAP evolution and other
common manipulations of parton distribution functions (PDFs). It also
includes C++ and Python interfaces to many of the key features.

Within HOPPET, PDFs are represented on a grid in x-space so as to
avoid limitations on the functional form of input
distributions. Good speed and accuracy are obtained through the
representation of splitting functions in terms of their convolution
with a set of piecewise polynomial basis functions, and Runge-Kutta
techniques are used for the evolution in Q.

Key features include:

* unpolarized evolution in the MSbar scheme to N3LO including heavy-quark thresholds
* QED evolution
* massless structure function evaluation up to N3LO
* polarized evolution up to NLO
* simple access (in Fortran 90) to the objects representing splitting
  function and PDFs, making it possible for a user to extend the
  facilities already provided.

The latest version can always be obtained from

    git clone https://github.com/hoppet-code/hoppet

Summaries of main new features in each release are in
[NEWS.md](NEWS.md). More detailed lists of changes are in the
[ChangeLog](ChangeLog) file.



## Note on changes in v2.0.0


Starting from v2.0.0, backwards compatibility with the v1 series of
HOPPET has been broken in one significant way: `hoppet_v1`
has become `hoppet` in all places in the code. This in particular
applies to `libhoppet_v1` which is now `libhoppet`, `module hoppet_v1`
which is now `module hoppet` and the C++ namespace `hoppetv1` which
is now `hoppet`.

Besides this change most users should be able to start using the v2
series without any difficulties.

## Note on Fortran compilers

HOPPET requires support for Fortran 2008. It is developed and tested
with recent versions of the gfortran compiler (v12 and later) on Linux
and Mac, and also also gets tested with the Intel (ifx) compiler, as
part of the CI.

# Build with `CMake`


```bash
mkdir build
cmake -S . -B build <extra flags>
cmake --build  build -j 
ctest --test-dir build  -j
cmake --install build
```

Key `<extra flags>` include:
- generic CMake flags, e.g. 
  * `-DCMAKE_INSTALL_PREFIX=/my/home/dir`
  * `-DCMAKE_Fortran_COMPILER=ifx`,
  * `-DCMAKE_BUILD_TYPE=RelWithDebInfo` (default is `Release`, which gives highest speed; an additional `-g` is included by default to help with debugging), 
  * `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`, to see the exact compilation commands as a json file
- flags pointing to the dependencies, `-DLHAPDF_DIR=/where/the/LHAPDF/is`,
- `-DHOPPET_USE_EXACT_COEF=ON`:    Compile-in exact coefficient functions.



# Build with the configure script

This is the v1.x build system which will be supported to some extent
through the v2.0.x releases, but not necessarily beyond. It does not
support all features (e.g. the Python interface).

For details see the INSTALL file. To get moving quickly, just specify
an installation prefix and, optionally, a fortran compiler (FC), and then do

    ./configure --prefix="..."  FC="..."
    make 
    make check
    make install      [if you're interested]
    make install-mod  [if you want the f90 module files installed too]

This is not autotools based: if you're used to more advanced usage of
autotools scripts, you'll be disappointed here...


# Example programs

## Streamlined interface

The simplest way to use HOPPET is through its streamlined interface.
There are examples for basic QCD evolution in 
- modern Fortran:
[example_f90/tabulation_example_streamlined.f90](example_f90/tabulation_example_streamlined.f90)
- Fortran 77: [example_f77/tabulation_example.f](example_f77/tabulation_example.f) 
- C++: [example_cpp/tabulation_example.cc](example_cpp/tabulation_example.cc)
- Python: [example_py/tabulation_example.py](example_py/tabulation_example.py)

Other features illustrated in the examples include

- accessing convolutions with splitting functions ([Fortran](example_f77/convolution_example.f))
- QCD+QED evolution ([Fortran90](example_f90/tabulation_example_qed_streamlined.f90), [C++](example_cpp/tabulation_example_qed.cc), [Python](example_py/tabulation_example_qed.py))
- massless structure functions ([Fortran90](example_f90/structure_functions_example.f90), [C++](example_cpp/structure_functions_example.cc), [Python](example_py/structure_function_example.py))
- comparisons with LHAPDF ([Fortran](example_f77/compare_lhapdf_hoppet.f), [Fortran90](example_f90/with-lhapdf/fast_pdf_evaluation.f90), [C++](example_cpp/with-lhapdf/fast_pdf_evaluation.cc), [Python](example_py/tabulation_example_lhapdf.py))


## Full interface  
The Fortran examples also illustrate usage of the full lower-level interface,
allowing fine-grained access to the underlying objects, which provide a
powerful and flexible framework for convolution and evolution.
See for example
[example_f90/tabulation_example.f90](example_f90/tabulation_example.f90).
An equivalent program with the startup, pdf initialisation and
evolution spread across different subroutines is given as
[example_f90/tabulation_example_2.f90](example_f90/tabulation_example_2.f90).

## Building the examples

For the CMake-based builds, the Fortran and C++ examples are all
compiled by default (unless `-DHOPPET_BUILD_EXAMPLES=OFF`) and
executables are to be found in `build/example_*/` directories. Python
examples require HOPPET to be installed and for the HOPPET Python module
to be in the PYTHONPATH. Assuming `hoppet-config` is in your path, you
can view the correct path with `hoppet-config --pythonpath`. 


For the old build system (configure) only the Fortran 90 examples are
built. For other languages, go into the relevant directory, look inside
the supplied Makefile and if need be edit it manually.

```bash
cd example_f77
# <edit the Makefile directly>
# compile
make tabulation_example
# run the program. Should give output identical to that from
# example_f90/tabulation_example
./tabulation_example
```

# Compilation within your own projects

## Manual compilation and in Makefiles

The `hoppet-config` scripts provides flags to aid with compilation. E.g.

```sh
gfortran -O3 -g tabulation_example.f90 -o tabulation_example $(hoppet-config --fflags --libs)
```

Run `hoppet-config -h` to see the full list of flags.

## Usage in CMake-based projects

If hoppet was built with CMake, it is possible to import its targets into
other projects. Namely, a minimal working project 

```CMake
cmake_minimum_required(VERSION 3.12)
project(mytest LANGUAGES Fortran)
find_package(hoppet)
add_executable(tabulation_example tabulation_example.f90)
target_link_libraries(tabulation_example PUBLIC hoppet::hoppet_static)
```

could be configured with 

```sh
cmake -S . -B BUILD -Dhoppet_DIR=hoppet/installation/share/hoppet/cmake
```

# Documentation

Detailed documentation is available as
[doc/HOPPET-v1-doc.tex](doc/HOPPET-v1-doc.tex). If you use HOPPET in a
scientific publication, please refer to https://arxiv.org/abs/0804.3755
(original v1 release) and https://arxiv.org/abs/2509.nnnnn (v2 features).


# Other notes

## Benchmarking

The benchmarking/ directory contains the programs used for the full
benchmarking, accuracy and precision testing. 

## Compiler warnings

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


## Branches


The [master](https://github.com/hoppet-code/hoppet/tree/master) branch
provides the latest development version of HOPPET.

