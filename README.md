[![build](https://github.com/hoppet-code/hoppet/actions/workflows/main.yml/badge.svg)](https://github.com/hoppet-code/hoppet/actions/workflows/main.yml)
<a href="https://pypi.org/project/HOPPET/"><img alt="PyPI" src="https://img.shields.io/pypi/v/HOPPET"/></a>
[![Documentation](https://img.shields.io/github/actions/workflow/status/hoppet-code/hoppet/build-wheels.yml?label=Documentation)](https://hoppet-code.github.io/hoppet/)


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
* simple access (in Fortran) to the objects representing splitting
  function and PDFs, making it possible for a user to extend the
  facilities already provided.

The latest version can always be obtained from

    git clone https://github.com/hoppet-code/hoppet

If you want a specific version, check out the relevant tag, e.g. for a
version 2.x.y, do

    cd hoppet
    git switch -d hoppet-2.x.y

Summaries of main new features in each release are in
[NEWS.md](NEWS.md). More detailed lists of changes are in the
[ChangeLog](ChangeLog) file.


## Documentation

Detailed documentation is available at
https://hoppet-code.github.io/hoppet/ and can be built also from source,
see [docs/README.md](docs/README.md).


If you use HOPPET in a scientific publication, please refer to
[arXiv:0804.3755](https://arxiv.org/abs/0804.3755) (original v1 release) and
[arXiv:2510.09310](https://arxiv.org/abs/2510.09310) (v2 features). The latter is now the
primary reference.



## Note on changes in v2 series

Starting from v2.0.0, backwards compatibility with the v1 series of
HOPPET has been broken in one significant way: `hoppet_v1`
has become `hoppet` in all places in the code. This in particular
applies to `libhoppet_v1` which is now `libhoppet`, `module hoppet_v1`
which is now `module hoppet` and the C++ namespace `hoppetv1` which
is now `hoppet`.

As of v2.1.0 the old hand-coded `./configure` script is no longer
supported, one must instead use CMake to build HOPPET.

Besides these changes most users should be able to start using the v2
series without any difficulties.


## Note on Fortran compilers

HOPPET requires support for Fortran 2008. It is developed and tested
with recent versions of the gfortran compiler (v10 and later) on Linux
and Mac, and also also gets tested with the Intel (ifx) compiler, as
part of the CI.

# Build


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

# Example programs

## Streamlined interface

The simplest way to use HOPPET is through its streamlined interface.
There are examples for basic QCD evolution in 
- modern Fortran:
[examples/f90/tabulation_example_streamlined.f90](examples/f90/tabulation_example_streamlined.f90)
- Fortran 77: [examples/f77/tabulation_example.f](examples/f77/tabulation_example.f) 
- C++: [examples/cpp/tabulation_example.cc](examples/cpp/tabulation_example.cc)
- Python: [examples/python/tabulation_example.py](examples/python/tabulation_example.py)

Other features illustrated in the examples include

- accessing convolutions with splitting functions ([Fortran](examples/f77/convolution_example.f))
- QCD+QED evolution ([Fortran90](examples/f90/tabulation_example_qed_streamlined.f90), [C++](examples/cpp/tabulation_example_qed.cc), [Python](examples/python/tabulation_example_qed.py))
- massless structure functions ([Fortran90](examples/f90/structure_functions_example.f90), [C++](examples/cpp/structure_functions_example.cc), [Python](examples/python/structure_function_example.py))
- comparisons with LHAPDF ([Fortran](examples/f77/compare_lhapdf_hoppet.f), [Fortran90](examples/f90/with-lhapdf/lhapdf_to_hoppet.f90), [C++](examples/cpp/with-lhapdf/lhapdf_to_hoppet.cc), [Python](examples/python/lhapdf_to_hoppet.py))


## Full interface  
The Fortran examples also illustrate usage of the full lower-level interface,
allowing fine-grained access to the underlying objects, which provide a
powerful and flexible framework for convolution and evolution.
See for example
[examples/f90/tabulation_example.f90](examples/f90/tabulation_example.f90).
An equivalent program with the startup, pdf initialisation and
evolution spread across different subroutines is given as
[examples/f90/tabulation_example_2.f90](examples/f90/tabulation_example_2.f90).

## Building the examples

The Fortran and C++ examples are all
compiled by default (unless `-DHOPPET_BUILD_EXAMPLES=OFF`) and
executables are to be found in `build/example_*/` directories. 

Python examples require a build with `-DHOPPET_BUILD_PYINTERFACE=ON`,
HOPPET to be installed and for the HOPPET Python module to be in the
PYTHONPATH. Assuming `hoppet-config` is in your path, you can view the
correct path with `hoppet-config --pythonpath`.  Alternatively you can 
simply `pip install hoppet` and then run the Python examples.

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


# Other notes

## Benchmarking

The benchmarking/ directory contains the programs used for the full
benchmarking, accuracy and precision testing. 

## Branches

- The [dev](https://github.com/hoppet-code/hoppet/tree/dev) branch
provides the latest development version of HOPPET. 
- The [master](https://github.com/hoppet-code/hoppet/tree/branch) will
usually be a little behind the dev branch and is intended to be suitable
for production use.


