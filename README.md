# libome

Implementation of massive operator matrix elements (OMEs) of the QCD twist-2
operators in x-space.

This library provides numerical implementations of the QCD corrections up to
three-loop order to the massive OMEs that enter both heavy-flavour Wilson
coefficients for deep-inelastic scattering and the matching relations for
variable flavour number scheme (VFNS) parton distribution functions (PDFs).

The results implemented here were obtained in a series of publications. If you
use this library, please cite those references, as provided in the file
[CITATION](CITATION). For your convenience, the BibTeX entries for these
references are collected in [CITATION.bib](CITATION.bib). The code is subject
to the GNU General Public License version 3.0 or later, see [LICENSE](LICENSE).

The library is written C++ and in addition offers a C interface that can also
be called from Fortran. See
[examples/use-c-interface.c](examples/use-c-interface.c) an
[examples/use-c-interface-from-fortran.f90](examples/use-c-interface-from-fortran.f90)
for examples of how to use it. Moreover, the library offers a couple of
convenience features that allow, e.g., to compute Mellin moments and Mellin
convolutions through numerical integration.

The results are implemented as generalised series expansions in the Bjorken
variable x around different expansion points, so that for most points the
evaluation should almost reach double precision.

## Dependencies

The library is written in C++17 and thus needs a C++ compiler that supports
this standard. In addition, it depends on:
- GNU Scientific Library (GSL) (tested with version 2.8)
- CMake (at least version 3.28)

If you want to compile the Fortran example, you also need a Fortran compiler,
but this is not built by default.

To build the documentation, you need Doxygen installed.

Finally, if you want to build the unit and integration tests, the GoogleTest
framework is automatically downloaded and built. Also this is not enabled by
default.

## Installation

To compile the library, first create a build directory and go there:
```
mkdir build
cd build
```
Then, run CMake to generate the build scripts. By default GNU Makefiles are
generated, but other build systems, like Ninja are also supported. In the
simplest case, running
```
cmake ..
```
should be enough. Details about what can be configured are discussed below.
Then actually compile the library using `make`
```
make
```
and finally install it using (may require root privileges; see below for how to
change the install directory)
```
make install
```

The CMake setup offers a couple of configuration options. They can be changed
using `-D...` command line options. The most important ones are
- `-DCMAKE_INSTALL_PREFIX=/path/to/install/to`
  Change the directory to install to. By default the library is installed to
  `/usr/local`, but if you want to install it elsewhere, this option allows
  you to change the install destination.
- `-DBUILD_EXAMPLES=true`
  Enable the compilation of the examples in the `examples/` directory. By
  default, examples are not built. Even when built, the example programs are
  never installed along with the library, but they can be found in the
  `build/examples` directory and can be run from there.
- `-DBUILD_FORTRAN_EXAMPLE=true`
  Enable the compilation of the example for how to use the C interface from
  Fortran (`examples/use-c-interface-from-fortran.f90`).
- `-DBUILD_TESTS=true`
  Enable building the unit and integration test suite. CMake will automatically
  download and build the GoogleTest framework that the test suite uses.
  The tests are built in `build/tests` and can be run by executing the command
  `ctest` in that directory.
- `-DCMAKE_BUILD_TYPE=Debug`
  If you need debug symbols, this option enables them. It also turns off
  compile time optimisation, which will make the code noticeably slower.
- `-DBUILD_SHARED_LIBS=Off`
  By default, the library is compiled as a shared library. If you use this
  option, a static library is built instead.
- `-DENABLE_FPES=On`
  If you plan to rely on floating point exceptions (FPEs) in your code, enable
  this to ensure that the code does not raise spurious FPEs due to compiler
  optimisation steps.

## Usage

Documentation about how to use the library in your own programs is available
via Doxygen and can be generated using
```
cd docs
doxygen
```
and viewed by pointing a browser to `docs/html/index.html`.

Simple usage examples can be found in the `examples/` directory. (See also the
`-DBUILD_EXAMPLES=true` CMake option explained above.)

The library is built as a shared library. To use it from another program, you
have to link to it. Assume you have written a program `example.cpp` with the
following contents
```
#include <iostream>
#include <cmath>
#include <ome/ome.h>

int main()
{
  double as = 0.118/(4.*M_PI), LM = std::log(1.59/100.), NF = 3., x = 0.2;
  std::cout << ome::AqqQPS_reg(as, LM, NF, x) << std::endl;
  return(0);
}
```

If the library is installed to a standard path (e.g,.
`/usr` or `/usr/local`) that is searched by default, it should be enough to just
specify the library (and its dependencies) at link time:
```
g++ -lome -o example example.cpp
```

If you have installed the library to a non-standard path, you will have to
tell the compiler where to find the header files and the linker where to find
the shared library. Suppose you have installed the library to `/opt/ome`. In
that case, you would have to use
```
g++ -I/opt/ome/include -L/opt/ome/lib -lome -o example example.cpp
```
In addition, you will need to specify the path to the shared library at runtime.
On Linux, this can, for example, be done using the `LD_LIBRARY_PATH` environment
variable.
```
LD_LIBRARY_PATH="/opt/lib:$LD_LIBRARY_PATH" ./example
```

Finally, the library also provides pkg-config information. If the installation
directory of libome is on the pkg-config search path (either because it is
installed in a standard location or by adding it to the `PKG_CONFIG_PATH`
environment variable) you can use the `pkg-config` command to obtain the
necessary compiler and linker flags. In the running example from above you can
use
```
export PKG_CONFIG_PATH="/opt/lib/pkgconfig:$PKG_CONFIG_PATH"
g++ `pkg-config --cflags --libs libome` -o example example.cpp
```
