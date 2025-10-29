# libome                                                             {#mainpage}

This library provides numerical implementations of the QCD corrections up to
three-loop order to the massive operator matrix elements (OMEs) of the QCD
twist-2 operators in x-space.

## Overview

There are seven unpolarised OMEs that contribute up to three-loop order:
| OME                                    | name used in libome |
|----------------------------------------|---------------------|
| \f$A_{Qg}\f$                           | `AQg`               |
| \f$A_{gg,Q}\f$                         | `AggQ`              |
| \f$A_{Qq}^\text{PS}\f$                 | `AQqPS`             |
| \f$A_{qq,Q}^\text{PS}\f$               | `AqqQPS`            |
| \f$A_{gq,Q}\f$                         | `AgqQ`              |
| \f$A_{qg,Q}\f$                         | `AqgQ`              |
| \f$A_{qq,Q}^{\text{NS},+}\f$           | `AqqQNSEven`        |
| \f$A_{qq,Q}^{\text{NS},-}\f$           | `AqqQNSOdd`         |
Regarding the analytic continuation, the implementation of all of these OMEs
correspond to the analytic continuation from even moments, except for
\f$A_{qq,Q}^{\text{NS},-}\f$, which is analytically continued from odd moments.

Similarly, there are seven polarised OMEs:
| OME                                 | name used in libome    |
|-------------------------------------|------------------------|
| \f$\Delta A_{Qg}\f$                 | `polAQg`               |
| \f$\Delta A_{gg,Q}\f$               | `polAggQ`              |
| \f$\Delta A_{Qq}^\text{PS}\f$       | `polAQqPS`             |
| \f$\Delta A_{qq,Q}^\text{PS}\f$     | `polAqqQPS`            |
| \f$\Delta A_{gq,Q}\f$               | `polAgqQ`              |
| \f$\Delta A_{qg,Q}\f$               | `polAqgQ`              |
| \f$\Delta A_{qq,Q}^{\text{NS},+}\f$ | `polAqqQNSEven`        |
| \f$\Delta A_{qq,Q}^{\text{NS},-}\f$ | `polAqqQNSOdd`         |
They are analytically continued from odd moments, except for
\f$A_{qq,Q}^{\text{NS},+}\f$, which is analytically continued from even
moments. The polarised OMEs were calculated in the Larin scheme.

The OMEs depend on four parameters:
- the strong coupling \f$a_s = \frac{\alpha_s}{4 \pi} = \frac{g_s}{(4 \pi)^2}\f$
  in the \f$\overline{\text{MS}}\f$ scheme for \f$N_F+1\f$ active flavours
- the logarithm \f$L_M = \log\left(\frac{m^2}{\mu^2}\right)\f$, where \f$m\f$
  is the mass of the heavy quark in the on-shell scheme and \f$\mu\f$ is the
  renormalisation scale
- the number of massless quark flavours \f$N_F\f$
- the Bjorken variable \f$x\f$

The OMEs are distributions in \f$x\f$ with the following general structure
\f[
A(x) = A_\text{reg}(x) + [A_+(x)]_+ + \delta(1-x) A_\delta
\f]
that is, there are regular, plus and delta parts. The library follows the
convention to append the suffix `_reg` to the regular part
\f$A_\text{reg}(x)\f$, the suffix `_plus` to plus distributions \f$A_+(x)\f$ and
the suffix `_delta` to the coefficient \f$A_\delta\f$ of delta distributions.
All OMEs have regular parts, but only the OMEs \f$(\Delta) A_{gg,Q}\f$ and
\f$(\Delta) A_{qq,Q}^{\text{NS},\pm}\f$ have plus and delta parts. 


## Basic usage from C++

To use it from your own C++ code, include the default header
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
#include <ome/ome.h>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
which will include all relevant further headers. Everything in this library
lives in the namespace `ome`. You will then have access to constant global
objects that represent the massive OMEs and can be evaluated with numerical
values for the parameters as arguments.

### Evaluation

Evaluation of all orders up to the highest available one (i.e. three loop) then
is possible by simply calling, e.g.,
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
double as = 0.01, LM = 5., NF = 3., x = 0.2;
double res_reg   = ome::AqqQNSEven_reg(as, LM, NF, x);
double res_plus  = ome::AqqQNSEven_plus(as, LM, NF, x);
double res_delta = ome::AqqQNSEven_delta(as, LM, NF);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Note that the delta part only depends on `as`, `LM` and `NF`, but not on `x`
since it is only the coefficient of the \f$\delta(1-x)\f$ distribution.

These constant global objects are implemented as nestings of
[laurent_polynomial] classes (and further classes at deeper nesting levels).
If you are interested in these implementation details, see the
[list of public data types](#ome-data-types) for a more detailed description
of how classes are nested to achive this.

### Truncation

Since it is sometimes desiable to evaluate an OME only up to a certain order in
the strong coupling constant, it is possible to truncate the series in
\f$a_s\f$ at a specified order. There is a [truncate] method that implements
this idea
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
auto trunc_as2 = ome::AQg_reg.truncate(2);
double res_trunc_as2 = trunc_as2(as, LM, NF, x);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One technical comment: truncation is a relatively cheap operation. The object
returned by the `truncate` method is a [laurent_polynomial_view] which is a
lightweight object that only refers to a [laurent_polynomial], but does not
contain a copy of the coefficient data itself.

### Coefficient access

It is also possible to isolate a specific coefficient of the (Laurent)
polynomials. For example, if you need only the \f$O(a_s^3)\f$ coefficient
of \f$A_{Qq}^{\text{PS}}\f$ you can extract it via the [access operator]
like so
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
auto coeff_as3 = ome::AQqPS_reg[3];
double res_coeff_as3 = coeff_as3(LM, NF, x);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To evaluate the coefficient, you of course no longer have to pass `as` as an
argument.

It is also possible to then extract the coefficient of, e.g., \f$O(LM^0)\f$ by
using the access operator again
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
auto coeff_as3_LM0 = ome_coeff_as3[0];
double res_coeff_as3_LM0 = coeff_as3_LM0(NF, x);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To find out the minimal and maximal exponents of a (Laurent) polynomial use the
[min_power] and [max_power] methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
int min_pow = ome::AQqPS_reg.min_power();
int max_pow = ome::AQqPS_reg.max_power();
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Accessing a coefficient outside this range will return an object that will
always evaluate to zero.

### `rpd_distribution` containers

In order to simplify the handling of OMEs as objects with regular, plus and
delta parts, there is a specialised container called [rpd_distribution] which
holds these three parts (or rather those of them that the OME actually has).
The idea is that instead of having to pass around three objects separately (and
having to implement special cases for when one of these parts is absent) you
can instead just pass one [rpd_distribution] container that holds all three.

Since not every OME has all three parts, they are stored in `std::optional`
objects. Thus, the absence of one of the three parts can be modeled using
`std::nullopt`.

For each of the OMEs, there is also a corresponding global [rpd_distribution]
object:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
bool has_reg = ome::AggQ.has_regular();
auto reg = ome::AggQ.get_regular(); // equivalent to ome::AggQ_reg object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The [rpd_distribution] container also supports truncation and coefficient
access
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
auto rpd_trunc_as2 = ome::AggQ.truncate(2);
auto rpd_coeff_as2 = ome::AggQ.get_coefficient(2);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Calculating Mellin moments and convolutions

Finally, the library provides simple interfaces to numerically calculate Mellin
moments and convolutions. The numerical integrator is wrapped by an interface
called [integration_engine] for which the library provides a default
implementation, [integration_engine_gsl] based on the GNU Scientific Libarary
and its CQUAD numerical integration routine.

The class [mellin_moment] implements the computation of Mellin moments. For OMEs
it is probably most straightforward to construct the object directly from a
[rpd_distribution] using the factory function [make_mellin_moment]. This
function takes the [rpd_distribution], an [integration_engine] and a list of
values for the numerical parameters that are not being integrated (i.e., `as`,
`LM` and `NF` for OMEs) as arguments. The moment can then be computed using the
[mellin_moment::integrate] method, providing the moment `n` and the absolute and
relative precision target for the numerical integration as arguments.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
ome::integration_engine_gsl engine;
ome::mellin_moment<double> mom
  = make_mellin_moment(ome::AgqQ, engine, as, LM, NF);

double eps_abs = 1.e-15, eps_rel = 1.e-10;
int n = 4;
auto res = mom.integrate(n, eps_abs, eps_rel);
if(std::get<0>(res) != ome::integration_status::sucess)
  std::cerr << "Warning: numerical integration did not succeed." << std::endl;

std::cout << std::get<1>(res) << " +- " << std::get<2>(res) << std::endl;
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The [mellin_moment::integrate] method returns an `std::tuple` of an
[integration_status], a double representing the value of the moment computed
and a double representing and error estimate from the numerical integration.

Similarly, Mellin convolutions can be computed using the [mellin_convolution]
class. It can also be constructed from an [rpd_distribution] using the
[make_mellin_convolution] factory function.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
ome::mellin_convolution<double> conv
  = make_mellin_convolution(ome::AgqQ, std::function(test_function)
                            engine, as, LM, NF);
double x_conv = 0.05;
auto res = conv.integrate(x_conv, eps_abs, eps_rel);
if(std::get<0>(res) != ome::integration_status::sucess)
  std::cerr << "Warning: numerical integration did not succeed." << std::endl;

std::cout << std::get<1>(res) << " +- " << std::get<2>(res) << std::endl;
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here, we pass the function we convolve with as an `std::function`.


## Basic usage from C

While the library is written in C++, it also offers a simplified interface that
can also be [used from C](#c-interface). While it does not
offer the [rpd_distribution] container or the interfaces to numerically
calculate Mellin moments and convolutions, it does allow to evaluate the OMEs,
truncated OMEs or OME coefficients.

In order to use the library in your C code, you still just include the default
header
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
#include <ome/ome.h>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
which only exposes the C interface if it is processed by a C compiler.

For every OME listed above, there exists a set of functions with C linkage
that return the evaluation of that OME. The signature for the functions
to evaluate the full OME follows the schema
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
double ome_<name>_reg(double as, double LM, double NF, double x);
double ome_<name>_plus(double as, double LM, double NF, double x);
double ome_<name>_delta(double as, double LM, double NF);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
where `<rpd>` is either `reg`, `plus` or `delta`. If the OME only has a regular
part, only the `reg` function is provided.

In addition, there is a set of functions for evaluating the OMEs truncated to
a particular order in \f$a_s\f$:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
double ome_<name>_reg_trunc_as(int order_as, double LM, double NF, double x);
double ome_<name>_plus_trunc_as(int order_as, double LM, double NF, double x);
double ome_<name>_delta_trunc_as(int order_as, double LM, double NF);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also evaluate only a specific coefficient of \f$a_s\f$, \f$L_M\f$ or
\f$N_F\f$:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
double ome_<name>_reg_coeff_as(int order_as, double LM, double NF, double x);
double ome_<name>_plus_coeff_as(int order_as, double LM, double NF, double x);
double ome_<name>_delta_coeff_as(int order_as, double LM, double NF);

double ome_<name>_reg_coeff_as_LM(int order_as, int order_LM, double NF, double x);
double ome_<name>_plus_coeff_plus_LM(int order_as, int order_LM, double NF, double x);
double ome_<name>_delta_coeff_delta_LM(int order_as, int order_LM, double NF);

double ome_<name>_reg_coeff_as_LM_NF(int order_as, int order_LM, int order_NF, double x);
double ome_<name>_plus_coeff_as_LM_NF(int order_as, int order_LM, int order_NF, double x);
double ome_<name>_delta_coeff_as_LM_NF(int order_as, int order_LM, int order_NF, double x);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, to determine the exponent ranges of non-vanishing coefficients of the
polynomials in \f$a_s\f$, \f$L_M\f$ or \f$N_F\f$, you can use the following
functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
/* minimum and maximum power in as for the whole OME */
int ome_<name>_reg_min_power();
int ome_<name>_reg_max_power();
int ome_<name>_plus_min_power();
int ome_<name>_plus_max_power();
int ome_<name>_delta_min_power();
int ome_<name>_delta_max_power();

/* minimum and maximum power in LM for the specified coefficient of as */
int ome_<name>_reg_coeff_as_min_power(int order_as);
int ome_<name>_reg_coeff_as_max_power(int order_as);
int ome_<name>_plus_coeff_as_min_power(int order_as);
int ome_<name>_plus_coeff_as_max_power(int order_as);
int ome_<name>_delta_coeff_as_min_power(int order_as);
int ome_<name>_delta_coeff_as_max_power(int order_as);

/* minimum and maximum power in NF for the specified coefficient of as and LM */
int ome_<name>_reg_coeff_as_LM_min_power(int order_as, int order_LM);
int ome_<name>_reg_coeff_as_LM_max_power(int order_as, int order_LM);
int ome_<name>_plus_coeff_as_LM_min_power(int order_as, int order_LM);
int ome_<name>_plus_coeff_as_LM_max_power(int order_as, int order_LM);
int ome_<name>_delta_coeff_as_LM_min_power(int order_as, int order_LM);
int ome_<name>_delta_coeff_as_LM_max_power(int order_as, int order_LM);
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


[laurent_polynomial]: @ref ome::laurent_polynomial "laurent_polynomial"
[laurent_polynomial_view]: @ref ome::laurent_polynomial_view "laurent_polynomial_view"
[truncate]: @ref ome::laurent_polynomial_view::truncate "truncate"
[access operator]: @ref ome::laurent_polynomial_view::operator[] "access operator"
[min_power]: @ref ome::laurent_polynomial_view::min_power "min_power"
[max_power]: @ref ome::laurent_polynomial_view::max_power "max_power"
[rpd_distribution]: @ref ome::rpd_distribution
[integration_engine]: @ref ome::integration_engine "integration_engine"
[integration_engine_gsl]: @ref ome::integration_engine_gsl "integration_engine_gsl"
[mellin_moment]: @ref ome::mellin_moment "mellin_moment"
[make_mellin_moment]: @ref ome::make_mellin_moment "make_mellin_moment"
[mellin_moment::integrate]: @ref ome::mellin_moment::integrate "integrate"
[integration_status]: @ref ome::integration_status "integration_status"
[mellin_convolution]: @ref ome::mellin_convolution "mellin_convolution"
[make_mellin_convolution]: @ref ome::make_mellin_convolution "make_mellin_convolution"
