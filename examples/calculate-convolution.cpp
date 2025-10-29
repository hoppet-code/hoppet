/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <tuple>
#include <ome/ome.h>

// Example for a test function we want to convolve the OMEs with
double test_function(double x)
{
  return(std::pow(1-x,5));
}

int main()
{
  // strong coupling, as = \alpha_s/(4 \pi)
  double as = 0.118/(4.*M_PI);
  // squared quark mass (on-shell scheme)
  double m2 = std::pow(1.59,2);
  // squared renormalisation scale
  double mu2 = std::pow(100.,2);
  // number of massless quark flavours
  double NF = 3;
  // precompute the mass logarithm
  double LM = std::log(m2/mu2);

  // lower and upper boundary for the scan in x
  double x_min = 1.e-5;
  double x_max = 1.;
  // number of steps in x
  int n_steps = 500;
  // precompute the logarithmic step in x to take in each iteration
  double x_log_step = std::log(x_max/x_min)/static_cast<double>(n_steps);

  // Handle for the numerical integrator we want to use
  ome::integration_engine_gsl engine;

  // Requested absolute and relative errors
  // This is passed on to the integrator. GSL interpretes this as when either
  // one of the bounds is reached, the integration terminates.
  double eps_abs = 0.;
  double eps_rel = 1.e-10;

  // Prepare Mellin convolution:
  // - select the OME (we directly pass the rpd_distribution object, which is
  //   a container for the regular, plus and delta part, so that we don't have
  //   to handle them individually)
  // - select the test function we want to convolve with
  // - specify the numerical integration engine
  // - pass the numerical values for the physics parameters
  ome::mellin_convolution<double> conv
    = ome::make_mellin_convolution(ome::AqqQNSEven,
        std::function(test_function), engine, as, LM, NF);

  std::cout << std::left
            << std::setw(25) << "# x"
            << std::setw(25) << "convolution"
            << std::setw(25) << "abs. integration error"
            << std::endl;
  for(int i = 0; i < n_steps; ++i)
  {
    // compute the value of x at which to evaluate the convolution
    double x = x_min * std::exp(static_cast<double>(i) * x_log_step);

    // Perform Mellin convolution numerically
    auto res = conv.integrate(x, eps_abs, eps_rel);
    if(std::get<0>(res) != ome::integration_status::success)
      std::cerr << "Warning: numerical integration did not succeed." << std::endl;

    std::cout << std::setprecision(15) << std::left
              << std::setw(25) << x
              << std::setw(25) << std::get<1>(res)
              << std::setw(25) << std::get<2>(res)
              << std::endl;
  }

  return(0);
}
