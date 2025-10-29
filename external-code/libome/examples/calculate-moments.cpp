/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <ome/ome.h>

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

  // We want to calculate the even moments n=2,4,...,10
  int n_min = 2;
  int n_max = 10;
  int n_step = 2;

  // Handle for the numerical integrator we want to use
  ome::integration_engine_gsl engine;

  // Requested absolute and relative errors
  // This is passed on to the integrator. GSL interpretes this as when either
  // one of the bounds is reached, the integration terminates.
  double eps_abs = 0.;
  double eps_rel = 1.e-10;

  // Prepare Mellin moment calculation:
  // - select the OME (we directly pass the rpd_distribution object, which is
  //   a container for the regular, plus and delta part, so that we don't have
  //   to handle them individually)
  // - specify the numerical integration engine
  // - pass the numerical values for the physics parameters
  ome::mellin_moment<double> mom
    = ome::make_mellin_moment(ome::AqqQNSEven, engine, as, LM, NF);

  std::cout << std::left
            << std::setw(5) << "# n"
            << std::setw(25) << "moment"
            << std::setw(25) << "abs. integration error"
            << std::endl;
  for(int n = n_min; n <= n_max; n += n_step)
  {
    // Perform Mellin transformation numerically
    auto res = mom.integrate(n, eps_abs, eps_rel);
    if(std::get<0>(res) != ome::integration_status::success)
      std::cerr << "Warning: numerical integration did not succeed." << std::endl;

    std::cout << std::setprecision(15) << std::left
              << std::setw(5)  << n
              << std::setw(25) << std::get<1>(res)
              << std::setw(25) << std::get<2>(res)
              << std::endl;
  }

  return(0);
}
