/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
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

  // lower and upper boundary for the scan in x
  double x_min = 1.e-5;
  double x_max = 1.;
  // number of steps in x
  int n_steps = 500;
  // precompute the logarithmic step in x to take in each iteration
  double x_log_step = std::log(x_max/x_min)/static_cast<double>(n_steps);

  // step through x from x_min to x_max in logarithmically distributed steps
  // and evaluate the omes at each point
  for(int i = 0; i < n_steps; ++i)
  {
    // compute the value of x at which to evaluate the OME
    double x = x_min * std::exp(static_cast<double>(i) * x_log_step);

    std::cout << std::setprecision(15) << std::left
              << std::setw(25) << x
              << std::setw(25) << ome::AqqQNSEven_reg(as,LM,NF,x)
              << std::setw(25) << ome::AqqQNSEven_plus(as,LM,NF,x)
              << std::endl;
  }

  return(0);
}
