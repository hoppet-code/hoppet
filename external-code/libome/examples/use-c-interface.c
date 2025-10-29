/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <ome/ome.h>

int main()
{
  // strong coupling, as = \alpha_s/(4 \pi)
  double as = 0.118/(4.*M_PI);
  // squared quark mass (on-shell scheme)
  double m2 = pow(1.59,2);
  // squared renormalisation scale
  double mu2 = pow(100.,2);
  // number of massless quark flavours
  double NF = 3;
  // precompute the mass logarithm
  double LM = log(m2/mu2);
  // Bjorken x
  double x = 0.2;

  printf("=== AqqQNSEven ===\n");

  printf("reg: %f\n", ome_AqqQNSEven_reg(as, LM, NF, x));
  printf("plus: %f\n", ome_AqqQNSEven_plus(as, LM, NF, x));
  printf("delta: %f\n\n", ome_AqqQNSEven_delta(as, LM, NF));

  printf("reg, trunc(as^2): %f\n", ome_AqqQNSEven_reg_trunc_as(2, as, LM, NF, x));
  printf("plus, trunc(as^2): %f\n", ome_AqqQNSEven_plus_trunc_as(2, as, LM, NF, x));
  printf("delta, trunc(as^2): %f\n\n", ome_AqqQNSEven_delta_trunc_as(2, as, LM, NF));

  printf("reg, coeff([as,3]): %f\n", ome_AqqQNSEven_reg_coeff_as(3, LM, NF, x));
  printf("plus, coeff([as,3]): %f\n", ome_AqqQNSEven_plus_coeff_as(3, LM, NF, x));
  printf("delta, coeff([as,3]): %f\n\n", ome_AqqQNSEven_delta_coeff_as(3, LM, NF));

  printf("reg, coeff([as,3],[LM,0]): %f\n", ome_AqqQNSEven_reg_coeff_as_LM(3, 0, NF, x));
  printf("plus, coeff([as,3],[LM,0]): %f\n", ome_AqqQNSEven_plus_coeff_as_LM(3, 0, NF, x));
  printf("delta, coeff([as,3],[LM,0]): %f\n\n", ome_AqqQNSEven_delta_coeff_as_LM(3, 0, NF));

  printf("reg, coeff([as,3],[LM,0],[NF,1]): %f\n", ome_AqqQNSEven_reg_coeff_as_LM_NF(3, 0, 1, x));
  printf("plus, coeff([as,3],[LM,0],[NF,1]): %f\n", ome_AqqQNSEven_plus_coeff_as_LM_NF(3, 0, 1, x));
  printf("delta, coeff([as,3],[LM,0],[NF,1]): %f\n\n", ome_AqqQNSEven_delta_coeff_as_LM_NF(3, 0, 1));

  printf("reg, min_power -- max_power: as^%d -- as^%d\n",
         ome_AqqQNSEven_reg_min_power(),
         ome_AqqQNSEven_reg_max_power());
  printf("plus, min_power -- max_power: as^%d -- as^%d\n",
         ome_AqqQNSEven_plus_min_power(),
         ome_AqqQNSEven_plus_max_power());
  printf("delta, min_power -- max_power: as^%d -- as^%d\n",
         ome_AqqQNSEven_delta_min_power(),
         ome_AqqQNSEven_delta_max_power());

  printf("reg, coeff([as,3]), min_power -- max_power: LM^%d -- LM^%d\n",
         ome_AqqQNSEven_reg_coeff_as_min_power(3),
         ome_AqqQNSEven_reg_coeff_as_max_power(3));
  printf("plus, coeff([as,3]), min_power -- max_power: LM^%d -- LM^%d\n",
         ome_AqqQNSEven_plus_coeff_as_min_power(3),
         ome_AqqQNSEven_plus_coeff_as_max_power(3));
  printf("delta, coeff([as,3]), min_power -- max_power: LM^%d -- LM^%d\n",
         ome_AqqQNSEven_delta_coeff_as_min_power(3),
         ome_AqqQNSEven_delta_coeff_as_max_power(3));

  printf("reg, coeff([as,3],[LM,0]), min_power -- max_power: NF^%d -- NF^%d\n",
         ome_AqqQNSEven_reg_coeff_as_LM_min_power(3,0),
         ome_AqqQNSEven_reg_coeff_as_LM_max_power(3,0));
  printf("plus, coeff([as,3],[LM,0]), min_power -- max_power: NF^%d -- MF^%d\n",
         ome_AqqQNSEven_plus_coeff_as_LM_min_power(3,0),
         ome_AqqQNSEven_plus_coeff_as_LM_max_power(3,0));
  printf("delta, coeff([as,3],[LM,0]), min_power -- max_power: NF^%d -- NF^%d\n",
         ome_AqqQNSEven_delta_coeff_as_LM_min_power(3,0),
         ome_AqqQNSEven_delta_coeff_as_LM_max_power(3,0));

  return(0);
}
