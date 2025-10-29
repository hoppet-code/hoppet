/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <ome/laurent_polynomial.h>

using namespace ome;
using poly = laurent_polynomial<double, double>;

TEST(LaurentPolynomialTest, DefaultConstructsZero)
{
  // Default construction for single-variable polynomials
  poly zero;
  EXPECT_EQ(zero(1.),0.)
    << "Default constructed laurent_polynomial doesn't evaluate to zero";

  // Default construction for nested polynomials
  laurent_polynomial<double, poly, double> zero_nested;
  EXPECT_EQ(zero_nested(1.,1.),0.)
    << "Default constructed nested laurent_polynomial doesn't evaluate to zero";
}

TEST(LaurentPolynomialTest, MinPowerConstruction)
{
  // Default minimum power is x^0
  poly default_min_power({1.});
  EXPECT_EQ(default_min_power(2.),1.)
    << "Default min_power is not x^0";

  // Negative minimum power
  poly negative_min_power({1.},-1);
  EXPECT_EQ(negative_min_power(2.),0.5)
    << "Negative min_power doesn't work";

  // Positive minimum power
  poly positive_min_power({1.},1);
  EXPECT_EQ(positive_min_power(2.),2.)
    << "Positive min_power doesn't work";
}

TEST(LaurentPolynomialTest, Evaluation)
{
  // Evaluation of an empty polynomial
  poly empty({},-2);
  EXPECT_EQ(empty(10.),0.)
    << "Evaluation of empty laurent_polynomial doesn't match";

  // Evaluation of a single-variable monomial
  poly mono({42.},2);
  EXPECT_EQ(mono(10.),4200.)
    << "Evaluation of single-variable monomial doesn't match";

  // Evaluation for single-variable polynomial
  // x+2*x^2+3*x^3+4*x^4
  // evaluated for x=10
  poly eval({1.,2.,3.,4.},1);
  EXPECT_DOUBLE_EQ(eval(10.),43210.)
    << "Evaluation of single-variable laurent_polynomial does not match";

  // Evaluation for nested polynomial
  // (1+y+y^2)/x + (-1/2/y + 2 + 3*y) + (4+y-y^2)*x
  // Evaluated for x=10,y=2
  laurent_polynomial<double, poly, double> eval_nested({
    poly({1.,1.,1.}),
    poly({-0.5,2.,3.},-1),
    poly({4.,1.,-1.})
  },-1);
  EXPECT_DOUBLE_EQ(eval_nested(10.,2.),28.45)
    << "Evaluation of nested laurent_polynomial does not match";
}

TEST(LaurentPolynomialTest, EvaluationWithPrecomputedMonomials)
{
  poly eval({1.,2.,3.,4.});
  std::vector<double> monomials({1000.,100.,10.,1.});

  EXPECT_DOUBLE_EQ(eval.eval_subst(monomials),1234.)
    << "Evaluation of single-variable polynomail with precomputed monomials does not match";

  std::vector<double> monomials2({1000.,1.});
  laurent_polynomial<double, poly, double> eval_nested({
    poly({3.,2.,1.}),
    poly({6.,5.,4.})
  });
  EXPECT_DOUBLE_EQ(eval_nested.eval_subst(monomials2,10.),123456.)
    << "Evaluation of nested laurent_polynomial with precomputed monomials does not match";
}

TEST(LaurentPolynomialTest, MinMaxPower)
{
  poly p({1.,2.,3.,4.,5.,6.},-2);

  // Query min power
  EXPECT_EQ(p.min_power(),-2)
    << "Minimum power doesn't match";

  // Query max power
  EXPECT_EQ(p.max_power(),3)
    << "Maximum power doesn't match";
}

TEST(LaurentPolynomialTest, MinMaxPowerOnEmpty)
{
  poly empty({},-5);

  // Query min power of empty polynomial
  EXPECT_EQ(empty.min_power(),-5)
    << "Minimum power of empty polynomial doesn't match";

  // Query max power of empty polynomial
  EXPECT_EQ(empty.max_power(),-6)
    << "Maximum power of empty polynomial doesn't match";
}

TEST(LaurentPolynomialTest, ElementAccess)
{
  poly p({1.,2.,3.,4.,5.,6.},-2);

  // Access too small power
  {
    SCOPED_TRACE("Coefficients before min_power do not match");
    EXPECT_EQ(p[-4], 0.);
    EXPECT_EQ(p[-3], 0.);
  }

  // Access first element
  EXPECT_EQ(p[-2], 1.)
    << "First element doesn't match";

  // Access middle element
  {
    SCOPED_TRACE("Middle elements do not match");
    EXPECT_EQ(p[-1], 2.);
    EXPECT_EQ(p[0], 3.);
    EXPECT_EQ(p[1], 4.);
    EXPECT_EQ(p[2], 5.);
  }

  // Access last element
  EXPECT_EQ(p[3], 6.)
    << "Last element doesn't match";

  // Access too large power
  {
    SCOPED_TRACE("Coefficients beyond max_power do not match");
    EXPECT_EQ(p[4], 0.);
    EXPECT_EQ(p[5], 0.);
  }
}

TEST(LauretPolynomialTest, ElementAccessOnNested)
{
  laurent_polynomial<double, poly, double> p_nested({
    poly({1.,1.,1.}),
    poly({-0.5,2.,3.},-1),
    poly({4.,1.,1.})
  },-1);

  SCOPED_TRACE("Element access on nested polynomial");

  // Below min_power
  EXPECT_EQ(p_nested[-2](10.),0.)
    << "Access below min_power doesn't match";

  // Access first element
  EXPECT_EQ(p_nested[-1](10.),111.)
    << "First element doesn't match";

  // Access last element
  EXPECT_EQ(p_nested[1](10.),114.)
    << "Last element doesn't match";

  // Above max_power
  EXPECT_EQ(p_nested[2](10.),0.)
    << "Access beyond max_power doesn't match";
}

TEST(LaurentPolynomialTest, Truncation)
{
  poly p({1.,2.,3.,4.},-2);

  auto p_trunc_m3 = p.truncate(-3);
  EXPECT_EQ(p_trunc_m3(10.),0.)
    << "Truncate below min_power doesn't match";

  auto p_trunc_m2 = p.truncate(-2);
  EXPECT_DOUBLE_EQ(p_trunc_m2(10.),0.01)
    << "Truncate at min_power doesn't match";

  auto p_trunc_m1 = p.truncate(-1);
  EXPECT_DOUBLE_EQ(p_trunc_m1(10.),0.21)
    << "Truncate at min_power+1 doesn't match";

  auto p_trunc_0 = p.truncate(0);
  EXPECT_DOUBLE_EQ(p_trunc_0(10.),3.21)
    << "Truncate at min_power+2 doesn't match";

  auto p_trunc_1 = p.truncate(1);
  EXPECT_DOUBLE_EQ(p_trunc_1(10.),p(10.))
    << "Truncate at max_power doesn't match";

  auto p_trunc_2 = p.truncate(2);
  EXPECT_DOUBLE_EQ(p_trunc_2(10.),p(10.))
    << "Truncate at max_power+1 doesn't match";
}


class eval_plus_int_haver
{
  public:
    eval_plus_int_haver(double value)
      : value_(value) {};

    double eval_plus_int(double x) const
    {
      return(value_/x);
    };

  private:
    double value_;
};

TEST(LaurentPolynomialTest, EvalPlusInt)
{
  using epih = eval_plus_int_haver;
  using poly_epih = laurent_polynomial<double, epih, double>;

  laurent_polynomial<double, poly_epih, double, double>
    p({poly_epih({epih(1.), epih(2.), epih(3.)}),
       poly_epih({epih(4.), epih(5.), epih(6.)})});

  constexpr double x = 10000., y = 10., z = 8.;
  constexpr double ref_val
    =  (1./z + 2./z * y + 3./z * y*y)
      +(4./z + 5./z * y + 6./z * y*y)*x;
  EXPECT_DOUBLE_EQ(p.eval_plus_int(x,y,z), ref_val);
}
