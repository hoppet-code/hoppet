/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <gtest/gtest.h>
#include <cmath>
#include <ome/functions.h>
#include <ome/laurent_polynomial.h>

using namespace ome;

template<typename Tnum>
class eval_id
{
  public:
    Tnum operator()(const Tnum x) const
    {
      return(x);
    };
};

template<typename Tnum>
class eval_sum
{
  public:
    Tnum operator()(const Tnum x, const Tnum y) const
    {
      return(x+y);
    };
};

// Tests for func_shift

TEST(FunctionsFuncShiftTest, DefaultConstructIdentity)
{
  func_shift<double, eval_id<double>> fshift;

  EXPECT_EQ(fshift(42.),42.);
}

TEST(FunctionsFuncShiftTest, EvaluateUnaryTarget)
{
  func_shift<double, eval_id<double>> fshift(1., eval_id<double>());

  EXPECT_EQ(fshift(1.), 2.);
}

TEST(FunctionsFuncShiftTest, EvaluateBinaryTarget)
{
  func_shift<double, eval_sum<double>, double> fshift(1., eval_sum<double>());

  EXPECT_EQ(fshift(1., 10.), 12.);
}


// Tests for func_apply

TEST(FunctionsFuncApplyTest, DefaultConstructIdentity)
{
  func_apply<double, eval_id<double>> fapply;

  EXPECT_EQ(fapply(42.), 42.);
}

TEST(FunctionsFuncApplyTest, EvaluateUnaryTarget)
{
  func_apply<double, eval_id<double>> fapply(f_omx<double>, eval_id<double>());

  EXPECT_EQ(fapply(2.), -1.);
}

TEST(FunctionsFuncApplyTest, EvaluateBinaryTarget)
{
  func_apply<double, eval_sum<double>, double> fapply(f_omx<double>, eval_sum<double>());

  EXPECT_EQ(fapply(0.25, 10.), 10.75);
}


// Tests for func_copy_and_log

TEST(FunctionsFuncCopyAndLogTest, DefaultConstruct)
{
  func_copy_and_log<double, eval_sum<double>> fcal;

  EXPECT_DOUBLE_EQ(fcal(2.), 2. + std::log(2.));
}

TEST(FunctionsFuncCopyAndLogTest, Evaluate)
{
  func_copy_and_log<double, eval_sum<double>> fcal{eval_sum<double>()};

  EXPECT_DOUBLE_EQ(fcal(2.), 2. + std::log(2.));
}


// Test for func_plusfunc_omx

TEST(FunctionsFuncPlusFuncOmxTest, DefaultConstruct)
{
  func_plusfunc_omx<double, eval_id<double>> fplus;

  EXPECT_DOUBLE_EQ(fplus(0.25), std::log(0.75)/0.75);
}

TEST(FunctionsFuncPlusFuncOmxTest, EvaluateUnaryTarget)
{
  func_plusfunc_omx<double, eval_id<double>> fplus{eval_id<double>()};

  EXPECT_DOUBLE_EQ(fplus(0.25), std::log(0.75)/0.75);
}

TEST(FunctionsFuncPlusFuncOmxTest, EvaluateBinaryTarget)
{
  func_plusfunc_omx<double, eval_sum<double>, double> fplus{eval_sum<double>()};

  EXPECT_DOUBLE_EQ(fplus(0.25,10.), (std::log(0.75)+10.)/0.75);
}

TEST(FunctionsFuncPlusFuncOmxTest, EvalPlusInt)
{
  // testing eval_plus_int unfortunately only makes sense in conjunction with a
  // wrapped laurent_polynomial class
  func_plusfunc_omx<double, laurent_polynomial<double, double>>
    fplus(laurent_polynomial<double, double>({4.0,5.0,6.0}));

  double x = 0.4;
  double lomx = std::log(1.-x);
  EXPECT_DOUBLE_EQ(fplus.eval_plus_int(x),
    -4.*lomx - 2.5*std::pow(lomx,2) - 2.*pow(lomx,3));
}
