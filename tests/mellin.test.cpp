/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <cmath>
#include <functional>
#include <gsl/gsl_sf_dilog.h>
#include <gtest/gtest.h>
#include <ome/rpd_distribution.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>

#define EXPECT_INT_RESULT(res, ref_val, eps_rel) \
  EXPECT_EQ(std::get<0>(res), integration_status::success); \
  EXPECT_NEAR(std::get<1>(res), ref_val, std::abs(ref_val) * eps_rel); \
  EXPECT_LE(std::get<2>(res), std::abs(ref_val) * eps_rel);


using namespace ome;

// Implementation of alternating harmonic sum S_{-1}(n) for reference values
double Sm1(int n)
{
  double res = 0.;

  for(int i = 1; i <= n; ++i)
  {
    res += (i%2 == 0 ? 1. : -1.)/static_cast<double>(i);
  }

  return(res);
}

// Implementation of nested harmonic sum S_{1,1}(n) for reference values
double S11(int n)
{
  double res = 0.;

  for(int i = 1; i <= n; ++i)
  {
    for(int j = 1; j <= i; ++j)
    {
      res += 1./static_cast<double>(i*j);
    }
  }

  return(res);
}

// Example for the regular part
class mellin_example_reg
{
  public:
    double operator()(double x) const
    {
      return(1./(1+x));
    };

    double operator()(double a, double x) const
    {
      return(a*operator()(x));
    };

    static
    double mellin_mom(int n)
    {
      // Regular part transforms to:
      // (-1)^(n-1)*(ln2 + S[-1,n-1])
      int nm1 = n-1;

      return((nm1%2==0 ? 1. : -1)*(std::log(2.) + Sm1(nm1)));
    };

    static
    double mellin_mom(double a, int n)
    {
      return(a * mellin_mom(n));
    };
};

// Example for the plus part
class mellin_example_plus
{
  public:
    double operator()(double x) const
    {
      return(-std::log(1.-x)/(1.-x));
    };

    double operator()(double a, double x) const
    {
      return(a * operator()(x));
    };

    static
    double mellin_mom(int n)
    {
      // Plus part transforms to:
      // -S[1,1,n-1]
      return(-S11(n-1));
    };
    
    static
    double mellin_mom(double a, int n)
    {
      return(a * mellin_mom(n));
    };
};

// Version of the plus part with explicit eval_plus_int
class mellin_example_plus_epi : public mellin_example_plus
{
  public:
    using has_eval_plus_int = std::true_type;

    double eval_plus_int(double x) const
    {
      return(0.5*std::pow(std::log(1-x),2));
    };
};

class mellin_example_plus_epi2 : public mellin_example_plus
{
  public:
    using has_eval_plus_int = std::true_type;

    double eval_plus_int(double a, double x) const
    {
      return(a * 0.5*std::pow(std::log(1-x),2));
    };
};

class mellin_example_delta
{
  public:
    static constexpr double value = 42.;

    double operator()(double a) const
    {
      return(a * value);
    };

    static
    double mellin_mom(int n)
    {
      // Delta part transforms to:
      // 42
      return(value);
    };

    static
    double mellin_mom(double a, int n)
    {
      return(a * mellin_mom(n));
    };
};

using rpd_distribution_unary = rpd_distribution<mellin_example_reg,
      mellin_example_plus, double>;
using rpd_distribution_binary = rpd_distribution<mellin_example_reg,
      mellin_example_plus, mellin_example_delta>;

using rpd_distribution_unary_epi = rpd_distribution<mellin_example_reg,
      mellin_example_plus_epi, double>;
using rpd_distribution_binary_epi = rpd_distribution<mellin_example_reg,
      mellin_example_plus_epi, mellin_example_delta>;

using rpd_distribution_unary_epi2 = rpd_distribution<mellin_example_reg,
      mellin_example_plus_epi2, double>;
using rpd_distribution_binary_epi2 = rpd_distribution<mellin_example_reg,
      mellin_example_plus_epi2, mellin_example_delta>;

// Test cases for mellin_moment

class MellinMomentTest : public testing::Test
{
  protected:
    integration_engine_gsl engine;

    double eps_abs = 0.;
    double eps_rel = 1.e-12;
    double a = 23.;
};

TEST_F(MellinMomentTest, IntegrateWithAll)
{
  mellin_moment<double> mom(
    mellin_example_reg(),
    mellin_example_plus(),
    mellin_example_delta::value,
    engine
  );

  for(int n = 1; n <= 10; ++n)
  {
    double ref_val = mellin_example_reg::mellin_mom(n)
      +mellin_example_plus::mellin_mom(n)
      +mellin_example_delta::mellin_mom(n);
    auto res = mom.integrate(n, eps_abs, eps_rel);

    EXPECT_INT_RESULT(res, ref_val, eps_rel);
  }
}

TEST_F(MellinMomentTest, IntegrateWithOnlyReg)
{
  mellin_moment<double> mom(
    mellin_example_reg(),
    std::nullopt,
    std::nullopt,
    engine
  );

  for(int n = 1; n <= 10; ++n)
  {
    double ref_val = mellin_example_reg::mellin_mom(n);
    auto res = mom.integrate(n, eps_abs, eps_rel);

    EXPECT_INT_RESULT(res, ref_val, eps_rel);
  }
}

TEST_F(MellinMomentTest, IntegrateWithOnlyPlus)
{
  mellin_moment<double> mom(
    std::nullopt,
    mellin_example_plus(),
    std::nullopt,
    engine
  );

  for(int n = 1; n <= 10; ++n)
  {
    double ref_val = mellin_example_plus::mellin_mom(n);
    auto res = mom.integrate(n, eps_abs, eps_rel);

    EXPECT_INT_RESULT(res, ref_val, eps_rel);
  }
}

TEST_F(MellinMomentTest, IntegrateWithOnlyDelta)
{
  mellin_moment<double> mom(
    std::nullopt,
    std::nullopt,
    mellin_example_delta::value,
    engine
  );

  for(int n = 1; n <= 10; ++n)
  {
    double ref_val = mellin_example_delta::mellin_mom(n);
    auto res = mom.integrate(n, eps_abs, eps_rel);

    EXPECT_INT_RESULT(res, ref_val, eps_rel);
  }
}

TEST_F(MellinMomentTest, MakeMellinMomentUnaryWithAll)
{
  rpd_distribution_unary rpd(
    std::make_optional(mellin_example_reg()),
    std::make_optional(mellin_example_plus()),
    std::make_optional(mellin_example_delta::value)
  );
  auto mom = make_mellin_moment(rpd, engine);

  double ref_val = mellin_example_reg::mellin_mom(2)
    +mellin_example_plus::mellin_mom(2)
    +mellin_example_delta::mellin_mom(2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentUnaryWithOnlyReg)
{
  rpd_distribution_unary rpd(
    std::make_optional(mellin_example_reg()),
    std::nullopt,
    std::nullopt
  );
  auto mom = make_mellin_moment(rpd, engine);

  double ref_val = mellin_example_reg::mellin_mom(2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentUnaryWithOnlyPlus)
{
  rpd_distribution_unary rpd(
    std::nullopt,
    std::make_optional(mellin_example_plus()),
    std::nullopt
  );
  auto mom = make_mellin_moment(rpd, engine);

  double ref_val = mellin_example_plus::mellin_mom(2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentUnaryWithOnlyDelta)
{
  rpd_distribution_unary rpd(
    std::nullopt,
    std::nullopt,
    std::make_optional(mellin_example_delta::value)
  );
  auto mom = make_mellin_moment(rpd, engine);

  double ref_val = mellin_example_delta::mellin_mom(2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentBinaryWithAll)
{
  rpd_distribution_binary rpd(
    std::make_optional(mellin_example_reg()),
    std::make_optional(mellin_example_plus()),
    std::make_optional(mellin_example_delta())
  );
  auto mom = make_mellin_moment(rpd, engine, a);

  double ref_val = mellin_example_reg::mellin_mom(a, 2)
    +mellin_example_plus::mellin_mom(a, 2)
    +mellin_example_delta::mellin_mom(a, 2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentBinaryWithOnlyReg)
{
  rpd_distribution_binary rpd(
    std::make_optional(mellin_example_reg()),
    std::nullopt,
    std::nullopt
  );
  auto mom = make_mellin_moment(rpd, engine, a);

  double ref_val = mellin_example_reg::mellin_mom(a, 2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentBinaryWithOnlyPlus)
{
  rpd_distribution_binary rpd(
    std::nullopt,
    std::make_optional(mellin_example_plus()),
    std::nullopt
  );
  auto mom = make_mellin_moment(rpd, engine, a);

  double ref_val = mellin_example_plus::mellin_mom(a, 2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinMomentTest, MakeMellinMomentBinaryWithOnlyDelta)
{
  rpd_distribution_binary rpd(
    std::nullopt,
    std::nullopt,
    std::make_optional(mellin_example_delta())
  );
  auto mom = make_mellin_moment(rpd, engine, a);

  double ref_val = mellin_example_delta::mellin_mom(a, 2);

  auto res = mom.integrate(2, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}


// Example test function for convolution
double testfunc_example(double x)
{
  return(1./(1.+x));
}

class mellin_example_conv
{
  public:
    // Mellin convolution of regular part with test function examples
    static
    double res_reg_testfunc(double x)
    {
      // Mellin convolution of reg_example (x) testfunc_example
      // = M^{-1}[M[1/(1+x)] * M[1/(1+x)]]
      return(std::log(4.*x/((1.+x)*(1.+x)))/(x-1.));
    };

    static
    double res_reg_testfunc(double a, double x)
    {
      return(a * res_reg_testfunc(x));
    };

    // Mellin convolution of plus part with test function examples
    static
    double res_plus_testfunc(double x)
    {
      // Mellin convolution of [plus_example]_+ (x) testfunc_example
      // = M^{-1}[M[(1/(-log(1-x)/(1-x))_+] * M[1/(1+x)]]
      double z2 = gsl_sf_dilog(1);
      double ln2 = std::log(2);
      double ln_x = std::log(x);
      double ln_omx = std::log(1.-x);
      double ln_opx = std::log(1.+x);
      double li2_x = gsl_sf_dilog(x);
      double li2_mx = gsl_sf_dilog(-x);
      double li2_opxh = gsl_sf_dilog((1.+x)*0.5);

      return((0.5*(ln2*ln2 + z2 - ln_omx*ln_omx - ln_opx*ln_opx)
             + li2_mx + li2_opxh - li2_x + ln_opx*ln_x)/(1.+x));
    };

    static
    double res_plus_testfunc(double a, double x)
    {
      return(a * res_plus_testfunc(x));
    };

    // Mellin convolution of delta part with test function examples
    static
    double res_delta_testfunc(double x)
    {
      // Mellin convolution of \delta(1-x)*delta_example (x) testfunc_example
      // = M^{-1}[M[42 \delta(1-x)] * M[1/(1+x)]]
      return(mellin_example_delta::value/(1.+x));
    };

    static    
    double res_delta_testfunc(double a, double x)
    {
      return(a * res_delta_testfunc(x));
    };

    // Mellin convolution of regular, plus and delta part with test function
    // examples
    static
    double res_rpd_testfunc(double x)
    {
      return(res_reg_testfunc(x) + res_plus_testfunc(x) + res_delta_testfunc(x));
    };

    static
    double res_rpd_testfunc(double a, double x)
    {
      return(a * res_rpd_testfunc(x));
    };
};


// Test cases for mellin_convolution

class MellinConvolutionTest : public testing::Test
{
  protected:
    double eps_abs = 0.0;
    double eps_rel = 1.e-12;
    double x = 0.7;
    double a = 23.;

    integration_engine_gsl engine;
};

TEST_F(MellinConvolutionTest, IntegrateWithAllWithoutExplicitExtraPlusInt)
{
  mellin_convolution<double> conv(
    mellin_example_reg(),
    mellin_example_plus(),
    std::nullopt,
    mellin_example_delta::value,
    testfunc_example,
    engine
  );

  double ref_val = mellin_example_conv::res_rpd_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, IntegrateWithAllWithExplicitExtraPlusInt)
{
  mellin_convolution<double> conv(
    mellin_example_reg(),
    mellin_example_plus_epi(),
    std::bind(&mellin_example_plus_epi::eval_plus_int,
      mellin_example_plus_epi(), std::placeholders::_1),
    mellin_example_delta::value,
    testfunc_example,
    engine
  );

  double ref_val = mellin_example_conv::res_rpd_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, IntegrateWithOnlyReg)
{
  mellin_convolution<double> conv(
    mellin_example_reg(),
    std::nullopt,
    std::nullopt,
    std::nullopt,
    testfunc_example,
    engine
  );

  double ref_val = mellin_example_conv::res_reg_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, IntegrateWithOnlyPlusWithoutExplicitExtraPlusInt)
{
  mellin_convolution<double> conv(
    std::nullopt,
    mellin_example_plus(),
    std::nullopt,
    std::nullopt,
    testfunc_example,
    engine
  );

  double ref_val = mellin_example_conv::res_plus_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_EQ(std::get<0>(res), integration_status::success);
  EXPECT_NEAR(std::get<1>(res), ref_val, std::abs(ref_val) * eps_rel);
  EXPECT_LE(std::get<2>(res), std::abs(ref_val) * eps_rel);
}

TEST_F(MellinConvolutionTest, IntegrateWithOnlyPlusWithExplicitExtraPlusInt)
{
  mellin_convolution<double> conv(
    std::nullopt,
    mellin_example_plus_epi(),
    std::bind(&mellin_example_plus_epi::eval_plus_int,
      mellin_example_plus_epi(), std::placeholders::_1),
    std::nullopt,
    testfunc_example,
    engine
  );

  double ref_val = mellin_example_conv::res_plus_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, IntegrateWithOnlyDelta)
{
  mellin_convolution<double> conv(
    std::nullopt,
    std::nullopt,
    std::nullopt,
    mellin_example_delta::value,
    testfunc_example,
    engine
  );

  double ref_val = mellin_example_conv::res_delta_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}


TEST_F(MellinConvolutionTest, MakeMellinConvolutionUnaryWithAllWithoutExplicitExtraPlusInt)
{
  rpd_distribution_unary rpd(
    std::make_optional(mellin_example_reg()),
    std::make_optional(mellin_example_plus()),
    std::make_optional(mellin_example_delta::value)
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine);

  double ref_val = mellin_example_conv::res_rpd_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionUnaryWithAllWithExplicitExtraPlusInt)
{
  rpd_distribution_unary_epi rpd(
    std::make_optional(mellin_example_reg()),
    std::make_optional(mellin_example_plus_epi()),
    std::make_optional(mellin_example_delta::value)
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine);

  double ref_val = mellin_example_conv::res_rpd_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionUnaryWithOnlyReg)
{
  rpd_distribution_unary rpd(
    std::make_optional(mellin_example_reg()),
    std::nullopt,
    std::nullopt
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine);

  double ref_val = mellin_example_conv::res_reg_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionUnaryWithOnlyPlusWithoutExplicitExtraPlusInt)
{
  rpd_distribution_unary rpd(
    std::nullopt,
    std::make_optional(mellin_example_plus()),
    std::nullopt
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine);

  double ref_val = mellin_example_conv::res_plus_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionUnaryWithOnlyDelta)
{
  rpd_distribution_unary rpd(
    std::nullopt,
    std::nullopt,
    std::make_optional(mellin_example_delta::value)
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine);

  double ref_val = mellin_example_conv::res_delta_testfunc(x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionBinaryWithAllWithoutExplicitExtraPlusInt)
{
  rpd_distribution_binary rpd(
    std::make_optional(mellin_example_reg()),
    std::make_optional(mellin_example_plus()),
    std::make_optional(mellin_example_delta())
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine, a);

  double ref_val = mellin_example_conv::res_rpd_testfunc(a, x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionBinaryWithAllWithExplicitExtraPlusInt)
{
  rpd_distribution_binary_epi2 rpd(
    std::make_optional(mellin_example_reg()),
    std::make_optional(mellin_example_plus_epi2()),
    std::make_optional(mellin_example_delta())
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine, a);

  double ref_val = mellin_example_conv::res_rpd_testfunc(a, x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionBinaryWithOnlyReg)
{
  rpd_distribution_binary rpd(
    std::make_optional(mellin_example_reg()),
    std::nullopt,
    std::nullopt
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine, a);

  double ref_val = mellin_example_conv::res_reg_testfunc(a, x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionBinaryWithOnlyPlusWithoutExplicitExtraPlusInt)
{
  rpd_distribution_binary rpd(
    std::nullopt,
    std::make_optional(mellin_example_plus()),
    std::nullopt
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine, a);

  double ref_val = mellin_example_conv::res_plus_testfunc(a, x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionBinaryWithOnlyPlusWithExplicitExtraPlusInt)
{
  rpd_distribution_binary_epi2 rpd(
    std::nullopt,
    std::make_optional(mellin_example_plus_epi2()),
    std::nullopt
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine, a);

  double ref_val = mellin_example_conv::res_plus_testfunc(a, x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}

TEST_F(MellinConvolutionTest, MakeMellinConvolutionBinaryWithOnlyDelta)
{
  rpd_distribution_binary rpd(
    std::nullopt,
    std::nullopt,
    std::make_optional(mellin_example_delta())
  );
  auto conv = make_mellin_convolution(rpd, std::function(testfunc_example), engine, a);

  double ref_val = mellin_example_conv::res_delta_testfunc(a, x);

  auto res = conv.integrate(x, eps_abs, eps_rel);

  EXPECT_INT_RESULT(res, ref_val, eps_rel);
}
