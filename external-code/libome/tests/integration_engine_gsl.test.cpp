/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <tuple>
#include <cmath>
#include <gtest/gtest.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_dilog.h>
#include <ome/integration_engine_gsl.h>

using namespace ome;

TEST(IntegrationEngineGSLTest, UseDirect)
{
  integration_engine_gsl engine;

  double x_min = -2., x_max = 2.;
  double eps_abs = 0.;
  double eps_rel = 1.e-12;
  auto res = engine.integrate(
    [] (double x) { return(x*x); },
    x_min, x_max, eps_abs, eps_rel);

  double ref_val = 16./3.;
  EXPECT_EQ(std::get<0>(res), integration_status::success);
  EXPECT_NEAR(std::get<1>(res), ref_val, ref_val * eps_rel);
  EXPECT_LE(std::get<2>(res), ref_val * eps_rel);
}

TEST(IntegrationEngineGSL, UseThroughInterface)
{
  const integration_engine<double>& engine = integration_engine_gsl();

  double x_min = -2., x_max = 2.;
  double eps_abs = 0.;
  double eps_rel = 1.e-12;
  auto res = engine.integrate(
    [] (double x) { return(x*x); },
    x_min, x_max, eps_abs, eps_rel);

  double ref_val = 16./3.;
  EXPECT_EQ(std::get<0>(res), integration_status::success);
  EXPECT_NEAR(std::get<1>(res), ref_val, ref_val * eps_rel);
  EXPECT_LE(std::get<2>(res), ref_val * eps_rel);
}

TEST(IntegrationEngineGSL, RequestUnreachablePrecision)
{
  integration_engine_gsl engine;

  double x_min = -2., x_max = 2.;
  double eps_abs = 0.;
  double eps_rel = 1.e-17;
  gsl_error_handler_t* handler = gsl_set_error_handler_off();
  auto res = engine.integrate(
    [] (double x) { return(x*x); },
    x_min, x_max, eps_abs, eps_rel);
  gsl_set_error_handler(handler);

  EXPECT_EQ(std::get<0>(res), integration_status::bad_tolerance);
}

TEST(IntegrationEngineGSL, IntegrateDivergence)
{
  integration_engine_gsl engine;

  double x_min = 0., x_max = 1.;
  double eps_abs = 0.;
  double eps_rel = 1.e-12;

  gsl_error_handler_t* handler = gsl_set_error_handler_off();
  auto res = engine.integrate(
    [] (double x) { return(1./x/x); },
    x_min, x_max, eps_abs, eps_rel);
  gsl_set_error_handler(handler);

  EXPECT_EQ(std::get<0>(res), integration_status::divergence_detected);
}

TEST(IntegrationEngineGSL, IntegrateSum)
{
  integration_engine_gsl engine;

  std::vector<std::tuple<std::function<double(double)>, double, double>> integrals{
    std::make_tuple([] (double x) { return(x*x); }, 0., 1.),
    std::make_tuple([] (double x) { return(1./(1.+x)); }, 0., 1.)
  };
  double offset = 3.;
  double eps_abs = 0.;
  double eps_rel = 1.e-10;
  auto res = engine.integrate_sum(integrals, offset, eps_abs, eps_rel);

  double ref_val = 1./3. + std::log(2.) + 3.;
  EXPECT_EQ(std::get<0>(res), integration_status::success);
  EXPECT_NEAR(std::get<1>(res), ref_val, std::abs(ref_val) * eps_rel);
  EXPECT_LE(std::get<2>(res), std::abs(ref_val) * eps_rel);
}

TEST(IntegrationEngineGSL, IntegrateSumPerfectCancellation)
{
  integration_engine_gsl engine;

  std::vector<std::tuple<std::function<double(double)>, double, double>> integrals{
    std::make_tuple([] (double x) { return(x*x); }, 0., 1.),
    std::make_tuple([] (double x) { return(-x*x); }, 0., 1.)
  };
  double offset = 0.;
  double eps_abs = 0.;
  double eps_rel = 1.e-10;
  gsl_error_handler_t* handler = gsl_set_error_handler_off();
  auto res = engine.integrate_sum(integrals, offset, eps_abs, eps_rel);
  gsl_set_error_handler(handler);

  EXPECT_EQ(std::get<0>(res), integration_status::bad_tolerance);
}
