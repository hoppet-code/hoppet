/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <utility>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <ome/integration_engine.h>
#include <ome/integration_engine_gsl.h>

namespace ome
{

  using numeric_type = integration_engine_gsl::numeric_type;
  using function_type = integration_engine_gsl::function_type;

  std::tuple<integration_status, numeric_type, numeric_type>
  integration_engine_gsl::integrate(function_type f,
            numeric_type x_min,
            numeric_type x_max,
            numeric_type eps_abs,
            numeric_type eps_rel) const
  {
    // Max. number of subintervals into which the integration range is split up
    constexpr size_t max_subintervals = 8*1024;
    // Result and absolute error of the integral
    double res = 0., res_abs_err = 0.;

    // Prepare integrand
    gsl_function integrand;
    integrand.function = &gsl_helper_kernel;
    integrand.params = &f;

    // Allocate integration workspace
    gsl_integration_cquad_workspace* workspace
      = gsl_integration_cquad_workspace_alloc(max_subintervals);
    if(!workspace)
      throw std::bad_alloc();

    // Run integrtion
    int ret = gsl_integration_cquad(
      &integrand,
      x_min,
      x_max,
      eps_abs,
      eps_rel,
      workspace,
      &res,
      &res_abs_err,
      nullptr
    );

    gsl_integration_cquad_workspace_free(workspace);


    // Check for errors in the integration
    if(ret != GSL_SUCCESS)
    {
      switch(ret)
      {
        case GSL_EMAXITER:
          return(std::make_tuple(integration_status::max_iterations_reached,
                                 res, res_abs_err));
        case GSL_EROUND:
          return(std::make_tuple(integration_status::rounding_error_detected,
                                 res, res_abs_err));
        case GSL_ESING:
          return(std::make_tuple(integration_status::singularity_detected,
                                 res, res_abs_err));
        case GSL_EDIVERGE:
          return(std::make_tuple(integration_status::divergence_detected,
                                 res, res_abs_err));
        case GSL_EDOM:
          return(std::make_tuple(integration_status::domain_error,
                                 res, res_abs_err));
        case GSL_EBADTOL:
          return(std::make_tuple(integration_status::bad_tolerance,
                                 res, res_abs_err));
        default:
          return(std::make_tuple(integration_status::other_error,
                                 res, res_abs_err));
      }
    }
  
    if(eps_abs < res_abs_err && (res == 0. || eps_rel < res_abs_err/fabs(res)))
    {
      return(std::make_tuple(integration_status::precision_not_reached,
                             res, res_abs_err));
    }
  
    return(std::make_tuple(integration_status::success, res, res_abs_err));        
  }


  std::tuple<integration_status, numeric_type, numeric_type>
  integration_engine_gsl::integrate_sum(std::vector<std::tuple<function_type,
                  numeric_type, numeric_type>> integrals,
                numeric_type offset,
                numeric_type eps_abs,
                numeric_type eps_rel) const
  {
    using std::sqrt;
    using std::min;
    using std::abs;

    std::vector<std::pair<numeric_type, numeric_type>> res_and_sqerrors;
    integration_status status = integration_status::success;

    if(integrals.empty())
      return(std::make_tuple(integration_status::success, offset,
                             static_cast<numeric_type>(0)));

    // First pass: try with naive eps/sqrt(n_ints)
    // This will succeed if all integrals yield the same contribution and
    // offset is zero
    numeric_type inv_sqrt_n_ints
      = static_cast<numeric_type>(1)/sqrt(integrals.size());
    numeric_type current_eps_rel = eps_rel * inv_sqrt_n_ints;
    numeric_type current_eps_abs = eps_abs * inv_sqrt_n_ints;
    for(const auto& [f, x_min, x_max] : integrals)
    {
      auto res = integrate(f, x_min, x_max,
                           current_eps_abs, current_eps_rel);
      if(std::get<0>(res) != integration_status::success)
        status = std::get<0>(res); // FIXME: only last status survives

      res_and_sqerrors.emplace_back(std::get<1>(res),
                                    std::get<2>(res)*std::get<2>(res));
    }

    // Iterative refinement
    while(true)
    {
      // Check termination criteria (error occurred or goals reached)
      numeric_type current_total = static_cast<numeric_type>(offset);
      numeric_type current_total_sq_abs_error = static_cast<numeric_type>(0);
      for(auto& [res, sq_abs_error] : res_and_sqerrors)
      {
        current_total += res;
        current_total_sq_abs_error += sq_abs_error;
      }
      numeric_type current_total_abs_error = sqrt(current_total_sq_abs_error);
      if(status != integration_status::success ||
         current_total_abs_error <= eps_abs ||
         current_total_abs_error <= eps_rel * abs(current_total))
      {
        return(std::make_tuple(status, current_total, current_total_abs_error));
      }

      // Select integral with largest error
      auto largest_error = std::max_element(
        res_and_sqerrors.cbegin(),
        res_and_sqerrors.cend(),
        [] (const auto& a, const auto& b) { return(a.second < b.second); }
      );
      auto largest_error_index
        = std::distance(res_and_sqerrors.cbegin(), largest_error);

      // Calculate refined error goal for this integral
      // We adjust the absolute error goal since we now have an estimate
      // for the total result.
      current_eps_abs = min(
        eps_rel * abs(current_total) * inv_sqrt_n_ints,
        0.5*sqrt(largest_error->second)
      );

      // Run refined integration
      auto &largest_error_integral = integrals[largest_error_index];
      auto res = integrate(
        std::get<0>(largest_error_integral),
        std::get<1>(largest_error_integral),
        std::get<2>(largest_error_integral),
        current_eps_abs,
        static_cast<numeric_type>(0)
      );
      if(std::get<0>(res) != integration_status::success)
        status = std::get<0>(res);
      res_and_sqerrors[largest_error_index]
        = std::make_pair(std::get<1>(res), std::get<2>(res) * std::get<2>(res));
    }

  }
  

  double integration_engine_gsl::gsl_helper_kernel(double x, void* params)
  {
    std::function<double(double)> *f =
      static_cast<std::function<double(double)>*>(params);
    return((*f)(x));
  }
}
