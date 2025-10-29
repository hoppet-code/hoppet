/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Abstract interface to different integration engines
 */

#ifndef LIBOME_INTEGRATION_ENGINE_H
#define LIBOME_INTEGRATION_ENGINE_H

#include <functional>
#include <tuple>
#include <vector>

namespace ome
{
  /**
   * \brief Status codes for the integration routines
   */
  enum class integration_status
  {
    /// Integration successful 
    success,
    /// Adaptive integration reached max number of subdivisions
    max_iterations_reached,
    /// Requested precision could not be reached due to rounding error
    rounding_error_detected,
    /// A non-integrable singularity or other bad integrand behavior was found
    /// in the integration interval
    singularity_detected,
    /// The integral is divergent, or too slowly convergent to be integrated
    /// numerically
    divergence_detected,
    /// The requested precision could not be reached
    precision_not_reached,
    /// One of the arguments to the integrator had an invalid value
    domain_error,
    /// Bad tolerace specification
    bad_tolerance,
    /// Unknown error from the GSL integrator
    other_error
  };


  /**
   * \brief Interface to different integration engines
   *
   * \details
   * An integration engine performes one dimensional integrals over univariate
   * functions. This interface unifies different possible implementations.
   *
   * \tparam Tnum Numerical type
   */
  template<typename Tnum>
  class integration_engine
  {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the univariate function type
      using function_type = std::function<numeric_type(numeric_type)>;

      /**
       * \brief Calculate a one-dimensional integral
       *
       * \details
       * Given a univariate callable \f$f(x)\f$ it calculate the numerical
       * integral
       * \f[
       *   \int\limits_{x_{\text{min}}}^{x_{\text{max}}} f(x) \, \mathrm{d}x
       * \f]
       *
       * \param f Univariate callable that represents the function to integrate
       * \param x_min Lower integration boundary
       * \param x_max Upper integration boundary
       * \param eps_abs Requested absolute error
       * \param eps_rel Requested relative error
       *
       * \returns A tuple of the status of the integration, the result of the
       *          integral and an estimate for the absolute error on the result.
       */
      virtual std::tuple<integration_status, numeric_type, numeric_type>
      integrate(function_type f,
                numeric_type x_min,
                numeric_type x_max,
                numeric_type eps_abs,
                numeric_type eps_rel) const = 0;

      /**
       * \brief Calculate a sum of integrals with a common relative error goal
       *
       * \details
       * Given tuples of univariate callables \f$f_i(x)\f$, and lower and upper
       * integration boundaries, \f$x_{\text{min},i}\f$ and
       * \f$x_{\text{max},i}\f$ as well as an offset \f$y\f$, calculate the sum
       * of numerical integrals
       * \f[
       *   y + \sum_i \int\limits_{x_{\text{min},i}}^{x_{\text{max},i}} f_i(x)
       *     \, \mathrm{d}x
       * \f]
       * while targeting a relative error goal for the sum of the integrals and
       * the offset.
       *
       * \param integrals Vector tuples of univariate callables that represent
       *        the functions to integrate, and lower and upper integration
       *        boundaries.
       * \param offset The offset \f$y\f$
       * \param eps_abs Requested absolute error
       * \param eps_rel Requested relative error
       *
       * \returns A tuple of the status of the integration, the result of the
       *          integral and an estimate for the absolute error on the result.
       */
      virtual std::tuple<integration_status, numeric_type, numeric_type>
      integrate_sum(std::vector<std::tuple<function_type, numeric_type,
                      numeric_type>> integrals,
                    numeric_type offset,
                    numeric_type eps_abs,
                    numeric_type eps_rel) const = 0;

      /// Virtual destructor
      virtual ~integration_engine() {};
  };

}

#endif

