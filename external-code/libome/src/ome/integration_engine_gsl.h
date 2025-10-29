/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Integration engine based on the GNU Scientific Libraray
 */

#ifndef LIBOME_INTEGRATION_ENGINE_GSL
#define LIBOME_INTEGRATION_ENGINE_GSL

#include <tuple>
#include <vector>
#include <ome/integration_engine.h>

namespace ome
{

  /**
   * \brief Double precision integration engine based on the GSL
   *
   * \details
   * This engine uses the CQUAD implementation of the GNU Scientific
   * Library.
   */
  class integration_engine_gsl : public integration_engine<double>
  {
    public:
      /**
       * \brief Calculate a one-dimensional integral
       *
       * \details
       * Given a univariate callable \f$f(x)\f$ it calculates the numerical
       * integral
       * \f[
       *   \int\limits_{x_{\text{min}}}^{x_{\text{max}}} f(x) \, \mathrm{d}x
       * \f]
       *
       * \param f Univariate callable that represents tha function to integrate
       * \param x_min Lower integration boundary
       * \param x_max Upper integration boundary
       * \param eps_abs Requested absolute error
       * \param eps_rel Requested relative error
       *
       * \returns A tuple of the status of the integration, the result of the
       *          integral and an estimate for the absolute error on the result.
       */
      std::tuple<integration_status, numeric_type, numeric_type>
      integrate(function_type f,
                numeric_type x_min,
                numeric_type x_max,
                numeric_type eps_abs,
                numeric_type eps_rel) const override;

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
                    numeric_type eps_rel) const override;

    private:
      /**
       * \brief GSL interand helper function
       *
       * \details
       * The function follows the conventions for integrands of GSL integrators
       * and takes a pointer to a std::function<double(double)> as a void pointer.
       * It is cast back to its origial type and evaluated on the first numeric
       * argument. The return value of the std::function is returned.
       *
       * \param x Point where the integrand should be evaluated
       * \param params void pointer to std::function<double(double)> that actually
       *        evaluates the integrand function
       *
       * \return Result of evaluating the integrand at x
       */
      static
      double gsl_helper_kernel(double x, void* params);
  };
}

#endif
