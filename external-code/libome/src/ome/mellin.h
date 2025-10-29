/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Mellin transformations and convolutions
 */

#ifndef LIBOME_MELLIN_H
#define LIBOME_MELLIN_H

#include <optional>
#include <functional>
#include <tuple>
#include <vector>
#include <cmath>
#include <ome/traits.h>
#include <ome/integration_engine.h>
#include <ome/rpd_distribution.h>

namespace ome
{
  /**
   * \brief Mellin transformation
   *
   * \details
   * Given univariate callables for the regular and plus parts
   * (\f$f_\mathrm{reg}(x)\f$ and \f$f_+(x)\f$) and the numerical value of the
   * delta part (\f$f_\delta\f$), als well as an integration eninge, this
   * class implements the calculation of Mellin moments, defined via
   * \f[
   *   \int_0^1 \mathrm{d}x \, x^{n-1} \left(f_\mathrm{reg}(x)
   *     +\left[f_+(x)\right]_+
   *     +\delta(1-x) f_\delta\right)
   *   = \int_0^1 \mathrm{d}x \left(x^{n-1} f_\mathrm{reg}(x)
   *     +(x^{n-1} - 1) f_+(x)\right) + f_\delta
   * \f]
   * The presence of the regular, plus and delta parts is optional.
   *
   * \tparam Tnum Numerical type
   */
  template<typename Tnum>
  class mellin_moment
  {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the univariate function type
      using function_type = std::function<numeric_type(numeric_type)>;

      /**
       * \brief Construct Mellin moment from a regular, plus and delta part
       *
       * \details
       * The individual parts are wrapped in std::optional so that their absence
       * can be signalled by passing std::nullopt.
       *
       * \param regular_part Regular part \f$f_\mathrm{reg}(x)\f$ (optional
       *        univariate callable)
       * \param plus_part Plus part \f$f_+(x)\f$ (optional univariate callable)
       * \param delta_part Delta part \f$f_\delta\f$ (optional numerical value)
       * \param engine Integration engine to use for numerical integrals
       */
      mellin_moment(std::optional<function_type> regular_part,
                    std::optional<function_type> plus_part,
                    std::optional<numeric_type> delta_part,
                    const integration_engine<numeric_type>& engine)
        : regular_part_(regular_part),
          plus_part_(plus_part),
          delta_part_(delta_part),
          engine_(engine) {};

      /**
       * \brief Compute the \f$n\f$th Mellin moment
       *
       * \details
       * Compute the \f$n\f$th Mellin moment via numerical integration to a
       * given precision. The absolute and relative precision goal are passed
       * on to the integration engine.
       *
       * \param n Mellin moment
       * \param eps_abs Absolute precision goal
       * \param eps_rel Relative precision goal
       *
       * \returns A tuple of the "last" integration status, the result for the
       *          Mellin moment and an uncertainty estimate.
       */
      std::tuple<integration_status, numeric_type, numeric_type>
        integrate(int n, numeric_type eps_abs, numeric_type eps_rel)
      {
        int nm1 = n-1;

        std::vector<std::tuple<function_type, numeric_type, numeric_type>>
          integrals;

        numeric_type offset = static_cast<numeric_type>(0);

        if(regular_part_.has_value())
        {
          auto reg_func = *regular_part_;
          integrals.emplace_back(
            [reg_func, nm1] (numeric_type x)
            {
              using std::pow;
              return(pow(x, nm1) * reg_func(x));
            },
            static_cast<numeric_type>(0),
            static_cast<numeric_type>(1)
          );
        }
        if(plus_part_.has_value())
        {
          auto plus_func = *plus_part_;
          integrals.emplace_back(
            [plus_func, nm1] (numeric_type x)
            {
              using std::pow;
              return((pow(x, nm1) - static_cast<numeric_type>(1))
                     * plus_func(x));
            },
            static_cast<numeric_type>(0),
            static_cast<numeric_type>(1)
          );
        }
        if(delta_part_.has_value())
        {
          offset += *delta_part_;
        }

        return(engine_.integrate_sum(integrals, offset, eps_abs, eps_rel));
      };

    private:
      std::optional<function_type> regular_part_;
      std::optional<function_type> plus_part_;
      std::optional<numeric_type> delta_part_;

      const integration_engine<numeric_type>& engine_;
  };


  /**
   * \brief Helper function to construct a mellin_moment object from an
   *        rpd_distribution
   *
   * \details
   * The rpd_distribution container holds the regular, plus and delta parts of
   * a function, but the integrand for numerically computing Mellin moments
   * requires a univariate function of the integration variable.
   * This helper function takes values for the parameters that are not
   * integrated over and constructs univariate functions by binding the free
   * arguments and passes them on to construct a mellin_moment object.
   *
   * \tparam Tnum Numerical type
   * \tparam Tfuncreg Callable type for the regular part
   * \tparam Tfuncplus Callable type for the plus part
   * \tparam Tfuncdelta Callable type for the delta part
   * \tparam Tcoeffarg Type parameter pack for the argument that are not
   *         integrated over
   *
   * \param rpd container for the mixed distribution
   * \param engine Integration engine to use for numerical integrals
   * \param args Function parameter pack that specifies the values for the
   *        remaining argument that are not integrated over
   *
   * \return Constructed mellin_moment object
   */
  template<typename Tnum,
           typename Tfuncreg,
           typename Tfuncplus,
           typename Tfuncdelta,
           typename... Tcoeffargs>
  mellin_moment<Tnum>
  make_mellin_moment(
    const rpd_distribution<Tfuncreg, Tfuncplus, Tfuncdelta>& rpd,
    const integration_engine<Tnum>& engine,
    Tcoeffargs... args
  )
  {
    using namespace std::placeholders;

    // We have to distinguish the case where there are additional args
    // from the one that does not since if there are no additonal args
    // the delta part is just a numerical value and not a callable.
    if constexpr (sizeof...(args) > 0)
    {
      return(mellin_moment<Tnum>(
        rpd.has_regular() ?
          std::make_optional(std::bind(*rpd.get_regular(), args..., _1))
          : std::nullopt,
        rpd.has_plus() ?
          std::make_optional(std::bind(*rpd.get_plus(), args..., _1))
          : std::nullopt,
        rpd.has_delta() ?
          std::make_optional((*rpd.get_delta())(args...))
          : std::nullopt,
        engine
      ));
    }
    else
    {
      return(mellin_moment<Tnum>(
        rpd.get_regular(),
        rpd.get_plus(),
        rpd.get_delta(),
        engine
      ));
    }
  }


  /**
   * \brief Mellin convolution
   *
   * \details
   * Given univariate callables for the regular and plus parts 
   * (\f$f_\mathrm{reg}(x)\f$ and \f$f_+(x)\f$) and the numerical value of the
   * delta part (\f$f_\delta\f$), a test function (\f$g(x)\f$) as well as an
   * integration eninge, this class implements the Mellin convolution, defined
   * via
   * \f[
   *   \left(f_\mathrm{reg}(x) + [f_+(x)]_+ + \delta(1-x) f_\delta\right)
   *     \otimes g(x)
   *     = \int_x^1 \frac{\mathrm{d}z}{z} f_\mathrm{reg}(x)
   *                g\left(\frac{x}{z}\right)
   *       +\int_x^1 \mathrm{d}z \, f_+(z) \left[
   *         \frac{1}{z} g\left(\frac{x}{z}\right) - g(x)
   *       \right]
   *       -\int_0^x \mathrm{d} \, f_+(z) g(x)
   *       +f_\delta g(x)
   * \f]
   * The presence of the regular, plus and delta parts is optional.
   *
   * \tparam Tnum Numerical type
   */
  template<typename Tnum>
  class mellin_convolution
  {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the univariate function type
      using function_type = std::function<numeric_type(numeric_type)>;

      /**
       * \brief Construct Mellin convolution from a regular, plus and delta part
       *        as well as the test function
       *
       * \details
       * The regular, plus and delta parts are wrapped in std::optional so that
       * their absence can be signalled by passing std::nullopt.
       * In addition, it is possible to pass a function as plus_part_extra which
       * directly evaluates to the result of the integral \f$\int_0^x \mathrm{d}z
       * \, f_+(z)\f$. If this function is not given but plus_part is present,
       * the integral is calculated numerically.
       *
       * \param regular_part Regular part \f$f_\mathrm{reg}(x)\f$ (optional
       *        univariate callable)
       * \param plus_part Plus part \f$f_+(x)\f$ (optional univariate callable)
       * \param plus_part_extra Integral of the plus part over \f$[0,x]\f$
       *        (optional univariate callable)
       * \param delta_part Delta part \f$f_\delta\f$ (optional numerical value)
       * \param testfunc Test function to convolve with
       * \param engine Integration engine to use for numerical integrals
       */
      mellin_convolution(std::optional<function_type> regular_part,
                         std::optional<function_type> plus_part,
                         std::optional<function_type> plus_part_extra,
                         std::optional<numeric_type> delta_part,
                         function_type testfunc,
                         const integration_engine<numeric_type>& engine)
        : regular_part_(regular_part),
          plus_part_(plus_part),
          plus_part_extra_(plus_part_extra),
          delta_part_(delta_part),
          testfunc_(testfunc),
          engine_(engine) {};

      /**
       * \brief Compute the Mellin convolution
       *
       * \details
       * Compute the Mellin convolution via numerical integration to a given
       * precision. The absolute and relative precision goal are passed
       * on to the integration engine. If there are both regular and plus parts
       * present, they are integrated separately and their error estimates are
       * added in quadrature. The returned integration status is indicates
       * success if all numerical integrations returned success. If one of them
       * fails, the returned status reflects the return value of the failed
       * integration. If more than one numerical integration fails, only the
       * last failure is returned.
       *
       * \param x Argument of the Mellin convolution
       * \param eps_abs Absolute precision goal
       * \param eps_rel Relative precision goal
       *
       * \returns A tuple of the "last" integration status, the result for the
       *          Mellin convolution and an uncertainty estimate.
       */
      std::tuple<integration_status, numeric_type, numeric_type>
        integrate(numeric_type x, numeric_type eps_abs, numeric_type eps_rel)
      {
        std::vector<std::tuple<function_type, numeric_type, numeric_type>>
          integrals;

        numeric_type offset = static_cast<numeric_type>(0);

        auto& testfunc = testfunc_;
        if(regular_part_.has_value())
        {
          auto reg_func = *regular_part_;
          integrals.emplace_back(
            [reg_func, testfunc, x] (numeric_type z)
            {
              return(reg_func(z) * testfunc(x/z)/z);
            },
            x,
            static_cast<numeric_type>(1)
          );
        }
        if(plus_part_.has_value())
        {
          // First part of the plus integral (over [x,1])
          auto plus_func = *plus_part_;
          integrals.emplace_back(
            [plus_func, testfunc, x] (numeric_type z)
            {
              return(plus_func(z) * (testfunc(x/z)/z - testfunc(x)));
            },
            x,
            static_cast<numeric_type>(1)
          );

          // Second part of the plus integral (over [0,x])
          if(plus_part_extra_.has_value())
          {
            offset -= (*plus_part_extra_)(x) * testfunc_(x);
          }
          else
          {
            numeric_type testfunc_eval = testfunc_(x);
            integrals.emplace_back(
              [plus_func, testfunc_eval] (numeric_type x)
              {
                return(-plus_func(x) * testfunc_eval);
              },
              static_cast<numeric_type>(0),
              x
            );
          }
        }
        if(delta_part_.has_value())
        {
          offset += *delta_part_ * testfunc_(x);
        }

        return(engine_.integrate_sum(integrals, offset, eps_abs, eps_rel));
      };

    private:
      std::optional<function_type> regular_part_;
      std::optional<function_type> plus_part_;
      std::optional<function_type> plus_part_extra_;
      std::optional<numeric_type> delta_part_;

      function_type testfunc_;

      const integration_engine<numeric_type>& engine_;
  };

  /**
   * \brief Helper function to construct a mellin_convolution object from an
   *        rpd_distribution
   *
   * \details
   * The rpd_distribution container holds the regular, plus and delta parts of
   * a function, but the integrand for numerically computing Mellin moments
   * requires a univariate function of the integration variable.
   * This helper function takes values for the parameters that are not
   * integrated over and constructs univariate functions by binding the free
   * arguments and passes them on to construct a mellin_convolution object.
   *
   * \tparam Tnum Numerical type
   * \tparam Tfuncreg Callable type for the regular part
   * \tparam Tfuncplus Callable type for the plus part
   * \tparam Tfuncdelta Callable type for the delta part
   * \tparam Tcoeffarg Type parameter pack for the argument that are not
   *         integrated over
   *
   * \param rpd container for the mixed distribution
   * \param testfunc Function with which to convolve
   * \param engine Integration engine to use for numerical integrals
   * \param args Function parameter pack that specifies the values for the
   *        remaining argument that are not integrated over
   *
   * \return Constructed mellin_moment object
   */

  template<typename Tnum,
           typename Tfuncreg,
           typename Tfuncplus,
           typename Tfuncdelta,
           typename... Tcoeffargs>
  mellin_convolution<Tnum>
  make_mellin_convolution(
    const rpd_distribution<Tfuncreg, Tfuncplus, Tfuncdelta>& rpd,
    std::function<Tnum(Tnum)> testfunc,
    const integration_engine<Tnum>& engine,
    Tcoeffargs... args
  )
  {
    using namespace std::placeholders;
    std::optional<std::function<Tnum(Tnum)>> eval_plus_int;

    // We have to distinguish the case where there are additional args
    // from the one that does not since if there are no additonal args
    // the delta part is just a numerical value and not a callable.
    if constexpr (sizeof...(args) > 0)
    {
      if constexpr(has_eval_plus_int<Tfuncplus>::value)
      {
        if(rpd.has_plus())
        {
          auto plus_part = *rpd.get_plus();
          eval_plus_int = std::make_optional(std::bind(&decltype(plus_part)::eval_plus_int, plus_part, args..., _1));
        }
      }
      return(mellin_convolution<Tnum>(
        rpd.has_regular() ?
          std::make_optional(std::bind(*rpd.get_regular(), args..., _1))
          : std::nullopt,
        rpd.has_plus() ?
          std::make_optional(std::bind(*rpd.get_plus(), args..., _1))
          : std::nullopt,
        eval_plus_int,
        rpd.has_delta() ?
          std::make_optional((*rpd.get_delta())(args...))
          : std::nullopt,
        testfunc,
        engine
      ));
    }
    else
    {
      if constexpr(has_eval_plus_int<Tfuncplus>::value)
      {
        if(rpd.has_plus())
        {
          auto plus_part = *rpd.get_plus();
          eval_plus_int = std::make_optional(std::bind(&decltype(plus_part)::eval_plus_int, plus_part, _1));
        }
      }
      return(mellin_convolution<Tnum>(
        rpd.get_regular(),
        rpd.get_plus(),
        eval_plus_int,
        rpd.get_delta(),
        testfunc,
        engine
      ));
    }
  }

}

#endif
