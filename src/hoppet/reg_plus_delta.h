
#include <tuple>
#include <utility>
#include <type_traits>
#include <stdexcept>
#include <functional>

namespace hoppet {

// this class has largely been written by ChatGPT based on the
// specification of how it should function and then lightly edited by
// hand.

/// @brief A class that combines regular, plus, and delta functions
/// @tparam R the type of the regular function
/// @tparam P the type of the plus function
/// @tparam D the type of the delta function
/// @tparam ...StoredArgs the types of any additional stored arguments
///
/// An object of this class can be used to initialise a
/// hoppet::grid_conv (aka hoppet::split_fn) object, providing the 
/// necessary interface `double operator()(double y, int piece)`.
///
/// Calling that function invokes the appropriate component functions
/// based on the value of `piece`, passing `x = exp(-y)` to relevant
/// combination of the regular, plus, and delta functions. Any additional
/// stored arguments are also passed to those functions if they accept them
/// (e.g. for use for an `nf` argument).
///
/// If a component function is a constant (not callable), it is treated
/// as a constant value. This is most useful when a component is absent, in
/// which case one can simply pass `0.0` as that component.
///
/// The functions should have measure dx (the extra factor of x needed
/// by the `double operator()(double y, int piece)` is provided by the
/// reg_plus_delta object).
///
/// ## Examples
///
/// This shows usage with simple lambdas for the three components,
/// assuming there is already a `grid` (hoppet::grid_def) object
/// defined, as well as a colour factor `cf`:
///
/// ```cpp
///  hoppet::split_fn pqq_v1(grid, hoppet::reg_plus_delta(
///    [](double x) {return -cf*(1+x); },
///    [](double x) {return 2*cf/(1-x); },
///    [](double x) {return cf * 3.0/2.0; }
///  ));
/// ```
///
/// Alternatively, if we express pqq in an equivalent form that involves
/// just a plus part, the regular and delta parts can simply be set to
/// zero:
///
/// ```cpp
///  hoppet::split_fn pqq_v2(grid, hoppet::reg_plus_delta(
///     0.0, 
///     [](double x) {return cf*(1+x*x)/(1-x);}, 
///     0.0
///  ));
/// ```
template<class R, class P, class D, class... StoredArgs>
class reg_plus_delta {
private:
  R reg_obj;
  P plus_obj;
  D delta_obj;

  std::tuple<StoredArgs...> stored_args;

public:

  /// @brief Constructor that accepts regular, plus, and delta functions along with additional arguments
  ///
  /// @tparam RR type of the regular function
  /// @tparam PP type of the plus function
  /// @tparam DD type of the delta function
  /// @tparam ...Args types of any additional arguments to be stored
  ///
  /// @param r regular function
  /// @param p plus function
  /// @param d delta function
  /// @param ...args additional arguments to be stored and passed to the regular, 
  ///        plus and delta functions in addition to the x argument
  ///
  template<typename RR, typename PP, typename DD, typename... Args>
  reg_plus_delta(RR&& r, PP&& p, DD&& d, Args&&... args)
      : reg_obj(std::forward<RR>(r)),
        plus_obj(std::forward<PP>(p)),
        delta_obj(std::forward<DD>(d)),
        stored_args(std::forward<Args>(args)...)
  {}

  double operator()(double y, int piece) const {
    double x = exp(-y);
    switch(piece) {
      case cc_REAL    : return x*(call(reg_obj, x) + call(plus_obj, x));
      case cc_REALVIRT: return x* call(reg_obj, x);
      case cc_VIRT    : return -x*call(plus_obj, x);
      case cc_DELTA   : return call(delta_obj, 0.0);
      default: throw std::invalid_argument("Invalid piece index");
    }
  }

private:
  template<typename T, typename... CallArgs, std::size_t... I>
  double invoke_helper(
    const T& obj,
    std::index_sequence<I...>,
    CallArgs&&... call_args
  ) const {
    // If callable, invoke with (y, stored_args...)
    if constexpr (std::is_invocable_v<T, CallArgs..., StoredArgs...>) {
        return std::invoke(obj,
                          std::forward<CallArgs>(call_args)...,
                          std::get<I>(stored_args)...);
    } else {
        // Otherwise obj is a constant → return it
        static_assert(std::is_convertible_v<T, double>,
                      "Constant must be convertible to double");
        return obj;
    }
  }

  template<typename T, typename... CallArgs>
  double call(const T& obj, CallArgs&&... args) const {
    return invoke_helper(obj,
                         std::index_sequence_for<StoredArgs...>{},
                         std::forward<CallArgs>(args)...);
  }
};

/// Deduction guide for reg_plus_delta
/// This tells the compiler:
/// 
/// When I see reg_plus_delta(r, p, d, extra...),
/// infer R = type of r, P = type of p, D = type of d, and
/// StoredArgs = pack of types of extra.
template<typename R, typename P, typename D, typename... Args>
reg_plus_delta(R, P, D, Args...)
    -> reg_plus_delta<R, P, D, Args...>;

}
