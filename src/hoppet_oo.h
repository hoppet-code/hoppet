#ifndef __HOPPET_OO__
#define __HOPPET_ OO__
#include "hoppet.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <tuple>
#include <functional>
//#include <concepts>

// Elements to think about
// - do we separate things out into different files
// - do we provide "physics" aliases: e.g. grid_quant -> pdf, grid_conv -> split_fn
// - do we 

// Things to perhaps add
// - [ ] decide take_view, etc., operator=
// - [ ] basic checks of grid compatibility, so as to get C++ errors rather than Fortran errors
// - [ ] mvv interface for splitting functions?


/// "forward" declaration of the Fortran grid_def type;
/// note that we only ever use pointers to it, so we do not
/// need to know anything about its actual structure, which
/// is hidden in the Fortran code.
class grid_def_f;
class grid_quant_f;
class grid_quant_2d_f;
class grid_conv_f;

/// grid_def function wrappers
extern "C" {
  grid_def_f * hoppet_cxx__grid_def__new(double dy, double ymax, int order, double eps);
  grid_def_f * hoppet_cxx__grid_def__new_from_grids(const grid_def_f ** griddefs, int ngrids, bool locked);
  grid_def_f * hoppet_cxx__grid_def__new_default(double dy, double ymax, int order = -6);
  grid_def_f * hoppet_cxx__grid_def__copy(const grid_def_f * griddef);
  void   hoppet_cxx__grid_def__delete(grid_def_f ** griddef);

  int    hoppet_cxx__grid_def__ny(const grid_def_f * griddef);
  void   hoppet_cxx__grid_def__y_values(const grid_def_f * griddef, double * yvals);
  void   hoppet_cxx__grid_def__x_values(const grid_def_f * griddef, double * xvals);
  bool   hoppet_cxx__grid_def__equiv(const grid_def_f * griddef1, const grid_def_f * griddef2);
}

/// grid_quant function wrappers
extern "C" {
  grid_quant_f * hoppet_cxx__grid_quant__new(const grid_def_f * griddef);
  void   hoppet_cxx__grid_quant__delete(grid_quant_f ** gridquant);
  double hoppet_cxx__grid_quant__at_y(const grid_def_f * griddef, const double * gq_data, double y);
  double   hoppet_cxx__grid_quant__trunc_mom(const grid_def_f * griddef, const double * data, double n, const double * ymax = 0);
  double * hoppet_cxx__grid_quant__data_ptr(grid_quant_f * gridquant);

  //void   hoppet_cxx__grid_quant__set_zero(void * gridquant);
  //void   hoppet_cxx__grid_quant__copy_from(void * gridquant, void * other);
}
inline void generic_delete(grid_quant_f * ptr) {if (ptr)hoppet_cxx__grid_quant__delete(&ptr);}

/// grid_quant_2d function wrappers
extern "C" {
  grid_quant_2d_f * hoppet_cxx__grid_quant_2d__new(const grid_def_f * griddef, int size);
  void   hoppet_cxx__grid_quant_2d__delete(grid_quant_2d_f ** gridquant);
  double * hoppet_cxx__grid_quant_2d__data_ptr(grid_quant_2d_f * gridquant);

  //void   hoppet_cxx__grid_quant__set_zero(void * gridquant);
  //void   hoppet_cxx__grid_quant__copy_from(void * gridquant, void * other);
}
inline void generic_delete(grid_quant_2d_f * ptr) {if (ptr) hoppet_cxx__grid_quant_2d__delete(&ptr);}


/// grid_conv function wrappers
extern "C" {
  grid_conv_f * hoppet_cxx_grid_conv__new_from_fn(const grid_def_f * grid_ptr, void * conv_ignd_c_fn_obj);
  grid_conv_f * hoppet_cxx_grid_conv__new_from_gc(const grid_conv_f * gc_other);
  void  hoppet_cxx_grid_conv__delete(grid_conv_f ** gridconv);
  void  hoppet_cxx_grid_conv__times_grid_quant(const grid_conv_f * conv,
                   const double * q_in_data, double * q_out_data);

  void hoppet_cxx_grid_conv__add(grid_conv_f * conv, const grid_conv_f * other, const double * factor = nullptr);
  void hoppet_cxx_grid_conv__multiply(grid_conv_f * conv, const double factor);
  grid_conv_f * hoppet_cxx_grid_conv__alloc_and_conv(const grid_conv_f * conv1, const grid_conv_f * conv2);
}

namespace hoppet {

typedef std::size_t size_type;


template <typename F>
using IsDoubleFunction = std::enable_if_t<
    std::is_invocable_r_v<double, F, double>, int
>;
//template <typename F>
//concept DoubleFunction =
//    std::is_invocable<F, double> &&
//    std::same_as<std::invoke_result_t<F, double>, double>;

//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_def Fortran type, non-owning
///
/// This version provides a "view" onto an existing Fortran grid_def object,
/// without taking ownership of it.
class grid_def_view {
public:

  grid_def_view() noexcept {}

  grid_def_view(grid_def_f * ptr) noexcept : _ptr(ptr) {}
  //grid_def_view(const grid_def_view & other) noexcept : _ptr(other._ptr) {}
  //grid_def_view & operator= (const grid_def_view & other) noexcept {
  //  _ptr = other._ptr;
  //  return *this;
  //}

  int ny() const {ensure_valid(); return hoppet_cxx__grid_def__ny(_ptr); }

  std::vector<double> y_values() const {
    std::vector<double> yvals(ny()+1);
    hoppet_cxx__grid_def__y_values(_ptr, yvals.data());
    return yvals;
  }
  std::vector<double> x_values() const {
    std::vector<double> xvals(ny()+1);
    hoppet_cxx__grid_def__x_values(_ptr, xvals.data());
    return xvals;
  }

  inline bool operator==(const grid_def_view & other) const noexcept {
    // equivalence is only true if both pointers are non-null
    // they are either the same pointer or equivalent grids
    return ptr() && other.ptr() && (
      _ptr == other._ptr || hoppet_cxx__grid_def__equiv(ptr(), other.ptr())
    );
  }
  inline bool operator!=(const grid_def_view & other) const noexcept {
    return !(*this == other);
  }

  void ensure_valid() const {if (!_ptr) throw std::runtime_error("hoppet::grid_def_view::ensure_valid: grid pointer is null");}
  void ensure_compatible(const grid_def_view & other) const {
    if (*this != other) {
      throw std::runtime_error("hoppet::grid_def_view::ensure_compatible: grids are not equivalent");
    }
  }

  const grid_def_f * ptr() const { return _ptr; }
  grid_def_f * ptr() { return _ptr; }
protected:
  /// @brief pointer to the underlying Fortran grid_def object
  grid_def_f * _ptr = nullptr;
};


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_def Fortran type, with ownership
///
/// This version takes ownership of the underlying Fortran grid_def object
class grid_def : public grid_def_view {
public:

  grid_def() {};

  /// @brief construct and allocate a new grid_def object
  ///
  /// @param dy      grid spacing
  /// @param ymax    maximum y value
  /// @param order   usual interpolation order parameter
  /// @param eps     accuracy for the adaptive integration
  //
  grid_def(double dy, double ymax, int order=-5, double eps=1e-7)
    : grid_def_view(hoppet_cxx__grid_def__new(dy, ymax, order, eps)) {}

  /// @brief construct a grid_def from multiple grid_def objects
  ///
  /// @param grids   vector of grid_def objects
  /// @param locked  whether the new grid should be "locked"
  grid_def(const std::vector<grid_def> & grids, bool locked=false) {
    int ngrids = static_cast<int>(grids.size());
    std::vector<const grid_def_f *> grid_ptrs(ngrids);
    for (int i=0; i<ngrids; ++i) {
      grid_ptrs[i] = grids[i].ptr();
    }
    _ptr = hoppet_cxx__grid_def__new_from_grids(grid_ptrs.data(), ngrids, locked);
  }


  grid_def(const grid_def & other)
    : grid_def_view(hoppet_cxx__grid_def__copy(other.ptr())) {}

  grid_def(grid_def && other) noexcept
    : grid_def_view(other.ptr()) {
    // null out the other pointer to avoid double deletion
    other._ptr = nullptr;
  }

  grid_def & operator=(const grid_def & other) {
    if (this != &other) {
      if (_ptr) hoppet_cxx__grid_def__delete(&_ptr);
      _ptr = hoppet_cxx__grid_def__copy(other.ptr());
    }
    return *this;
  }
    
  ~grid_def() {if (_ptr) hoppet_cxx__grid_def__delete(&_ptr); }

protected:
  grid_def(grid_def_f * ptr) : grid_def_view(ptr) {}
  friend grid_def grid_def_default(double dy, double ymax, int order);
};  

/// @brief  create a grid_def object with the default choice of nested, locked grids
/// @param dy      spacing of the coarsest grid
/// @param ymax    maximum y value for the coarsest grid
/// @param order   usual interpolation order parameter
/// @return 
inline grid_def grid_def_default(double dy, double ymax, int order) {
  return grid_def(hoppet_cxx__grid_def__new_default(dy, ymax, order));
}

/// @brief Base class for views of objects with an associated data pointer
template<typename T>
class data_view {
public:

  typedef T extras_type;

  data_view() noexcept {}

  data_view(double * data_ptr, std::size_t size, const T & extras) noexcept
    : _data(data_ptr), _size(size), _extras(extras) {}

  explicit data_view(const data_view<T> & other) noexcept {
    take_view(other);
  }

  data_view<T> & operator=(const data_view<T> & other) {
    this->copy_data(other);
    return *this;
  }

  void take_view(const data_view<T> & other) noexcept {
    _data = other._data;
    _size = other._size;
    _extras = other._extras;
  }

  double       * data()       {return _data;}
  const double * data() const {return _data;}
  std::size_t    size() const {return _size;}

  data_view & set_data_ptr(double * data_in) noexcept {_data = data_in; return *this;}
  data_view & set_size(std::size_t size_in) noexcept {_size = size_in; return *this;}
  data_view & set_extras(const T & extras_in) noexcept {_extras = extras_in; return *this;}
  /// compound assignment arithmetic operators
  ///@{
  data_view<T> & operator+=(const data_view<T> & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] += other_data[iy];
    return *this;
  }

  data_view<T> & operator-=(const data_view<T> & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] -= other_data[iy];
    return *this;
  }

  data_view<T> & operator*=(double val) {
    std::size_t sz = size();
    double * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] *= val;
    return *this;
  }
  data_view<T> & operator/=(double val) {
    std::size_t sz = size();
    double * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] /= val;
    return *this;
  }
  ///@}

  /// copy the data from other, assuming *this is initialised and of the correct size, etc.
  void copy_data(const data_view<T> & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    std::copy(other_data, other_data + sz, this_data);
  }

  const T & extras() const {return _extras;}
  T & extras() {return _extras;}

  void reset() {
    _data = nullptr;
    _size = 0;
    _extras = T();
  }

protected:
  double *      _data = nullptr;
  std::size_t   _size = 0;
  T _extras = T();

  inline std::tuple<std::size_t, double *, const double *> prepare_compound(const data_view<T> & b ) {
    extras().ensure_compatible(b.extras());
    return std::make_tuple(size(), data(), b.data());
  }
};


template<typename V, typename P>
class data_owner : public V {
public:
  data_owner() {}
  virtual ~data_owner() {del();}
  void del() {generic_delete(_ptr); reset();}

  virtual void alloc_virtual(const typename V::extras_type & extras) = 0;
  void reset() {
    V::reset();
    _ptr = nullptr;
  }

  data_owner & operator=(const data_owner & other) {
    if (_ptr == other._ptr) return *this; // self-assignment check
    copy(other);
    return *this;
  }

  data_owner(data_owner && other) noexcept  {move_no_del(other);}

  data_owner & operator=(data_owner && other) noexcept {
    if (this != &other) {
      del();
      move_no_del(other);
    }
    return *this;
  }

  /// @brief  make a copy of other, including allocating new storage if needed
  /// @param  other the other grid_quant to copy from
  /// @return a reference to the current object
  void copy(const V & other) {
    if (_ptr && this->extras() == other.extras()) {
      //std::cout << " reusing existing storage in grid_quant::copy\n";
      this->copy_data(other);
    } else {
      //std::cout << " allocating new storage in grid_quant::copy, size = " << other.size() << "\n";
      del();
      alloc_virtual(other.extras());
      this->copy_data(other);
      //std::cout << " data[50] = " << data()[50] << " v " << other.data()[50] << "\n";
    }
  }

  const P* ptr() const { return _ptr; }
  P* ptr()  { return _ptr; }

protected:

  /// @brief  move the contents from other into this, without deleting existing data
  /// @param  other: the other grid_quant to move from
  /// @return a reference to the current object
  ///
  /// Note that this does not delete any existing data in *this, so be careful to call del() 
  /// first if needed.
  void move_no_del(data_owner & other) noexcept {
    //std::cout << "actually moving " << other.ptr() << "\n"; 
    // move the semantics
    this->take_view(other);
    _ptr = other._ptr;
    other.reset();
  }

  P * _ptr = nullptr;
};

/// @brief Object-oriented wrapper around the grid_quant Fortran type, non-owning
///
/// This version provides a "view" onto an existing Fortran grid_quant object,
/// without taking ownership of it.
class grid_quant_view : public data_view<grid_def_view> {
public:

  grid_quant_view() {}
  grid_quant_view(double * data_ptr, std::size_t size, const grid_def_view & grid) 
    : data_view<grid_def_view>(data_ptr, size, grid) {}

  template<typename F, IsDoubleFunction<F> = 0>
  void assign(const F & fn) {
    if (!data()) {
      throw std::runtime_error("grid_quant_view::assign(fn): grid_quant_view object not associated");
    }
    // there is a design choice here: do we package whatever function
    // we have received into a Fortran-callable function, or do we
    // just do the filling on the C++ side? The latter is simpler
    // to implement, so let's do that for now, though it implies an extra allocation
    std::vector<double> yvals = grid().y_values();
    double * my_data = data();
    for (std::size_t iy=0; iy<size(); ++iy) {
      my_data[iy] = fn(yvals[iy]);
    }
  } 

  template<typename F, IsDoubleFunction<F> = 0>
  grid_quant_view & operator=(const F & fn) {assign(fn); return *this;}

  //template<> 
  grid_quant_view & operator=(const grid_quant_view & other) = default;// {copy_data(other); return *this;}

 
  //explicit grid_quant_view (const grid_quant_view & other) noexcept : data_view<grid_def_view>(other) {}   
//
//  /// @brief assignment operator, where the data from other is copied into this
//  /// @param other 
//  /// @return 
//  grid_quant_view & operator=(const grid_quant_view & other) {
//    this->copy_data(other);
//    return *this;
//  }
//
  template<typename T>
  double & operator[](T i) {return data()[i];}

  template<typename T>
  const double & operator[](T i) const {return data()[i];}

  /// return a ref to the associated grid definition
  const grid_def_view & grid() const {return _extras;}

  double at_y(double y) const {
    return hoppet_cxx__grid_quant__at_y(grid().ptr(), _data, y);
  }

  double at_x(double x) const {
    double y = std::log(1.0/x);
    return hoppet_cxx__grid_quant__at_y(grid().ptr(), _data, y);
  }

  /// @brief  compute the truncated moment up to grid's ymax
  /// @param n the moment index, where n=1 corresponds to momentum, n=0 to number
  /// @return the moment
  double truncated_moment(double n) const {
    return hoppet_cxx__grid_quant__trunc_mom(grid().ptr(), data(), n);
  }

  /// @brief  compute the truncated moment up to the specified y
  /// @param n the moment index, where n=1 corresponds to momentum, n=0 to number
  /// @param y the y value up to which to compute the moment
  /// @return the moment
  double truncated_moment(double n, double y) const {
    return hoppet_cxx__grid_quant__trunc_mom(grid().ptr(), data(), n, &y);
  }

};


//-----------------------------------------------------------------------------
class grid_quant : public data_owner<grid_quant_view, grid_quant_f> {
public:

  typedef grid_quant_view view_type;

  grid_quant() {}

  /// construct and allocate a grid_quant for the given grid
  grid_quant(const grid_def_view & grid) {alloc(grid);}

  // make sure we have the move constructor, move assignment and copy assignment
  // (they tend to get removed if other allocators are present, but we just want
  // them to do the default thing, since the move support etc is already in the
  // data_owner<...> base class.)
  grid_quant(grid_quant && other) noexcept = default;
  grid_quant & operator=(grid_quant && other) noexcept = default;
  grid_quant & operator=(const grid_quant & other) noexcept = default;

  void alloc(const grid_def_view & grid) {
    _extras = grid;
    _ptr    = hoppet_cxx__grid_quant__new(grid.ptr());
    _data   = hoppet_cxx__grid_quant__data_ptr(_ptr);
    _size   = static_cast<std::size_t>(grid.ny()+1);
  }
  void alloc_virtual(const grid_def_view & grid) override {alloc(grid);}
  

  /// construct and allocate a grid_quant for the given grid, and fill it
  /// with the specified function
  template<class F>
  grid_quant(const grid_def_view & grid, const F & fn) : grid_quant(grid) {
    assign(fn);
  }

  template<typename F>
  grid_quant & operator=(const F & fn) {assign(fn); return *this;}

  //grid_quant_view & view() {return *this;} // is this needed?

  /// copy constructor
  grid_quant(const grid_quant & other) {copy(other);}

  grid_quant(const grid_quant_view & other) {copy(other);}

  /// move constructor
  //grid_quant(grid_quant && other) noexcept {move_no_del(other);


//  /// @brief  assign from a function 
//  /// @tparam T generic function type
//  /// @param  fn any double(double) callable
//  /// @return a reference to the current object
//  template<typename F>
//  grid_quant & operator=(const F & fn) {
//    if (!_ptr) {
//      throw std::runtime_error("grid_quant::operator=(fn): grid_quant object not allocated");
//    }
//    // there is a design choice here: do we package whatever function
//    // we have received into a Fortran-callable function, or do we
//    // just do the filling on the C++ side? The latter is simpler
//    // to implement, so let's do that for now.
//    int ny = grid().ny();
//    std::vector<double> yvals = grid().y_values();
//    double * my_data = data();
//    for (std::size_t iy=0; iy<size(); ++iy) {
//      my_data[iy] = fn(yvals[iy]);
//    }
//    return *this;
//  }  
};

/// binary arithmetic operators
///@{

// these all use a copy-and-modify strategy, exploiting the move semantics
// for copy elision
inline grid_quant operator+(grid_quant a, const grid_quant_view & b) {a += b; return a;}
inline grid_quant operator-(grid_quant a, const grid_quant_view & b) {a -= b; return a;}
inline grid_quant operator*(grid_quant a, double b) {a *= b; return a;}
inline grid_quant operator*(double b, grid_quant a) {a *= b; return a;}
inline grid_quant operator/(grid_quant a, double b) {a /= b; return a;}

template<typename F> 
inline grid_quant operator*(const grid_def_view & grid, F && fn) {
  // create an output grid_quant
  grid_quant result(grid);
  result = fn;
  return result;
}

///@}

//-----------------------------------------------------------------------------
struct gq2d_extras {
  grid_def_view grid;
  std::size_t   stride;
  std::size_t   dim1_sz;
  gq2d_extras() : grid(), stride(0), dim1_sz(0) {}
  gq2d_extras(const gq2d_extras & other) : grid(other.grid), stride(other.stride), dim1_sz(other.dim1_sz) {}
  gq2d_extras(const grid_def_view & grid, std::size_t dim1_sz) : grid(grid), stride(grid.ny() + 1), dim1_sz(dim1_sz) {}
};
//-----------------------------------------------------------------------------
class grid_quant_2d_view : public data_view<gq2d_extras> {
public:  
  grid_quant_2d_view() {}
  grid_def_view grid() const {return extras().grid;}
  grid_quant_view operator[](std::size_t i) {
    grid_quant_view result(this->data() + i * extras().stride, extras().stride, extras().grid);
    return result;
  }
};


//-----------------------------------------------------------------------------
class grid_quant_2d : public data_owner<grid_quant_2d_view, grid_quant_2d_f> {

public:
  grid_quant_2d() {}
  grid_quant_2d(const grid_def_view & grid, std::size_t dim1_size) {
    alloc(gq2d_extras(grid, dim1_size));
  }
  void alloc_virtual(const gq2d_extras & extras_in) override {
    alloc(extras_in);
  }

  // make sure we have the move constructor, move assignment and copy assignment
  grid_quant_2d            (      grid_quant_2d && other) noexcept = default;
  grid_quant_2d & operator=(      grid_quant_2d && other) noexcept = default;
  grid_quant_2d & operator=(const grid_quant_2d &  other) noexcept = default;

  void alloc(gq2d_extras extras_in) {
    _extras = extras_in;
    _ptr    = hoppet_cxx__grid_quant_2d__new(extras_in.grid.ptr(), static_cast<int>(extras_in.dim1_sz));
    _data   = hoppet_cxx__grid_quant_2d__data_ptr(_ptr);
    _size   = _extras.stride * extras_in.dim1_sz;
  }
};


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_conv Fortran type, non-owning
class grid_conv_view {
public:
  grid_conv_view() {}
  grid_conv_view(const grid_def_view & grid) : _grid(grid) {}
  grid_conv_view(const grid_def_view & grid, grid_conv_f * ptr) : _grid(grid), _ptr(ptr) {}
  //grid_conv_view(grid_conv_f * ptr) : _ptr(ptr) {}
  const grid_def_view & grid() const {return _grid;}
  const grid_conv_f * ptr() const { return _ptr; }
  grid_conv_f * ptr() { return _ptr; }

  /// compound assignment arithmetic operators
  ///@{
  grid_conv_view & operator+=(const grid_conv_view & other) {
    grid().ensure_compatible(other.grid());
    hoppet_cxx_grid_conv__add(_ptr, other.ptr());
    return *this;
  }
  grid_conv_view & operator-=(const grid_conv_view & other) {
    grid().ensure_compatible(other.grid());
    double minus_one = -1.0;
    hoppet_cxx_grid_conv__add(_ptr, other.ptr(), &minus_one);
    return *this;
  }
  grid_conv_view & operator*=(double factor) {
    grid().ensure_valid();
    hoppet_cxx_grid_conv__multiply(_ptr, factor);
    return *this;
  }
  grid_conv_view & operator/=(double factor) {
    grid().ensure_valid();
    hoppet_cxx_grid_conv__multiply(_ptr, 1.0/factor);
    return *this;
  }
  ///@}

protected:
  grid_def_view  _grid;
  grid_conv_f *  _ptr = nullptr;
};

class grid_conv : public grid_conv_view {
public:
  grid_conv() {}

  /// construct a grid_conv object and initialise it with the given function
  ///
  /// @tparam Func  any callable type matching double(double y, int piece)
  /// @param grid   the grid definition
  /// @param conv_ignd_fn  the convolution integrand function
  ///
  template<typename Func>
  grid_conv(const grid_def_view & grid, Func && conv_ignd_fn) : grid_conv_view(grid) {

    // under the hood, the fortran calls hoppet_grid_conv_f__wrapper
    // (defined in hoppet_oo.cc), with a pointer to the function object
    std::function<double(double,int)> fn_ptr = std::forward<Func>(conv_ignd_fn);
    _ptr = hoppet_cxx_grid_conv__new_from_fn(grid.ptr(), &fn_ptr);
  }

  grid_conv(const grid_def_view & grid, grid_conv_f * ptr) : grid_conv_view(grid, ptr) {}

  grid_conv(const grid_conv & other) {
    _ptr  = hoppet_cxx_grid_conv__new_from_gc(other.ptr());
    _grid = other.grid();
    //std::cout << "copy constructing grid_conv, _ptr = " << _ptr << "\n";
  }

  grid_conv(grid_conv && other) noexcept
    : grid_conv_view(other) {
    other._ptr = nullptr;
  }


  void del() {if (_ptr) {hoppet_cxx_grid_conv__delete(&_ptr);}}
  ~grid_conv() {del();}

  grid_conv & operator=(grid_conv && other) noexcept {
    if (this != &other) {
      del();
      _ptr = other._ptr;
      other._ptr = nullptr;
    }
    return *this;
  }
};

inline grid_quant operator*(const grid_conv_view & conv, const grid_quant_view & q) {
  // create an output grid_quant
  if (conv.grid().ptr() != q.grid().ptr()) {
    throw std::runtime_error("grid_conv * grid_quant: grids do not match");
  }
  grid_quant result(q.grid());

  hoppet_cxx_grid_conv__times_grid_quant(conv.ptr(), q.data(), result.data());
  return result;
}

// these all use a copy-and-modify strategy, exploiting the move semantics
// for copy elision
inline grid_conv operator+(grid_conv a, const grid_conv_view & b) {a += b; return a;}
inline grid_conv operator-(grid_conv  a, const grid_conv_view & b) {a -= b; return a;}
inline grid_conv operator*(grid_conv a, double b) {a *= b; return a;}
inline grid_conv operator*(double b, grid_conv a) {a *= b; return a;}
inline grid_conv operator/(grid_conv a, double b) {a /= b; return a;}

inline grid_conv operator*(grid_conv_view const & a, grid_conv_view const & b) {
  a.grid().ensure_compatible(b.grid());
  grid_conv_f * res = hoppet_cxx_grid_conv__alloc_and_conv(a.ptr(), b.ptr());
  return grid_conv(a.grid(), res);
}

} // end namespace hoppet


/// objects globally defined in the streamlined interface
namespace hoppet {
namespace sl {
  /// a view of the grid_def object being used in the streamlined interface
  extern grid_def_view grid;
}
}

#endif // __HOPPET_OO__