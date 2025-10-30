#ifndef __HOPPET_OO__
#define __HOPPET_ OO__
#include "hoppet.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <tuple>

/// "forward" declaration of the Fortran grid_def type;
/// note that we only ever use pointers to it, so we do not
/// need to know anything about its actual structure, which
/// is hidden in the Fortran code.
class grid_def_f;
class grid_quant_f;


extern "C" {
  grid_def_f * hoppet_cxx__grid_def__new(double dy, double ymax, int order, double eps);
  void   hoppet_cxx__grid_def__delete(grid_def_f ** griddef);

  int    hoppet_cxx__grid_def__ny(grid_def_f * griddef);
  void   hoppet_cxx__grid_def__y_values(grid_def_f * griddef, double * yvals);
  void   hoppet_cxx__grid_def__x_values(grid_def_f * griddef, double * xvals);
}

extern "C" {
  void * hoppet_cxx__grid_quant__new(grid_def_f * griddef);
  void   hoppet_cxx__grid_quant__delete(void ** gridquant);
  double hoppet_cxx__grid_quant__at_y(void * gridquant, double y);
  double * hoppet_cxx__grid_quant__data_ptr(void * gridquant);

  //void   hoppet_cxx__grid_quant__set_zero(void * gridquant);
  //void   hoppet_cxx__grid_quant__copy_from(void * gridquant, void * other);
}

namespace hoppet {


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_def Fortran type, non-owning
///
/// This version provides a "view" onto an existing Fortran grid_def object,
/// without taking ownership of it.
class grid_def_view {
public:

  grid_def_view() {}

  grid_def_view(grid_def_f * ptr) : _ptr(ptr) {}
  grid_def_view(const grid_def_view & other) : _ptr(other._ptr) {}

  int ny() const {return hoppet_cxx__grid_def__ny(_ptr); }

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

  grid_def_f * ptr() const { return _ptr; }
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

  grid_def(double dy, double ymax, int order=-5, double eps=1e-7)
    : grid_def_view(hoppet_cxx__grid_def__new(dy, ymax, order, eps)) {}

  grid_def(grid_def_f * ptr) : grid_def_view(ptr) {}

  /// for now, disable copy construction to avoid double deletions
  grid_def(const grid_def & other) = delete;

  ~grid_def() {if (_ptr) hoppet_cxx__grid_def__delete(&_ptr); }

};  



//-----------------------------------------------------------------------------
class grid_quant {
public:
  /// construct and allocate a grid_quant for the given grid
  grid_quant(const grid_def_view & grid)
    : _ptr(hoppet_cxx__grid_quant__new(grid.ptr())), 
      _grid(grid), _owns_ptr(true) {}

  /// construct and allocate a grid_quant for the given grid, and fill it
  /// with the specified function
  template<class T>
  grid_quant(const grid_def_view & grid, const T & fn)
    : grid_quant(grid) {
    *this = fn;
  }

  /// copy constructor
  grid_quant(const grid_quant & other) {
    std::cout << "copy constructing\n";
    if (other._is_tmp && other._owns_ptr) move(other);
    else                     copy(other);
  }

  grid_quant(grid_quant && other) {
    std::cout << "move constructing\n";
    move(other);
  }

  /// @brief delete the underlying Fortran object if allocated and owned
  void del() {if (_ptr && _owns_ptr) hoppet_cxx__grid_quant__delete(&_ptr); _ptr=nullptr;}

  /// @brief destructor
  ~grid_quant() {del();}

  std::size_t size() const { return _grid.ny()+1; }

// move assignment
  grid_quant& operator=(grid_quant&& other) noexcept {
  std::cout << "move assigning\n";
  if (this != &other) {
    del();
    _grid = std::move(other._grid);
    _ptr = other._ptr;
    _owns_ptr = other._owns_ptr;
    _is_tmp = false;

    other._ptr = nullptr;
    other._owns_ptr = false;
    other._is_tmp = false;
  }
  return *this;
  }

  grid_quant & move(const grid_quant & other) {
    std::cout << "moving " << other.ptr() << "\n"; 
    del();
    // move the semantics
    _ptr = other._ptr;
    _grid = other._grid;
    _owns_ptr = other._owns_ptr;
    _is_tmp = false;
    other._ptr = nullptr;
    other._owns_ptr = false;
    return *this;
  }
  

  grid_quant & copy(const grid_quant & other) {
    std::cout << "copying " << other.ptr() << "\n";
    if (_ptr && grid().ptr() == other.grid().ptr()) return copy_data(other);
    else {
      del();
      _ptr = hoppet_cxx__grid_quant__new(other._grid.ptr());
      _grid = other.grid();
      _owns_ptr = true;
      _is_tmp = false;
      return copy_data(other);
    }
  }

  grid_quant & operator=(const grid_quant & other) {
    if (_ptr == other._ptr) return *this; // self-assignment check
    if (other._is_tmp)      return move(other);
    return copy(other);
  }

  /// @brief  assign from a function 
  /// @tparam T generic function type
  /// @param  fn any double(double) callable
  /// @return a reference to the current object
  template<class T>
  grid_quant & operator=(const T & fn) {
    // there is a design choice here: do we package whatever function
    // we have received into a Fortran-callable function, or do we
    // just do the filling on the C++ side? The latter is simpler
    // to implement, so let's do that for now.
    int ny = _grid.ny();
    std::vector<double> yvals = _grid.y_values();
    double * data = hoppet_cxx__grid_quant__data_ptr(_ptr); // hack to get data pointer
    for (int iy=0; iy<=ny; ++iy) {
      data[iy] = fn(yvals[iy]);
    }
    return *this;
  }


  const double * data() const {return hoppet_cxx__grid_quant__data_ptr(_ptr);}
  double       * data()       {return hoppet_cxx__grid_quant__data_ptr(_ptr);}

  template<typename T>
  double & operator[](T i) {return data()[i];}

  template<typename T>
  const double & operator[](T i) const {return data()[i];}

  /// return a ref to the associated grid definition
  const grid_def_view & grid() const {return _grid;}

  double at_y(double y) const {
    return hoppet_cxx__grid_quant__at_y(_ptr, y);
  }
  double at_x(double x) const {
    double y = std::log(1.0/x);
    return hoppet_cxx__grid_quant__at_y(_ptr, y);
  }

  grid_quant & mk_tmp() {
    _is_tmp = true;
    return *this;
  }
  bool is_tmp() const { return _is_tmp; }

  void * ptr() const { return _ptr; }
  
  /// unary arithmetic operators
  ///@{

  grid_quant & operator+=(const grid_quant & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] += other_data[iy];
    return *this;
  }

  grid_quant & operator-=(const grid_quant & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] -= other_data[iy];
    return *this;
  }

  grid_quant & operator*=(double val) {
    std::size_t sz = size();
    double * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] *= val;
    return *this;
  }
  grid_quant & operator/=(double val) {
    std::size_t sz = size();
    double * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] /= val;
    return *this;
  }

  ///@}

  /// binary arithmetic operators
  ///@{
  grid_quant operator+(const grid_quant & other) const {
    grid_quant new_gq(_grid);    
    new_gq.mk_tmp();
    double * new_gq_data = new_gq.data();
    const double * this_data = data();
    const double * other_data = other.data();
    std::size_t sz = size();
    //auto [sz, new_grid, new_grid_data, this_data, other_data] = prepare_binary(other);
    for (std::size_t iy=0; iy<sz; ++iy) new_gq_data[iy] = this_data[iy] + other_data[iy];
    return new_gq;
  }
  ///@}



protected:  

  inline std::tuple<std::size_t, double *, const double *> prepare_compound(const grid_quant & other) {
    assert( _grid.ptr() == other._grid.ptr() && "grid_quant unary op with another grid_quant: grids do not match");
    double * this_data = data();
    const double * other_data = other.data();
    return std::make_tuple(size(), this_data, other_data);
  }

  //inline std::tuple<std::size_t, grid_quant, double *, const double *, const double *> prepare_binary(const grid_quant & other) const {
  //  assert( _grid.ptr() == other._grid.ptr() && "grid_quant binary op with another grid_quant: grids do not match");
  //  grid_quant new_gq(_grid);
  //  new_gq._is_tmp = true;
  //  double * new_gq_data = new_gq.data();
  //  const double * this_data = data();
  //  const double * other_data = other.data();
  //  return std::make_tuple(size(), new_gq, new_gq_data, this_data, other_data);
  //}


  /// copy the data from other, assuming *this is initialised
  grid_quant & copy_data(const grid_quant & other) {
    int ny = _grid.ny();
    const double * src_data = other.data();
    double * dst_data = data();
    for (int iy=0; iy<=ny; ++iy) {
      dst_data[iy] = src_data[iy];
    }
    return *this;
  }


  grid_def_view _grid;
  mutable void * _ptr = nullptr;
  mutable bool _owns_ptr = false;
  mutable bool _is_tmp = false;
};

} // end namespace hoppet
#endif // __HOPPET_OO__