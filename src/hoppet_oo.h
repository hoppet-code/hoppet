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
  grid_def_f * hoppet_cxx__grid_def__new_from_grids(grid_def_f ** griddefs, int ngrids, bool locked);
  grid_def_f * hoppet_cxx__grid_def__new_default(double dy, double ymax, int order = -6);
  grid_def_f * hoppet_cxx__grid_def__copy(grid_def_f * griddef);
  void   hoppet_cxx__grid_def__delete(grid_def_f ** griddef);

  int    hoppet_cxx__grid_def__ny(grid_def_f * griddef);
  void   hoppet_cxx__grid_def__y_values(grid_def_f * griddef, double * yvals);
  void   hoppet_cxx__grid_def__x_values(grid_def_f * griddef, double * xvals);
}

extern "C" {
  grid_quant_f * hoppet_cxx__grid_quant__new(const grid_def_f * griddef);
  void   hoppet_cxx__grid_quant__delete(grid_quant_f ** gridquant);
  double hoppet_cxx__grid_quant__at_y(const grid_quant_f * gridquant, double y);
  double * hoppet_cxx__grid_quant__data_ptr(grid_quant_f * gridquant);

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
    std::vector<grid_def_f *> grid_ptrs(ngrids);
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

class grid_quant; // forward declaration

class grid_quant_view {
public:

  grid_quant_view() {}
  grid_quant_view(const grid_def_view & grid, grid_quant_f * ptr)
    : _grid(grid), _ptr(ptr) {}

  grid_quant_view operator=(const grid_quant & other) = delete; // prevent assignment
  //grid_quant_view(const grid_quant & other) = delete; // prevent copying from owning class

  //grid_quant_view view(const grid_quant_view & other) {} 

  std::size_t size() const { return _grid.ny()+1; }

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

  void * ptr() const { return _ptr; }
  
  /// unary arithmetic operators
  ///@{

  grid_quant_view & operator+=(const grid_quant_view & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] += other_data[iy];
    return *this;
  }

  grid_quant_view & operator-=(const grid_quant_view & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] -= other_data[iy];
    return *this;
  }

  grid_quant_view & operator*=(double val) {
    std::size_t sz = size();
    double * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] *= val;
    return *this;
  }
  grid_quant_view & operator/=(double val) {
    std::size_t sz = size();
    double * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] /= val;
    return *this;
  }

  ///@}


protected:
  /// @brief pointer to the underlying Fortran grid_quant object
  grid_def_view _grid;
  grid_quant_f * _ptr = nullptr;

  inline std::tuple<std::size_t, double *, const double *> prepare_compound(const grid_quant_view & other) {
    assert( _grid.ptr() == other._grid.ptr() && "grid_quant unary op with another grid_quant: grids do not match");
    double * this_data = data();
    const double * other_data = other.data();
    return std::make_tuple(size(), this_data, other_data);
  }


};

//-----------------------------------------------------------------------------
class grid_quant : public grid_quant_view {
public:

  grid_quant() {}

  /// construct and allocate a grid_quant for the given grid
  grid_quant(const grid_def_view & grid) : 
      grid_quant_view(grid, hoppet_cxx__grid_quant__new(grid.ptr())) {}

  /// construct and allocate a grid_quant for the given grid, and fill it
  /// with the specified function
  template<class T>
  grid_quant(const grid_def_view & grid, const T & fn) : grid_quant(grid) {
    *this = fn;
  }

  //grid_quant_view & view() {return *this;} // is this needed?

  /// copy constructor
  grid_quant(const grid_quant & other) {
    //std::cout << "copy constructing\n";
    copy(other);
  }

  grid_quant(const grid_quant_view & other) {
    //std::cout << "copy constructing\n";
    copy(other);
  }

  /// move constructor
  grid_quant(grid_quant && other) noexcept {
    move_no_del(other);
  }

  /// @brief delete the underlying Fortran object if allocated
  void del() {if (_ptr) hoppet_cxx__grid_quant__delete(&_ptr); _ptr=nullptr;}

  /// @brief destructor
  ~grid_quant() {del();}

  // move assignment
  grid_quant& operator=(grid_quant&& other) noexcept {
    //std::cout << "move assigning\n";
    if (this != &other) {
      del();
      move_no_del(other);
    }
    return *this;
  }


  grid_quant & copy(const grid_quant_view & other) {
    if (_ptr && grid().ny() == other.grid().ny()) {
      return copy_data(other);
    }
    else {
      del();
      _ptr = hoppet_cxx__grid_quant__new(other.grid().ptr());
      _grid = other.grid();
      return copy_data(other);
    }
  }

  grid_quant & operator=(const grid_quant & other) {
    if (_ptr == other._ptr) return *this; // self-assignment check
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




protected:  

  /// copy the data from other, assuming *this is initialised and of the correct size, etc.
  grid_quant & copy_data(const grid_quant_view & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    std::copy(other_data, other_data + sz, this_data);
    return *this;
  }

  /// @brief  move the contents from other into this, without deleting existing data
  /// @param  other the other grid_quant to move from
  /// @return a reference to the current object
  ///
  /// Note that this does not delete any existing data in *this, so be careful to call del() 
  /// first if needed.
  grid_quant & move_no_del(grid_quant & other) noexcept {
    //std::cout << "actually moving " << other.ptr() << "\n"; 
    // move the semantics
    _ptr = other._ptr;
    _grid = other._grid;
    other._ptr = nullptr;
    return *this;
  }
  

  //xgrid_def_view _grid;
  //xgrid_quant_f * _ptr = nullptr;
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

///@}


} // end namespace hoppet


/// objects globally defined in the streamlined interface
extern grid_def_f * hoppet_sl_grid_ptr;
namespace hoppet {
namespace sl {
  /// a view of the grid_def object being used in the streamlined interface
  extern grid_def_view grid;
}
}

#endif // __HOPPET_OO__