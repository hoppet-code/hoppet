#ifndef __HOPPET_OO__
#define __HOPPET_OO__
#include "hoppet.h"
#include <vector>

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
/// @brief Object-oriented wrapper around the grid_def Fortran type
///
/// @todo allow for nested grids
/// @todo 
class grid_def {
public:


  grid_def() {};

  grid_def(double dy, double ymax, int order=-5, double eps=1e-7)
    : _ptr(hoppet_cxx__grid_def__new(dy, ymax, order, eps)) {}

  grid_def(grid_def_f * ptr, bool owns_ptr=false)
    : _ptr(ptr), _owns_ptr(owns_ptr) {}

  ~grid_def() {if (_ptr && _owns_ptr) hoppet_cxx__grid_def__delete(&_ptr); }

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
private:
  grid_def_f * _ptr = nullptr;
  bool _owns_ptr = true;
};  

//-----------------------------------------------------------------------------
class grid_quant {
public:
  /// construct and allocate a grid_quant for the given grid
  grid_quant(const grid_def & grid)
    : _ptr(hoppet_cxx__grid_quant__new(grid.ptr())), 
      _grid(grid.ptr(),false), _owns_ptr(true) {}

  /// construct and allocate a grid_quant for the given grid, and fill it
  /// with the specified function
  template<class T>
  grid_quant(const grid_def & grid, const T & fn)
    : grid_quant(grid) {
    *this = fn;
  }


  grid_quant(const grid_quant & other) {copy(other);}

  void del() {if (_ptr && _owns_ptr) hoppet_cxx__grid_quant__delete(&_ptr); _ptr=nullptr;}
  ~grid_quant() {del();}

  grid_quant & move(const grid_quant & other) {
    del();
    // move the semantics
    _ptr = other._ptr;
    _grid = other._grid;
    _owns_ptr = other._owns_ptr;
    other._ptr = nullptr;
    other._owns_ptr = false;
    return *this;
  }

  // copy the data, assuming this is initialised
  grid_quant & copy_data(const grid_quant & other) {
    int ny = _grid.ny();
    const double * src_data = other.data();
    double * dst_data = data();
    for (int iy=0; iy<=ny; ++iy) {
      dst_data[iy] = src_data[iy];
    }
    return *this;
  }

  grid_quant & copy(const grid_quant & other) {
    if (_ptr && grid().ptr() == other.grid().ptr()) return copy_data(other);
    else {
      del();
      _ptr = hoppet_cxx__grid_quant__new(other._grid.ptr());
      _grid = grid_def(other._grid.ptr(), false);
      _owns_ptr = true;
      _is_tmp = false;
      return copy_data(other);
    }
  }

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

  grid_quant & operator=(const grid_quant & other) {
    if (_ptr == other._ptr) return *this; // self-assignment check
    if (other._is_tmp)      return move(other);
    return copy(other);
  }

  const double * data() const {
    return hoppet_cxx__grid_quant__data_ptr(_ptr);
  }
  double * data() {
    return hoppet_cxx__grid_quant__data_ptr(_ptr);
  }

  template<typename T>
  double & operator[](T i) {return data()[i];}

  template<typename T>
  const double & operator[](T i) const {return data()[i];}
    
  /// return a ref to the associated grid definition
  const grid_def & grid() const {return _grid;}

  double at_y(double y) const {
    return hoppet_cxx__grid_quant__at_y(_ptr, y);
  }

  void * ptr() const { return _ptr; }


private:
  grid_def _grid;
  mutable void * _ptr = nullptr;
  mutable bool _owns_ptr = false;
  mutable bool _is_tmp = false;
};

} // end namespace hoppet
#endif // __HOPPET_OO__