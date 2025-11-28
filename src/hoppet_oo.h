#ifndef __HOPPET_OO__
#define __HOPPET_OO__
#include <vector>
#include <cmath>
#include <cassert>
#include <tuple>
#include <functional>
#include <concepts>
#include <optional>
#include "hoppet.h"
#include "hoppet/base_types.h"
// see also hoppet_cxx_oo.f90 for the Fortran prototypes
#include "hoppet/fortran_prototypes.h"

// Elements to think about
// - [~] do we separate things out into different files (done to some extent)
// - [~] do we provide "physics" aliases: e.g. grid_quant -> pdf_flav, grid_conv -> split_fn
// - [x] make sure hoppet/qcd.h (etc) gets installed properly
// - [ ] think about hoppet__qcd v hoppet_cxx__qcd naming conventions
// - [x] making sure we treat default epsilon correctly
// - [ ] where does hoppet_oo.h go?
// - [x] check everything gets installed
// - [ ] nomenclature: pdf_flav, v. pdf v. pdf_set (old hoppet, set=allflav; LHAPDF set: = error set)
// - [ ] downgrade from concepts to C++17 is_invocable_r_v (e.g. for swig compatibility)?

// Next steps:
// - [x] add the running_coupling class
// - [x] add the dglap_holder    
// - [x] generic string conversion (hoppet_c_f_string_utils.f90)
// - [~] add the pdf_table
//       - [ ] grid_quant_3d?
//       - [x] evolution functions
//       - [x] PDF access functions
//       - [ ] does copying also copy the evolution operators? **NO**
//       - [ ] applying an evolution operator from a different dglap_holder?
//       - [ ] moving a pdf_table
//       - [x] fill from LHAPDF type functions
//       - [x] write to LHAPDF
//       - [x] access to lambda_eff
// - [~] array of tables
//       - [x] view creation and indexing
//       - [x] array creation
//       - [ ] access to efficient at_xQ functions
// - [ ] add the evln_operator class
// - [ ] the evolution function
// - [x] add split option when creating splitting functions
// - [ ] functions to query & change flavour representations
// - [x] add the mass thresholds 
// - [ ] grid_quant_3d?
// - [~] add streamlined interface functions
// - [ ] support for probes
// - [ ] add qed support, including pdf_table allocators
// - [x] add access to things like the beta function coefficients, qcd group constants, etc.
// - [ ] add structure function support
// - [ ] test installed version (re include paths)
// - [~] add documentation

// Things to perhaps add
// - [x] add take_view to the obj_view class
// - [x] basic checks of grid compatibility, so as to get C++ errors rather than Fortran errors
// - [x] mvv-like interface for splitting functions?



#define RETURN_INT_MEMBER(classname, membername)           inline int    membername() const {return hoppet_cxx__##classname##__##membername(valid_ptr());} 
#define RETURN_LOG_MEMBER(classname, membername)           inline bool   membername() const {return hoppet_cxx__##classname##__##membername(valid_ptr());} 
#define RETURN_DBL_MEMBER(classname, membername)           inline double membername() const {return hoppet_cxx__##classname##__##membername(valid_ptr());} 
#define RETURN_INT_MEMBER_I(classname, membername)         inline int    membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_DBL_MEMBER_I(classname, membername)         inline double membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_OBJ_MEMBER(classname, membername, typename) inline typename##_view membername() const {return typename##_view(hoppet_cxx__##classname##__##membername(valid_ptr()));} 
#define RETURN_OBJ_MEMBER_I( classname, membername, typename) inline typename##_view membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_OBJ_MEMBER_IJ(classname, membername, typename) inline typename##_view membername(int i, int j) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i, j);} 

#define RETURN_OBJ_MEMBER_REF(classname, membername, typename) inline typename##_ref membername() const {return typename##_ref(hoppet_cxx__##classname##__##membername(valid_ptr()));} 
#define RETURN_OBJ_MEMBER_REF_I( classname, membername, typename) inline typename##_ref membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_OBJ_MEMBER_REF_IJ(classname, membername, typename) inline typename##_ref membername(int i, int j) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i, j);} 

/// Namespace hoppet contains the object-oriented C++ interface to HOPPET,
/// as well as various integer constants that are useful in the streamlined
/// interface (e.g. hoppet::iflv_g for the gluon flavour index).
///
namespace hoppet {

typedef std::size_t size_type;

//template <typename F>
//using IsDoubleFunction = std::enable_if_t<
//    std::is_invocable_r_v<double, F, double>, int
//>;
//template <typename F>
//concept DoubleFunction =
//    std::is_invocable<F, double> &&
//    std::same_as<std::invoke_result_t<F, double>, double>;
//template <typename F>
//concept DoubleFnDouble =
//    std::invocable<F, double> &&
//    std::same_as<std::invoke_result_t<F, double>, double>;
//
//template <typename F>
//concept DoubleFnDoubleInt =
//    std::invocable<F, double, int> &&
//    std::same_as<std::invoke_result_t<F, double, int>, double>;
//template <typename F>
//concept VoidFnDoubleDoubleDoubleptr =
//    std::invocable<F, double, double, double *> &&
//    std::same_as<std::invoke_result_t<F, double, double, double *>, void >;
//
//template <typename F>
//concept VoidFnDoubleDoubleDoublevec =
//    std::invocable<F, double, double, std::vector<double> &> &&
//    std::same_as<std::invoke_result_t<F, double, double, std::vector<double> &>, void >;



class grid_def_view; // forward declaration

//-----------------------------------------------------------------------------
/// @brief A reference to an underlying Fortran object that holds a grid definition
///
class grid_def_ref : public obj_view<grid_def_f> {
public:

  using base_t = obj_view<grid_def_f>;
  using base_t::base_t; // ensures that constructors are inherited
  using view_t = grid_def_view;
  using ref_t = grid_def_ref;

  /// return the upper limit of the iy index (i.e. highest valid index)
  int ny() const {return hoppet_cxx__grid_def__ny(valid_ptr()); }

  /// return the size of the grid (i.e. total number of points)
  std::size_t size() const {return static_cast<std::size_t>(ny())+1;}

  /// @brief return a vector of the y=ln(1/x) values for the grid points 
  ///
  /// The returned vector has size size()=ny()+1, and for nested grids
  /// some y points are likely to be repeated.
  std::vector<double> y_values() const {
    std::vector<double> yvals(size());
    hoppet_cxx__grid_def__y_values(valid_ptr(), yvals.data());
    return yvals;
  }
  /// @brief return a vector of the x values for the grid points 
  std::vector<double> x_values() const {
    std::vector<double> xvals(size());
    hoppet_cxx__grid_def__x_values(valid_ptr(), xvals.data());
    return xvals;
  }

  /// @brief return true if both grids non-null and `other` is equivalent to this
  inline bool operator==(const grid_def_ref & other) const noexcept {
    // equivalence is only true if both pointers are non-null and
    // they are either the same pointer or equivalent grids
    return ptr() && other.ptr() && (
      ptr() == other.ptr() || hoppet_cxx__grid_def__equiv(ptr(), other.ptr())
    );
  }
  /// @brief return true if either grid is null or `other` is not equivalent to this
  inline bool operator!=(const grid_def_ref & other) const noexcept {
    return !(*this == other);
  }

  /// throw a runtime_error if the other grid is not equivalent to this
  void ensure_compatible(const grid_def_ref & other) const {
    if (*this != other) {
      throw std::runtime_error(
        "hoppet::grid_def_view::ensure_compatible: grids are not equivalent");
    }
  }

  
  RETURN_DBL_MEMBER(grid_def,dy)     ///< return the grid's dy spacing (coarsest if using subgrids)
  RETURN_DBL_MEMBER(grid_def,ymax)   ///< return the grid's maximum y value (largest among subgrids)
  RETURN_DBL_MEMBER(grid_def,eps)    ///< return the grid's  adaptive integration precision (only for individual grids)  
  RETURN_INT_MEMBER(grid_def,nsub)   ///< return the number of subgrids (0 if this in an individual grid)
  RETURN_INT_MEMBER(grid_def,order)  ///< return the interpolation order (meaninful only for individual grids)
  RETURN_LOG_MEMBER(grid_def,locked) ///< return true if the subgrids are locked (only valid with subgrids)

  /// @brief return the grid_def_view for subgrid i (index runs from 1 to nsub()).
  /// 
  /// If the subgrids are locked, subgrids are ordered in increasing dy and ymax,
  /// otherwise the order is that in which they were provided at construction.
  RETURN_OBJ_MEMBER_REF_I(grid_def,subgd,grid_def) 

  /// @brief return the iy index at which subgrid i starts
  RETURN_INT_MEMBER_I(grid_def,subiy) 

};

//-----------------------------------------------------------------------------
/// @brief A view of a grid_def object, i.e. a definition of a grid in y=ln(1/x)
class grid_def_view : public grid_def_ref {
public:
  using base_t = grid_def_ref;
  using base_t::base_t; // ensures that constructors are inherited
  using view_t = grid_def_view;
  using ref_t = grid_def_ref;

};

//-----------------------------------------------------------------------------
/// @brief A definition of a grid in \f$y=ln(1/x)\f$.
///
/// Most users will want to construct the grid_def object with the function
/// hoppet::grid_def_default, e.g. 
///
/// ```cpp
/// double dy   = 0.1;  /// grid spacing, sufficient for about $10^{}
/// double ymax = 12.0;
/// int order   = -6;
/// hoppet::grid_def grid = hoppet::grid_def_default(dy, ymax, order);
/// ```
/// which directly constructs a nested grid that provides good
/// accuracy over a wide range of x values.
///
/// See grid_def_view for additional member functions.
///
class grid_def : public obj_owner<grid_def_view> {
public:

  using base_t = obj_owner<grid_def_view>;
  using base_t::base_t; // ensures that constructors are inherited
  using view_t = grid_def_view;
  using ref_t  = grid_def_ref;

  /// @brief construct and allocate a new grid_def object
  ///
  /// @param dy      grid spacing
  /// @param ymax    maximum y value
  /// @param order   usual interpolation order parameter
  /// @param eps     accuracy for the adaptive integration
  //
  explicit grid_def(double dy, double ymax, int order=-5, double eps=hoppet_default_conv_eps)
    : grid_def(hoppet_cxx__grid_def__new(dy, ymax, order, eps)) {}

  /// @brief construct a grid_def from multiple grid_def objects
  ///
  /// @param grids   vector of grid_def objects
  /// @param locked  whether the new grid should be "locked"
  explicit grid_def(const std::vector<grid_def_view> & grids, bool locked=false) {
    int ngrids = static_cast<int>(grids.size());
    std::vector<const grid_def_f *> grid_ptrs(ngrids);
    for (int i=0; i<ngrids; ++i) {
      grid_ptrs[i] = grids[i].ptr();
    }
    _ptr = hoppet_cxx__grid_def__new_from_grids(grid_ptrs.data(), ngrids, locked);
  }

  /// the copy constructor and copy assignment operator are disabled to avoid
  /// accidental copies of grid_def objects; while the copies themselves
  /// would work correctly, if you construct a grid_quant (or similar) from
  /// a copied grid, if that grid goes out of scope, the grid_quant would
  /// be left with a dangling pointer. Eliminating the copy constructor
  /// encourages users to keep a single grid_def object around and -- if needed
  /// -- to pass grid_def_view to pass references to other objects.
  grid_def(const grid_def &) = delete; // disable copy constructor
  grid_def & operator=(const grid_def &) = delete; // disable copy assignment operator

  grid_def & operator=(grid_def && other) noexcept {
    if (this != &other) {
      _ptr = other._ptr;
      other._ptr = nullptr;
    }
    return *this;
  }
  grid_def(grid_def && other) noexcept
    : obj_owner<grid_def_view>(other._ptr) {
    other._ptr = nullptr;
  }

};  

/// @brief return a grid_def object with HOPPET's default choice of nested, locked grids
///
/// @param dy      spacing of the coarsest grid
/// @param ymax    maximum y value for the coarsest grid
/// @param order   usual interpolation order parameter
/// @return the grid_def object
inline grid_def grid_def_default(double dy, double ymax, int order) {
  return grid_def(hoppet_cxx__grid_def__new_default(dy, ymax, order));
}

class grid_quant; // forward declaration

//-----------------------------------------------------------------------------
/// @brief A view of a grid_quant object, i.e. a function of y=ln(1/x) stored at fixed grid points in y
///
/// This does not own the underlying data, but can still manipulate it.
///
class grid_quant_view : public data_view<grid_def_view> {
public:

  grid_quant_view() {}

  /// @brief  construct a grid_quant_view from existing data and grid (internal use only)
  ///
  /// @param data_ptr pointer to the data array
  /// @param size     size of the data array
  /// @param grid     associated grid definition
  grid_quant_view(double * data_ptr, std::size_t size, const grid_def_view & grid) 
    : data_view<grid_def_view>(data_ptr, size, grid) {}


  /// @brief  assuming the grid_quant_view's storage has been set up, assign from a function 
  /// @tparam T generic function type
  /// @param  fn any double(double) callable
  ///
  template<typename F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double>>>
  void assign(const F & fn) {
  //void assign(const DoubleFnDouble auto & fn) {
    ensure_valid();
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

  /// @brief assign a constant value to all elements
  /// @param val   the value to assign
  void assign(double val) {
    ensure_valid(); 
    std::fill(data(), data()+size(), val);
  }

  /// @brief  ensure that the grid_quant_view is valid (i.e. associated with data), otherwise throw an exception
  void ensure_valid() const {
    if (!data()) {throw std::runtime_error("hoppet::grid_quant_view::ensure_valid: data pointer is null");}
  }

  /// @brief  assign from a function, similar to assign()
  template<typename F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double>>>
  grid_quant_view & operator=(const F & fn) {assign(fn); return *this;}

  /// @brief  copy assignment operator, which copies data from other, assuming
  ///         that this is already allocated and of the correct size, etc.
  ///
  /// @param other the other grid_quant_view to copy from
  //template<> 
  grid_quant_view & operator=(const grid_quant_view & other) = default;// {copy_data(other); return *this;}

  // /// @brief  assign a constant value to all elements
  // /// @param val the value to assign
  grid_quant_view & operator=(double val) {assign(val); return *this;}

  /// @brief  indexing operator
  /// @param i index of the element to access
  /// @return reference to the element at index i
  ///
  /// Note that y =log(1/x) do _not_ monotonically increase with i, because
  /// of hoppet's nested grid structure.
  double & operator[](std::size_t i) {return data()[i];}

  const double & operator[](std::size_t i) const {return data()[i];}

  /// return a ref to the associated grid definition
  const grid_def_view & grid() const {return _extras;}

  /// @brief    return the interpolated value at the specified y=ln(1/x)
  /// @param y  the y value to interpolate
  /// @return   the interpolated value
  double at_y(double y) const {
    return hoppet_cxx__grid_quant__at_y(grid().ptr(), _data, y);
  }

  /// @brief    return the interpolated value at the specified x
  /// @param x  the x value to interpolate
  /// @return   the interpolated value
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

  /// @brief  compute the luminosity with another grid_quant_view
  /// @param other   the other grid_quant_view
  /// @return the luminosity grid_quant (in a dimensionless definition)
  grid_quant luminosity(const grid_quant_view & other) const;
  // {
  //  grid().ensure_compatible(other.grid());
  //  grid_quant result(grid());
  //  hoppet_cxx__grid_quant__luminosity(
  //    grid().ptr(),
  //    data(),
  //    other.data(),
  //    result.data()
  //  );
  //  return result;
  //}

};


//-----------------------------------------------------------------------------
/// @brief A function of y=ln(1/x) stored at fixed grid points in y (see grid_quant_view for core member functions)
///
/// Note that Fortran hoppet does not use grid_quant objects by default, instead it just 
/// uses `real(dp) :: gq(0:ny)` arrays. This C++ wrapper instead explicitly
/// goes via a fortran grid_quant object, which manages the data array internally.
class grid_quant : public data_owner<grid_quant_view, grid_quant_f> {
public:

  using view_t = grid_quant_view;

  grid_quant() {}

  /// construct and allocate a grid_quant for the given grid
  grid_quant(const grid_def_view & grid) {alloc(grid);}

  // make sure we have the move constructor, move assignment and copy assignment
  // (they tend to get removed if other allocators are present, but we just want
  // them to do the default thing, since the move support etc is already in the
  // data_owner<...> base class.)
  grid_quant(const grid_quant & other)      {copy(other);}
  grid_quant(const grid_quant_view & other) {copy(other);}
  grid_quant(grid_quant && other) noexcept = default;
  grid_quant & operator=(grid_quant && other) noexcept = default;
  grid_quant & operator=(const grid_quant & other) noexcept = default;

  void alloc(const grid_def_view & grid) {
    _extras = grid;
    _ptr    = hoppet_cxx__grid_quant__new(grid.ptr());
    _data   = hoppet_cxx__grid_quant__data_ptr(_ptr);
    _size   = static_cast<std::size_t>(grid.size());
  }
  void alloc_virtual(std::size_t /*sz*/, const grid_def_view & grid) override {alloc(grid);}
  

  /// construct and allocate a grid_quant for the given grid, and fill it
  /// with the specified function
  template<class F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double>>>
  grid_quant(const grid_def_view & grid, const F & fn) : grid_quant(grid) {
    assign(fn);
  }

  /// assign from a function, similar to assign()
  template<class F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double>>>
  grid_quant & operator=(const F & fn) {assign(fn); return *this;}

  grid_quant & operator=(double val) {assign(val); return *this;}

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
inline grid_quant operator-(const grid_quant_view & a) {return -1.0 * a;}

//template<typename F> 
template<typename F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double>>>
inline grid_quant operator*(const grid_def_view & grid, const F & fn) {
  // create an output grid_quant
  grid_quant result(grid);
  result = fn;
  return result;
}

///@}

// actual implementation of luminosity, needs to be defined out of the grid_quant_view
// class, because it returns a grid_quant (not yet defined there)
inline grid_quant grid_quant_view::luminosity(const grid_quant_view & other) const {
  grid().ensure_compatible(other.grid());
  grid_quant result(grid());
  hoppet_cxx__grid_quant__luminosity(
    grid().ptr(),
    data(),
    other.data(),
    result.data()
  );
  return result;
}


/// @brief  an alias for grid_quant_view, i.e. representing a single PDF flavour across y=ln(1/x)
typedef grid_quant_view pdf_flav_view;
/// @brief  an alias for grid_quant, i.e. representing a single PDF flavour across y=ln(1/x)
typedef grid_quant      pdf_flav;


//-----------------------------------------------------------------------------
/// @brief Internal class related to storage of properties of grid_quant_2d objects
struct gq2d_extras {
  grid_def_view grid;
  std::size_t   size_dim1;
  std::size_t   size_dim0;
  gq2d_extras() : grid(), size_dim1(0), size_dim0(0) {}
  gq2d_extras(const gq2d_extras & other) : grid(other.grid), size_dim1(other.size_dim1), size_dim0(other.size_dim0) {}
  gq2d_extras(const grid_def_view & grid_in, std::size_t size_dim0_in) : grid(grid_in), size_dim1(grid_in.size()), size_dim0(size_dim0_in) {}
  void ensure_compatible(const gq2d_extras & other) const {
    if (size_dim0 != other.size_dim0) throw std::runtime_error("hoppet::gq2d_extras::ensure_compatible: incompatible grid_quant_2d dim1_sz");
    grid.ensure_compatible(other.grid);
  }
  bool operator==(const gq2d_extras & other) const {
    return size_dim0 == other.size_dim0 && grid == other.grid;
  }
  bool operator!=(const gq2d_extras & other) const { return !(*this == other); }
};

//-----------------------------------------------------------------------------
/// @brief A non-owning view of a grid_quant_2d object, i.e. a function of
///        y=ln(1/x) stored at fixed grid points in y for each of several flavour
///        indices
class grid_quant_2d_view : public data_view<gq2d_extras> {
public:  
  using base_t = data_view<gq2d_extras>;
  using base_t::base_t; // ensures that constructors are inherited

  grid_quant_2d_view() {}
  grid_def_view grid() const {return extras().grid;}
  grid_quant_view operator[](std::size_t i) {
    grid_quant_view result(this->data() + i * extras().size_dim1, extras().size_dim1, extras().grid);
    return result;
  }
  grid_quant_view operator()(std::size_t i) {return (*this)[i];}
  double & operator()(std::size_t iflv, std::size_t iy) {return this->data()[iflv * extras().size_dim1 + iy];}
  const double & operator()(std::size_t iflv, std::size_t iy) const {return this->data()[iflv * extras().size_dim1 + iy];}  


  /// assign using a function f(x,Q, xpdf_array), similar to the LHAPDF
  /// "evolve" fortran interface, where xpdf_array is a double* pointing 
  /// to an array with at least 13 entries
  template<typename F, std::enable_if_t<std::is_invocable_r_v<void, F, double, double, double*>, int> = 0>
  void assign_xQ_into(const F & fn, double Q) {
    ensure_valid();
    if (size_flv() < iflv_max+1) {
      throw std::runtime_error("grid_quant_2d_view::assign(fn): grid_quant_2d_view size_flv() = " 
                  + std::to_string(size_flv()) + " < iflv_max+1 = " + std::to_string(iflv_max+1));
    }
    std::vector<double> xvals = grid().x_values();
    std::vector<double> xpdf (size_flv());
    for (std::size_t iy=0; iy < grid().size(); ++iy) {
      fn(xvals[iy], Q, xpdf.data());
      for (std::size_t iflv=0; iflv < size_flv(); ++iflv) {
        (*this)(iflv,iy) = (iflv <= iflv_max) ? xpdf[iflv] : 0.0;
      }
    }
  }

  /// assign using a function f(x,Q, vector<double> & pdf_vec), similar to the LHAPDF
  /// C++ `PDF::xfxQ(x,Q,pdf_vec)` member function.
  template<typename F, std::enable_if_t<std::is_invocable_r_v<void, F, double, double, std::vector<double> &>, int> = 0>
  void assign_xQ_into(const F & fn, double Q) {
    ensure_valid();
    std::vector<double> xvals = grid().x_values();
    std::vector<double> xpdf (size_flv());
    for (std::size_t iy=0; iy < grid().size(); ++iy) {
      fn(xvals[iy], Q, xpdf);
      // check only the first time around
      if (iy == 0 && xpdf.size() > size_flv()) {
        throw std::runtime_error("grid_quant_2d_view::assign(fn): supplied pdf_vec size = " 
                  + std::to_string(xpdf.size()) + " > grid_quant_2d_view size_flv() = " + std::to_string(size_flv()));
      }
      // up to the flavours that have been supplied, copy the values
      for (std::size_t i=0; i < xpdf.size(); ++i) {
        this->data()[i*extras().size_dim1 + iy] = (i <= iflv_max) ? xpdf[i] : 0.0;
      }
      // for any remaining flavours, set to zero
      for (std::size_t i=xpdf.size(); i < size_flv(); ++i) {
        this->data()[i*extras().size_dim1 + iy] = 0.0;
      }
    }
  }


  inline std::size_t size_dim0() const {return extras().size_dim0;}
  inline std::size_t size_dim1() const {return extras().size_dim1;}
  /// return the size of the flavour dimension (dim0)
  inline std::size_t size_flv() const {return size_dim0();}
  /// return the size of the y dimension (dim1)
  inline std::size_t size_y() const {return size_dim1();}

  void ensure_valid() const {
    if (!data()) {
      throw std::runtime_error("hoppet::grid_quant_2d_view::ensure_valid: data pointer is null");
    }
  }

  void assign(double val) {
    ensure_valid();
    std::fill(data(), data() + size(), val);
  } 

  grid_quant_2d_view & operator=(double value) {assign(value); return *this;}
};


//-----------------------------------------------------------------------------
/// @brief A function of y=ln(1/x) stored at fixed grid points in y for
///        each of several flavour indices

class grid_quant_2d : public data_owner<grid_quant_2d_view, grid_quant_2d_f> {

public:

  using view_t = grid_quant_2d_view;

  grid_quant_2d() {}
  grid_quant_2d(const grid_def_view & grid, std::size_t dim0_size) {
    alloc(gq2d_extras(grid, dim0_size));
  }
  void alloc_virtual(std::size_t /*sz*/, const gq2d_extras & extras_in) override {
    alloc(extras_in);
  }

  // make sure we have the move constructor, move assignment and copy assignment
  // explicit copy constructor to perform a deep copy
  grid_quant_2d(const grid_quant_2d      & other) {copy(other); }
  grid_quant_2d(const grid_quant_2d_view & other) {copy(other); }
  grid_quant_2d & operator=(const grid_quant_2d &  other) noexcept = default;
  grid_quant_2d            (      grid_quant_2d && other) noexcept = default;
  grid_quant_2d & operator=(      grid_quant_2d && other) noexcept = default;

  void alloc(gq2d_extras extras_in) {
    _extras = extras_in;
    _ptr    = hoppet_cxx__grid_quant_2d__new(extras_in.grid.ptr(), static_cast<int>(extras_in.size_dim0));
    _data   = hoppet_cxx__grid_quant_2d__data_ptr(_ptr);
    _size   = _extras.size_dim1 * extras_in.size_dim0;
  }

  grid_quant_2d & operator=(double value) {assign(value); return *this;}

};

inline grid_quant_2d operator+(grid_quant_2d a, const grid_quant_2d_view & b) {a += b; return a;}
inline grid_quant_2d operator-(grid_quant_2d a, const grid_quant_2d_view & b) {a -= b; return a;}
inline grid_quant_2d operator*(grid_quant_2d a, double b) {a *= b; return a;}
inline grid_quant_2d operator*(double b, grid_quant_2d a) {a *= b; return a;}
inline grid_quant_2d operator/(grid_quant_2d a, double b) {a /= b; return a;}
inline grid_quant_2d operator-(const grid_quant_2d_view & a) {return -1.0 * a;}


typedef grid_quant_2d_view pdf_view;
typedef grid_quant_2d      pdf;
inline  grid_quant_2d pdf_qcd(grid_def_view const & grid) {return grid_quant_2d(grid, ncompmax+1);}




//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_conv Fortran type, non-owning
class grid_conv_view : public obj_view<grid_conv_f> {
public:
  using base_t = obj_view<grid_conv_f>;
  using base_t::base_t; // ensures that constructors are inherited
  using view_t = grid_conv_view;

  grid_conv_view() = default;
  //grid_conv_view(const grid_def_view & grid) :  base_t(nullptr, grid) {}
  grid_conv_view(grid_conv_f * ptr) : base_t(ptr) {}

  grid_def_ref grid() const {return grid_def_ref(hoppet_cxx__grid_conv__grid(ptr()));}

  grid_conv_view & operator=(const grid_conv_view & other) {
    if (!ptr()) throw std::runtime_error("grid_conv_view::operator=: grid_conv_view object not associated");
    hoppet_cxx__grid_conv__copy_contents(ptr(), other.ptr());
    return *this;
  }

  /// compound assignment arithmetic operators
  ///@{
  grid_conv_view & operator+=(const grid_conv_view & other) {
    grid().ensure_compatible(other.grid());
    hoppet_cxx__grid_conv__add(_ptr, other.ptr());
    return *this;
  }
  grid_conv_view & operator-=(const grid_conv_view & other) {
    grid().ensure_compatible(other.grid());
    double minus_one = -1.0;
    hoppet_cxx__grid_conv__add(_ptr, other.ptr(), &minus_one);
    return *this;
  }
  grid_conv_view & operator*=(double factor) {
    grid().ensure_valid();
    hoppet_cxx__grid_conv__multiply(_ptr, factor);
    return *this;
  }
  grid_conv_view & operator/=(double factor) {
    grid().ensure_valid();
    hoppet_cxx__grid_conv__multiply(_ptr, 1.0/factor);
    return *this;
  }
  ///@}

  /// @brief  return an inefficiently computed (truncated) moment of the convolution operator
  /// @param n the moment index, where n=1 corresponds to momentum, n=0 to number
  /// @return the moment  
  ///
  /// Note that this is inefficient because it involves constructing
  /// a temporary grid_quant object to perform the integration and
  /// performing a full (grid.size()^2 operations) convolution.
  double truncated_moment(double n) const {
    return hoppet_cxx__grid_conv__moment(ptr(), n);
  }
};


//-----------------------------------------------------------------------------
/// @brief An object representing a convolution operator (splitting
///        function) on a grid in y=ln(1/x)
class grid_conv : public obj_owner<grid_conv_view> {
public:

  typedef obj_owner<grid_conv_view> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  grid_conv() {}

  /// construct a grid_conv object and initialise it with the given function
  ///
  /// @param grid          the grid definition
  /// @param conv_ignd_fn  the convolution integrand function, double(double y, int piece)
  /// @param split_array   optional array of points in y at which to split the adaptive integrations
  ///
  /// The function conv_ignd_fn should take two arguments:
  ///
  /// - y: the value of ln(1/x) at which to evaluate the integrand
  /// - piece: an integer identifying which sub-interval of the grid is
  ///   being integrated (one of hoppet::cc_REAL, hoppet::cc_REALVIRT, hoppet::cc_VIRT, hoppet::cc_DELTA)
  /// - return value: the value of the convolution integrand at (y, piece), including
  ///   a factor of x = exp(-y), reflecting the fact that the integration is over dy = dx/x
  ///  
  /// The split_array argument is only needed if the integrand has discontinuities
  /// or non-smoothness at specific y values
  template<typename F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double, int>>>
  grid_conv(const grid_def_ref & grid, 
            F && conv_ignd_fn, 
            const std::vector<double> & split_array = {}) : base_t(nullptr) {

    //std::cout << "grid_conv: constructing from function object, grid = " << grid.ptr() <<" " << this->grid().ptr() << "\n";
    using FuncType = decltype(conv_ignd_fn);
    // under the hood, the fortran calls hoppet_grid_conv_f__wrapper
    // (defined in hoppet_oo.cc), with a pointer to the function object
    std::function<double(double,int)> fn_ptr = std::forward<FuncType>(conv_ignd_fn);
    if (split_array.size() > 0) {
      int split_array_size = split_array.size();
      _ptr = hoppet_cxx__grid_conv__new_from_fn(grid.valid_ptr(), &fn_ptr, split_array.data(), &split_array_size);
    } else {
      _ptr = hoppet_cxx__grid_conv__new_from_fn(grid.valid_ptr(), &fn_ptr);
    }
  }
};

inline grid_quant operator*(const grid_conv_view & conv, const grid_quant_view & q) {
  // create an output grid_quant
  conv.grid().ensure_compatible(q.grid());
  grid_quant result(q.grid());

  hoppet_cxx__grid_conv__times_grid_quant(conv.ptr(), q.data(), result.data());
  return result;
}

template<typename F, typename = std::enable_if_t<std::is_invocable_r_v<double, F, double, int>>>
inline grid_conv operator*(const grid_def_view & grid, F && conv_ignd_fn) {
  return grid_conv(grid, std::forward<decltype(conv_ignd_fn)>(conv_ignd_fn));
}

// these all use a copy-and-modify strategy, exploiting the move semantics
// for copy elision
inline grid_conv operator+(grid_conv a, const grid_conv_view & b) {a += b; return a;}
inline grid_conv operator-(grid_conv a, const grid_conv_view & b) {a -= b; return a;}
inline grid_conv operator*(grid_conv a, double b) {a *= b; return a;}
inline grid_conv operator*(double b, grid_conv a) {a *= b; return a;}
inline grid_conv operator/(grid_conv a, double b) {a /= b; return a;}
inline grid_conv operator-(const grid_conv_view & a) {return -1.0 * a;}


inline grid_conv operator*(grid_conv_view const & a, grid_conv_view const & b) {
  a.grid().ensure_compatible(b.grid());
  grid_conv_f * ptr = hoppet_cxx__grid_conv__alloc_and_conv(a.ptr(), b.ptr());
  return grid_conv(ptr);
}

typedef grid_conv split_fn;
typedef grid_conv_view split_fn_view;


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the split_mat Fortran type, non-owning
class split_mat_view : public obj_view<split_mat_f> {
public:

  typedef obj_view<split_mat_f> base_t;
  //typedef grid_def_view extra_t;
  using base_t::base_t; // ensures that constructors are inherited

  split_mat_view & operator=(const split_mat_view & other) {
    if (!ptr()) throw std::runtime_error("split_mat_view::operator=: split_mat_view object not associated");
    hoppet_cxx__split_mat__copy_contents(ptr(), other.ptr());
    return *this;
  }

  const grid_def_ref grid() const { return qq().grid(); }
  int nf() const { return hoppet_cxx__split_mat__nf(ptr()); }

  /// views of the individual components of the splitting matrix
  ///@{
  grid_conv_view qq      () const { return grid_conv_view(hoppet_cxx__split_mat__qq      (ptr())); }
  grid_conv_view qg      () const { return grid_conv_view(hoppet_cxx__split_mat__qg      (ptr())); }
  grid_conv_view gq      () const { return grid_conv_view(hoppet_cxx__split_mat__gq      (ptr())); }
  grid_conv_view gg      () const { return grid_conv_view(hoppet_cxx__split_mat__gg      (ptr())); }
  grid_conv_view ns_plus () const { return grid_conv_view(hoppet_cxx__split_mat__ns_plus (ptr())); }
  grid_conv_view ns_minus() const { return grid_conv_view(hoppet_cxx__split_mat__ns_minus(ptr())); }
  grid_conv_view ns_v    () const { return grid_conv_view(hoppet_cxx__split_mat__ns_v    (ptr())); }
  ///@}


  /// compound assignment arithmetic operators
  ///@{
  split_mat_view & operator+=(const split_mat_view & other) {
    ensure_compatible(other);
    hoppet_cxx__split_mat__add(_ptr, other.ptr());
    return *this;
  }
  split_mat_view & operator-=(const split_mat_view & other) {
    ensure_compatible(other);
    double minus_one = -1.0;
    hoppet_cxx__split_mat__add(_ptr, other.ptr(), &minus_one);
    return *this;
  }
  split_mat_view & operator*=(double factor) {
    grid().ensure_valid();
    hoppet_cxx__split_mat__multiply(_ptr, factor);
    return *this;
  }
  split_mat_view & operator/=(double factor) {
    grid().ensure_valid();
    hoppet_cxx__split_mat__multiply(_ptr, 1.0/factor);
    return *this;
  }
  ///@}

  /// throws an exception if other is not compatible with *this
  void ensure_compatible(const split_mat_view & other) const {
    if (nf() != other.nf()) {
      throw std::runtime_error("hoppet::split_mat_view::ensure_compatible: incompatible split_mat nf");
    }
    grid().ensure_compatible(other.grid());
  }

};

//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the split_mat Fortran type, owning
class split_mat : public obj_owner<split_mat_view> {

public:

  typedef obj_owner<split_mat_view> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  split_mat() {}

  /// construct and allocate a split_mat object for the given number of flavours
  split_mat(int nf) {
    _ptr = hoppet_cxx__split_mat__new(nf);
  }
}; 

// these all use a copy-and-modify strategy, exploiting the move semantics
// for copy elision
inline split_mat operator+(split_mat a, const split_mat_view & b) {a += b; return a;}
inline split_mat operator-(split_mat a, const split_mat_view & b) {a -= b; return a;}
inline split_mat operator*(split_mat a, double b) {a *= b; return a;}
inline split_mat operator*(double b, split_mat a) {a *= b; return a;}
inline split_mat operator/(split_mat a, double b) {a /= b; return a;}
inline split_mat operator-(const split_mat_view & a) {return -1.0 * a;}

inline split_mat operator*(split_mat_view const & a, split_mat_view const & b) {
  a.grid().ensure_compatible(b.grid());
  split_mat_f * ptr = hoppet_cxx__split_mat__alloc_and_conv(a.ptr(), b.ptr());
  return split_mat(ptr);
}

/// Return commutator of two splitting matrices, i.e. [a,b] = a*b - b*a.
/// Note that this make use of the underlying structure of the splitting matrices
/// and is better than explicitly writing a*b - b*a
inline split_mat commutator(split_mat_view const & a, split_mat_view const & b) {
  a.grid().ensure_compatible(b.grid());
  split_mat_f * ptr = hoppet_cxx__split_mat__alloc_and_commutate(a.ptr(), b.ptr());
  return split_mat(ptr);
}


inline grid_quant_2d operator*(const split_mat_view & split, const grid_quant_2d_view & q) {
  split.grid().ensure_compatible(q.grid());
  if (q.extras().size_dim0 <= ncompmax) throw std::runtime_error("split_fn * grid_quant_2d: grid_quant_2d dim1_sz too small");
  grid_quant_2d result(q.grid(), q.extras().size_dim0);
  hoppet_cxx__split_mat__times_grid_quant_2d(split.ptr(), q.data(), result.data());
  // zero out any components beyond ncompmax
  for (size_t i = ncompmax+1; i < result.extras().size_dim0; ++i) {result[i].assign(0);}
  return result;
}

//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the mass_threshold_mat Fortran type, non-owning
class mass_threshold_mat_view : public obj_view<mass_threshold_mat_f> {
public:

  typedef obj_view<mass_threshold_mat_f> base_t;
  typedef grid_def_view extra_t;
  using base_t::base_t; // ensures that constructors are inherited

  mass_threshold_mat_view & operator=(const mass_threshold_mat_view & other) {
    if (!ptr()) throw std::runtime_error("mass_threshold_mat_view::operator=: mass_threshold_mat_view object not associated");
    hoppet_cxx__mass_threshold_mat__copy_contents(ptr(), other.ptr());
    return *this;
  }

  const grid_def_ref grid() const { return pshq().grid(); }
  int nf_heavy() const { return hoppet_cxx__mass_threshold_mat__nf_int(ptr()); }

  #define MTM_MEMBER(NAME) RETURN_OBJ_MEMBER(mass_threshold_mat,NAME,grid_conv)
  /// views of the individual components of the splitting matrix
  ///@{
  MTM_MEMBER(pshq      ) //< A^PS_Qq    Q+Qbar from singlet(nflight)
  MTM_MEMBER(pshg      ) //< A^PS_Qg    Q+Qbar from gluon  (nflight)
  MTM_MEMBER(nsqq_h    ) //< A^NS_qq,Q  ΔNS(nfheavy) from NS(nflight)
  MTM_MEMBER(sgg_h     ) //< A^S_gg,Q   Δg(nfheavy) from g(nflight)
  MTM_MEMBER(sgq_H     ) //< A^S_gq,Q   Δg(nfheavy) from singlet(nflight)
  MTM_MEMBER(psqq_h    ) //< A^PS_qq,Q  Δsinglet(nfheavy) from singlet(nflight)
  MTM_MEMBER(sqg_h     ) //< A^S_qg,Q   Δsinglet(nfheavy) from gluon(nflight)
  MTM_MEMBER(nsmqq_h   ) //< A^{NSm}_qq,Q ΔNSminus(1:nflight) from NSminus(1:nflight)
  MTM_MEMBER(pshg_msbar) //< replaces PShg when masses are MSbar (not yet supported at N3LO)
  ///@}
  #undef MTM_MEMBER

  /// compound assignment arithmetic operators
  ///@{
  mass_threshold_mat_view & operator+=(const mass_threshold_mat_view & other) {
    ensure_compatible(other);
    hoppet_cxx__mass_threshold_mat__add(_ptr, other.ptr());
    return *this;
  }
  mass_threshold_mat_view & operator-=(const mass_threshold_mat_view & other) {
    ensure_compatible(other);
    double minus_one = -1.0;
    hoppet_cxx__mass_threshold_mat__add(_ptr, other.ptr(), &minus_one);
    return *this;
  }
  mass_threshold_mat_view & operator*=(double factor) {
    grid().ensure_valid();
    hoppet_cxx__mass_threshold_mat__multiply(_ptr, factor);
    return *this;
  }
  mass_threshold_mat_view & operator/=(double factor) {
    grid().ensure_valid();
    hoppet_cxx__mass_threshold_mat__multiply(_ptr, 1.0/factor);
    return *this;
  }
  ///@}

  /// throws an exception if other is not compatible with *this
  void ensure_compatible(const mass_threshold_mat_view & other) const {
    if (nf_heavy() != other.nf_heavy()) {
      throw std::runtime_error("hoppet::mass_threshold_mat_view::ensure_compatible: incompatible mass_threshold_mat nf");
    }
    grid().ensure_compatible(other.grid());
  }

};

//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the mass_threshold_mat Fortran type, owning
class mass_threshold_mat : public obj_owner<mass_threshold_mat_view> {

public:

  typedef obj_owner<mass_threshold_mat_view> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  mass_threshold_mat() {}

  /// construct and allocate a mass_threshold_mat object for the given number of flavours
  mass_threshold_mat(int nf_heavy) {
    _ptr = hoppet_cxx__mass_threshold_mat__new(nf_heavy);
  }
}; 

// these all use a copy-and-modify strategy, exploiting the move semantics
// for copy elision
inline mass_threshold_mat operator+(mass_threshold_mat a, const mass_threshold_mat_view & b) {a += b; return a;}
inline mass_threshold_mat operator-(mass_threshold_mat a, const mass_threshold_mat_view & b) {a -= b; return a;}
inline mass_threshold_mat operator*(mass_threshold_mat a, double b) {a *= b; return a;}
inline mass_threshold_mat operator*(double b, mass_threshold_mat a) {a *= b; return a;}
inline mass_threshold_mat operator/(mass_threshold_mat a, double b) {a /= b; return a;}
inline mass_threshold_mat operator-(const mass_threshold_mat_view & a) {return -1.0 * a;}


inline grid_quant_2d operator*(const mass_threshold_mat_view & mtm, const grid_quant_2d_view & q) {
  mtm.grid().ensure_compatible(q.grid());
  if (q.extras().size_dim0 <= ncompmax) throw std::runtime_error("mass_threshold * grid_quant_2d: grid_quant_2d dim1_sz too small");
  grid_quant_2d result(q.grid(), q.extras().size_dim0);
  hoppet_cxx__mass_threshold_mat__times_grid_quant_2d(mtm.ptr(), q.data(), result.data());
  // zero out any components beyond ncompmax
  for (size_t i = ncompmax+1; i < result.extras().size_dim0; ++i) {result[i].assign(0.0);}
  return result;
}


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the running_coupling Fortran type, non-owning
///
class running_coupling_view : public obj_view<running_coupling_f> {
public:
  typedef obj_view<running_coupling_f> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  /// @brief Evaluate the running coupling at a given scale Q
  /// @param Q 
  /// @return the value of alpha_s(Q)
  double operator()(double Q) const {return hoppet_cxx__running_coupling__value(valid_ptr(), Q);}

  /// @brief Evaluate the running coupling at a given scale Q, forcing it to correspond to the specific fixed nf
  /// @param Q 
  /// @param fixnf fixed number of flavours
  /// @return the value of alpha_s(Q)
  ///
  /// Note that this may be slower than the version without fixnf, since it may
  /// bypass fast interpolation and might need to solve the evolution beyond the 
  /// cached range
  double operator()(double Q, int fixnf) const {return hoppet_cxx__running_coupling__value(valid_ptr(), Q, &fixnf);}

  /// @brief  Return the number of loops (nloop) for this running coupling
  /// @return nloop
  ///
  /// nloop=1 means the first term in the beta function, etc.
  int nloop() const { return hoppet_cxx__running_coupling__num_loops(valid_ptr()); }

  /// @brief  Return the range of active flavours (nflcl) for this running coupling
  /// @return  a tuple (nflcl_lo, nflcl_hi)
  ///
  /// use this as `auto [nflcl_lo, nflcl_hi] = rc.nf_range();`
  /// or (if nflcl_lo, nflcl_hi are already defined) `std::tie(nflcl_lo, nflcl_hi) = rc.nf_range();`
  std::tuple<int,int> nf_range() const {
    int lo=0, hi=0;
    hoppet_cxx__running_coupling__nf_range(valid_ptr(), &lo, &hi);
    return {lo, hi};
  }

  /// @brief  Return the number of active flavours at a given scale Q  
  /// @param Q 
  /// @return 
  int nf_at_Q(double Q) const {
    return hoppet_cxx__running_coupling__nf_at_q(valid_ptr(), Q);
  }

  /// @brief  Return the number of active flavours at a given scale Q, 
  ///         along with the range of Q for which this number of flavours is valid
  ///
  /// @param Q 
  /// @param Qlo is set to the lower edge of the Q range
  /// @param Qhi is set to the upper edge of the Q range
  /// @return number of active flavours at scale Q
  int nf_at_Q(double Q, double &Qlo, double & Qhi) const {
    return hoppet_cxx__running_coupling__nf_at_q(valid_ptr(), Q, &Qlo, &Qhi);
  }

  /// @brief  Return the number of active flavours at a given scale Q, 
  ///         along with the range of Q for which this number of flavours is valid
  ///

  /// @param Q 
  /// @param Qlo is set to the lower edge of the Q range
  /// @param Qhi is set to the upper edge of the Q range
  /// @param muM_mQ a non-default choice for the matching scale / quark mass ratio
  /// @return number of active flavours at scale Q
  int nf_at_Q(double Q, double &Qlo, double & Qhi, double muM_mQ) const {
    return hoppet_cxx__running_coupling__nf_at_q(valid_ptr(), Q, &Qlo, &Qhi, &muM_mQ);
  }

  /// @brief Returns the quark mass for a given flavour 
  /// @param iflv 
  /// @return 
  double quark_mass(int iflv) const { 
    // there is some ambiguity here iflv could be usual 4,5,6 or 
    // iflv_c, iflv_b, iflv_t, which are offset by iflv_g; be tolerant
    // and accept both
    int iflv_lcl = iflv > 6 ? iflv - iflv_g : iflv;
    return hoppet_cxx__running_coupling__quark_mass(valid_ptr(), iflv_lcl); 
  }

  bool quark_masses_are_msbar() const { return hoppet_cxx__running_coupling__quark_masses_are_msbar(valid_ptr()); }

  /// @brief returns the Q range for which the coupling has the corresponding number of flavours nflcl
  ///
  /// @param nflcl  the number of flavours
  /// @return       a tuple (Qlo, Qhi) giving the range in Q
  ///
  /// use this as `auto [Qlo, Qhi] = rc.Q_range_at_nf(nflcl);`
  /// or (if Qlo, Qhi are already defined) `std::tie(Qlo, Qhi) = rc.Q_range_at_nf(nflcl);`
  std::tuple<double,double> Q_range_at_nf(int nflcl) const {
    double Qlo=0.0, Qhi=0.0;
    hoppet_cxx__running_coupling__q_range_at_nf(valid_ptr(), nflcl, &Qlo, &Qhi);
    return {Qlo, Qhi};
  }

};

//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the running_coupling Fortran type, owning
///
class running_coupling : public obj_owner<running_coupling_view> {
public:
  typedef obj_owner<running_coupling_view> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  running_coupling() {}

  running_coupling(double alphas_at_Q, double Q, int nloop, int fixnf, const std::optional<double> & Qmax = std::nullopt) {
    _ptr = hoppet_cxx__running_coupling__new_fixnf(alphas_at_Q, Q, nloop, fixnf, 
                                                   Qmax.has_value() ? &(*Qmax) : nullptr);
  }

  running_coupling(double alphas_at_Q,    ///< value of alpha_s at scale Q
                   double Q,              ///< scale at which the coupling is supplied
                   int nloop,             ///< number of loops (1=LO, 2=NLO, etc)
                   double mc,             ///< charm quark mass
                   double mb,             ///< bottom quark mass
                   double mt,             ///< top quark mass
                   bool masses_are_MSbar = false,  ///< false: masses are pole masses; true: masses are MSbar
                   double muMatch_mQuark = 1.0,    ///< matching scale / quark mass ratio
                   const std::optional<double> & Qmax = std::nullopt ///< optional override of maximum Q value
                  ) {
    _ptr = hoppet_cxx__running_coupling__new_varnf(alphas_at_Q, Q, nloop,
                                                  mc, mb, mt, masses_are_MSbar, muMatch_mQuark, 
                                                  Qmax.has_value() ? &(*Qmax) : nullptr);
  }

};



//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the dglap_holder Fortran type, non-owning
class dglap_holder_view : public obj_view<dglap_holder_f> {
public:
  typedef obj_view<dglap_holder_f> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  /// @brief  set the number of active flavours (nf) for this dglap_holder
  /// @param nf  the new number of active flavours
  ///
  /// Note: this also sets the global nf variable, affecting all
  /// the global variables in the hoppet::qcd namespace
  void set_nf(int nf) {
    hoppet_cxx__dglap_holder__set_nf(valid_ptr(), nf);
  }

  /// @brief  return a view of the splitting function matrix for the given
  ///         evolution loop and number of flavours
  ///
  /// @param iloop  number of loops for the splitting function (1=LO, 2=NLO, etc)
  /// @param nf     the nf value for which to return the splitting matrix
  /// @return       a split_mat_view object corresponding to the requested splitting matrix
  inline split_mat_view p(int iloop, int nf) {
    return split_mat_view(hoppet_cxx__dglap_holder__allp(valid_ptr(), iloop, nf));
  }
  inline mass_threshold_mat_view mtm(int iloop, int nf_heavy) {
    // only have these checks when debugging is turned on
    if (iloop < 3 || iloop > nloop()) {
      throw std::runtime_error("dglap_holder_view::mtm(iloop,nf_heavy) requested iloop " + std::to_string(iloop) +
                               " outside supported range [3," + std::to_string(nloop()) + "]");
    }

    return mass_threshold_mat_view(hoppet_cxx__dglap_holder__allmtm(valid_ptr(), iloop, nf_heavy));
  }

  /// @brief return the grid definition used in this dglap_holder
  ///
  /// Note that if you use grid to construct other objects, you must
  /// ensure that those objects do not outlive the underlying
  /// dglap_holder that owns the grid.
  RETURN_OBJ_MEMBER(dglap_holder,grid,grid_def)

  /// @brief return the maximum number of loops supported by this dglap_holder
  RETURN_INT_MEMBER(dglap_holder,nloop)

  /// @brief return nf value that is currently set in this dglap_holder
  RETURN_INT_MEMBER(dglap_holder,nf)
  
  /// @brief return the factorization scheme used in this dglap_holder
  RETURN_INT_MEMBER(dglap_holder,factscheme)

  RETURN_OBJ_MEMBER(dglap_holder,p_lo,split_mat)
  RETURN_OBJ_MEMBER(dglap_holder,p_nlo,split_mat)
  RETURN_OBJ_MEMBER(dglap_holder,p_nnlo,split_mat)
  RETURN_OBJ_MEMBER(dglap_holder,p_n3lo,split_mat)

};

/// @brief Object-oriented wrapper around the dglap_holder Fortran type, owning
class dglap_holder : public obj_owner<dglap_holder_view> {
public:
  typedef obj_owner<dglap_holder_view> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  dglap_holder() {}
  dglap_holder(const hoppet::grid_def_view & grid, int factscheme, int nloop, int nflo, int nfhi) {
    _ptr = hoppet_cxx__dglap_holder__new(grid.ptr(), &factscheme, &nloop, &nflo, &nfhi);
  }
};


/// @brief  A view of the segmentation information for PDF tables
class pdfseginfo_view : public obj_view<pdfseginfo_f> {
public:
  typedef obj_view<pdfseginfo_f> base_t;
  using base_t::base_t; // ensures that constructors are inherited
  RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_lo)  ///< return the lowest  iQ index for this segment
  RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_hi)  ///< return the highest iQ index for this segment
  RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_lo)   ///< return the lowest  ln(ln(Q/Λ)) value for this segment
  RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_hi)   ///< return the highest ln(ln(Q/Λ)) value for this segment
  RETURN_DBL_MEMBER(pdfseginfo,dlnlnQ)     ///< return the spacing in ln(ln(Q/Λ)) value for this segment
  RETURN_DBL_MEMBER(pdfseginfo,inv_dlnlnQ) ///< return the inverse spacing in ln(ln(Q/Λ)) value for this segment
};

//-----------------------------------------------------------------------------
/// @brief  A view of pdf_table, i.e. a tabulation of PDFs on a (Q,iflv,y) grid
///
class pdf_table_view : public obj_view<pdf_table_f> {
public:
  typedef obj_view<pdf_table_f> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  /// @brief copy assignment operator
  ///
  /// Note that currently this does _not_ copy the cached evolution operators
  pdf_table_view & operator=(const pdf_table_view & other) {
    if (!ptr()) throw std::runtime_error("pdf_table_view::operator=: pdf_table_view object not associated");
    hoppet_cxx__pdf_table__copy_contents(ptr(), other.ptr());
    return *this;
  }

  /// @name Functions to fill the table by DGLAP evolution
  ///@{
  /// @brief Fill the PDF table by evolving from an input PDF at scale Q0
  ///
  /// Note that if (untie_nf) is true [default false] then alphas is
  /// allowed to have its natural value for nf, even if the dglap_holder
  /// cannot reach this nf value.
  void evolve(double                        Q0,        ///< scale of input PDF
              const grid_quant_2d_view &    pdf_at_Q0, ///< input PDF at scale Q0
              const dglap_holder_view &     dh,        ///< DGLAP holder to use for evolution
              const running_coupling_view & coupling,  ///< running coupling to use for evolution
              double muR_over_Q = 1.0,                 ///< ratio of renormalization to factorization scale
              int    nloop = -1,                       ///< number of loops to use; if <0, use coupling.nloop()
              bool   untie_nf = false                  ///< whether to untie nf in evolution
            ) {
    ensure_compatible(pdf_at_Q0);
    if (nloop < 0) nloop = coupling.nloop();
    hoppet_cxx__pdf_table__evolve(valid_ptr(), Q0, pdf_at_Q0.data(), dh.ptr(), coupling.ptr(),
                                 muR_over_Q, nloop, untie_nf);
  }  

  /// @brief Fill pre-evolution information for this table, to used by calling evolve()
  ///
  /// Note that if (untie_nf) is true [default false] then alphas is
  /// allowed to have its natural value for nf, even if the dglap_holder
  /// cannot reach this nf value.
  void pre_evolve(
              double                        Q0,        ///< scale of input PDF
              const dglap_holder_view &     dh,        ///< DGLAP holder to use for evolution
              const running_coupling_view & coupling,  ///< running coupling to use for evolution
              double muR_over_Q = 1.0,                 ///< ratio of renormalization to factorization scale
              int    nloop = -1,                       ///< number of loops to use; if <0, use coupling.nloop()
              bool   untie_nf = false                  ///< whether to untie nf in evolution
            ) {
    if (nloop < 0) nloop = coupling.nloop();
    hoppet_cxx__pdf_table__pre_evolve(valid_ptr(), Q0, dh.ptr(), coupling.ptr(),
                                 muR_over_Q, nloop, untie_nf);
  }  

  /// @brief Fill the PDF table by evolution of an initial condition, using pre_evolve() information
  ///
  /// @param pdf_at_Q0  the initial condition PDF at scale Q0 (as set in a pre_evolve(...) call)
  ///
  /// This is several times faster than the evolve(...) call without
  /// pre-evolution, and is the preferred option when evolving multiple
  /// initial conditions with the same evolution settings.
  void evolve(const grid_quant_2d_view &    pdf_at_Q0 //< input PDF at scale Q0
             ) {
    ensure_compatible(pdf_at_Q0);
    hoppet_cxx__pdf_table__evolve_frompre(valid_ptr(), pdf_at_Q0.data());
  }  
  ///@}


  /// @name Functions to fill the table from existing information across x,Q,iflv
  ///@{
  /// @brief assign to (i.e. fill) this table using a function fn(x,Q, xpdf_array)
  /// 
  /// @param fn Should have the signature `void(double x, double Q, double * xpdf_array)`;
  ///           xpdf_array is a double* pointing to an array with 13 entries that is 
  ///           filled by fn.  The function signature is compatible with the LHAPDF
  ///           "evolve" fortran interface;  
  ///
  /// Any flavours beyond 13 (i.e. beyond iflv_top) are set to zero.
  template<typename F, std::enable_if_t<std::is_invocable_r_v<void, F, double, double, double*>, int> = 0>
  void assign_xQ_into(const F & fn) {
    ensure_valid();
    if (size_flv() < hoppet::iflv_max+1) {
      throw std::runtime_error("grid_quant_2d_view::assign(fn): pdf_table size_flv() = " 
                  + std::to_string(size_flv()) + " < iflv_max+1 = " + std::to_string(hoppet::iflv_max+1));
    }
    std::vector<double> xvals = grid().x_values();
    std::vector<double> xpdf (size_flv());
    for (std::size_t iQ = 0 ; iQ < size_Q(); ++iQ) {
      double Q = Q_vals(iQ);
      pdf_view pdf_at_Q = at_iQ(iQ);
      for (std::size_t iy=0; iy < xvals.size(); ++iy) {
        fn(xvals[iy], Q, xpdf.data());
        for (std::size_t iflv=0; iflv < size_flv(); ++iflv) {
          pdf_at_Q(iflv,iy) = (iflv <= hoppet::iflv_max) ? xpdf[iflv] : 0.0;
        }
      }
    }
  }

  /// @brief assign to (i.e. fill) this table using a function fn(x,Q, xpdf_vector)
  /// 
  /// @param fn Should have the signature `void(double x, double Q,
  ///           std::vector<double> xpdf_vector)`; xpdf_vector is resized and
  ///           filled by fn. filled by fn.  The function signature is compatible
  ///           LHAPDF's C++ `PDF::xfxQ(x,Q,pdf_vec)` member function.
  ///
  /// Any flavours beyond the size of the xpdf_vector are set to zero.
  template<typename F, std::enable_if_t<std::is_invocable_r_v<void, F, double, double, std::vector<double> &>, int> = 0>
  void assign_xQ_into(const F & fn) {
    ensure_valid();
    if (size_flv() < hoppet::iflv_max+1) {
      throw std::runtime_error("grid_quant_2d_view::assign(fn): pdf_table size_flv() = " 
                  + std::to_string(size_flv()) + " < iflv_max+1 = " + std::to_string(hoppet::iflv_max+1));
    }
    std::vector<double> xvals = grid().x_values();
    std::vector<double> xpdf (size_flv());
    for (std::size_t iQ = 0 ; iQ < size_Q(); ++iQ) {
      double Q = Q_vals(iQ);
      pdf_view pdf_at_Q = at_iQ(iQ);
      for (std::size_t iy=0; iy < xvals.size(); ++iy) {
        fn(xvals[iy], Q, xpdf);
        // run a check only for the first y point
        if (iy == 0 && xpdf.size() > size_flv()) {
          throw std::runtime_error("grid_quant_2d_view::assign(fn): xpdf_vector size = " 
                  + std::to_string(xpdf.size()) + " > pdf_table size_flv() = " + std::to_string(size_flv()));
        }
        for (std::size_t iflv=0; iflv < xpdf.size(); ++iflv) {pdf_at_Q(iflv,iy) = xpdf[iflv];}
        for (std::size_t iflv=xpdf.size(); iflv < size_flv(); ++iflv) {pdf_at_Q(iflv,iy) = 0.0;}
      }
    }
  }
  ///@}

  /// @name Access to interpolated PDF values
  ///@{
  ///
  /// @brief return the value of xf(x,Q,iflv) at x = exp(-y)
  /// @param y     log(1/x)
  /// @param Q     the factorisation scale at which to evaluate the PDF
  /// @param iflv  the index of the flavour (in C++ numbering)
  double at_yQf(double y, double Q, int iflv) const {
    return hoppet_cxx__pdf_table__at_yqf(ptr(), y, Q, iflv);
  }

  /// @brief return the value of xf(x,Q,iflv)
  /// @param x     the momentum fraction x
  /// @param Q     the factorisation scale at which to evaluate the PDF
  /// @param iflv  the index of the flavour (in C++ numbering)
  double at_xQf(double x, double Q, int iflv) const {return at_yQf(-std::log(x), Q, iflv);}

  /// @brief fill the result array with the values of xf(x,Q,:) at the given y and Q
  /// @param y       y = log(1/x)
  /// @param Q       the factorisation scale at which to evaluate the PDF
  /// @param result  pointer to an array of size at least iflv_max()+1 to be filled with the results
  void at_yQ_into(double y, double Q, double * result) const {
    hoppet_cxx__pdf_table__at_yq_into(ptr(), y, Q, result);
  }

  /// @brief fill the result array with the values of xf(x,Q,:) 
  /// @param x       the momentum fraction x
  /// @param Q       the factorisation scale at which to evaluate the PDF
  /// @param result  pointer to an array of size at least iflv_max()+1 to be filled with the results
  void at_xQ_into(double x, double Q, double * result) const {
    at_yQ_into(-std::log(x), Q, result);
  }

  /// @brief return a grid_quant_2d (i.e. full set of pdf flavours across y) at the specified Q
  /// @param Q 
  /// @return the grid_quant_2d at the specified Q
  grid_quant_2d at_Q(double Q) const {
    grid_quant_2d result(grid(), size_flv());
    hoppet_cxx__pdf_table__at_q_into(valid_ptr(), Q, result.data());
    return result;
  }
  ///@}

  
  /// @name Direct access to the table data
  ///@{
  //-------------------------------------------------------------

  /// @brief return a grid_quant_2d_view (i.e. pdf slice) at the specified iQ
  ///
  /// @param iQ     the index in the Q dimension
  /// @return       a view of the table's full (y,iflv) structure at iQ
  ///
  /// Note that the view can be read from and written to. In the latter case, this
  /// effectively (silently) discards the function's const qualifier
  grid_quant_2d_view at_iQ(size_t iQ) const {
    double * tab_ptr = hoppet_cxx__pdf_table__tab_ptr(valid_ptr());
    // these are the sizes of the two dimensions of result (shifted by 1 index wrt table)
    size_t size_dim0 = size_flv();
    size_t size_dim1 = static_cast<size_t>(grid().size());
    size_t iQ_size = size_dim0 * size_dim1;
    double * iQ_ptr  = tab_ptr + iQ * iQ_size;
    return grid_quant_2d_view(iQ_ptr, iQ_size, gq2d_extras(grid(), size_dim0));
  }

  /// @brief return a grid_quant_view (i.e. pdf_flav slice) at the specified iQ, iflv
  ///
  /// @param iQ     the index in the Q dimension
  /// @param iflv   the index in the flavour dimension
  /// @return       a view of the table's full (y) structure at iQ,iflv
  ///
  /// Note that the view can be read from and written to. In the latter case, this
  /// effectively (silently) discards the function's const qualifier
  grid_quant_view at_iQf(size_t iQ, size_t iflv) const {
    double * tab_ptr = hoppet_cxx__pdf_table__tab_ptr(valid_ptr());
    // these are the sizes of the two dimensions of result (shifted by 1 index wrt table)
    size_t size_dim0 = size_flv();
    size_t size_dim1 = static_cast<size_t>(grid().size());
    size_t iQ_size = size_dim0 * size_dim1;
    double * iQf_ptr  = tab_ptr + iQ * iQ_size + iflv * size_dim1;
    return grid_quant_view(iQf_ptr, size_dim1, grid());
  }


  /// @brief alias for at_iQ(iQ)
  grid_quant_2d_view operator[](size_t iQ) const {return at_iQ(iQ);}
  /// @brief alias for at_iQ(iQ)
  grid_quant_2d_view operator()(size_t iQ) const {return at_iQ(iQ);}
  /// @brief alias for at_iQf(iQ,iflv)
  grid_quant_view operator()(size_t iQ, size_t iflv) const {return at_iQf(iQ,iflv);}

  //-------------------------------------------------------------
  ///@}

  /// @name Querying the table properties
  ///@{
  ///
  /// @brief return the maximum valid iflv index, in C++ numbering; 
  size_t iflv_max() const {return static_cast<size_t>(tab_iflv_max() - iflv_min_fortran);  }

  /// @brief return the size of the flavour (dim1=2nd) dimension of the table 
  ///
  /// For pure QCD tables, this is iflv_max() + 2, because it includes
  /// hoppet::iflv_info, which is at index iflv_max() + 1. 
  size_t size_flv() const {return static_cast<size_t>(hoppet_cxx__pdf_table__size_flv(valid_ptr())); }

  /// return the size of the Q (dim0=1st) dimension of the table
  size_t size_Q() const {return nQ()+1;}


  /// @brief write the table to a file in LHAPDF(6) grid format
  void write_lhapdf(const running_coupling_view & coupling, 
                    const std::string &           basename,
                    int                           pdf_index = 0,
                    const std::optional<int> &    iy_increment = std::nullopt,
                    const std::optional<std::vector<int>> &  flav_indices = std::nullopt,
                    const std::optional<std::vector<int>> &  flav_pdg_ids = std::nullopt,
                    const std::optional<std::vector<double>> & flav_rescale = std::nullopt
                   ) const {
    auto opt_sclptr = [](const auto & optin) {return optin ? &(optin.value())     : nullptr;};
    auto opt_vecptr = [](const auto & optin) {return optin ? optin.value().data() : nullptr;};

    int n_flav = flav_indices ? flav_indices->size() : 0;

    hoppet_cxx__pdf_table__write_lhapdf(this->valid_ptr(), coupling.valid_ptr(),
                                        basename.c_str(), pdf_index,
                                        opt_sclptr(iy_increment),
                                        flav_indices ? & n_flav : nullptr,
                                        opt_vecptr(flav_indices),
                                        opt_vecptr(flav_pdg_ids),
                                        opt_vecptr(flav_rescale)
                                       );
  }

  /// @brief  throw an error if some_pdf is not compatible with *this (grid & flavour dimension size)
  /// @param some_pdf 
  void ensure_compatible(const grid_quant_2d_view & some_pdf) const {
    grid().ensure_compatible(some_pdf.grid());
    if (some_pdf.size_flv() != size_flv()) {
      throw std::runtime_error(
        "hoppet::pdf_table_view::ensure_compatible: incompatible pdf_table.size_flv()="+std::to_string(size_flv())+
        " vs pdf_at_Q0.size_flv()="+std::to_string(some_pdf.size_flv()));
    }
  }

  // think carefully which of these interfaces should be public, which perhaps renamed
  RETURN_OBJ_MEMBER(pdf_table,grid,grid_def)      ///< return the table's grid definition
  RETURN_INT_MEMBER(pdf_table,nQ)                 ///< return nQ, the highest Q index (size_Q() = nQ()+1)
  RETURN_INT_MEMBER(pdf_table,lnlnQ_order)        ///< return the order of the ln(ln(Q)) interpolation
  RETURN_LOG_MEMBER(pdf_table,freeze_at_Qmin)     ///< if true, access functions return pdf at Qmin for Q<Qmin; otherwise they return zero
  RETURN_LOG_MEMBER(pdf_table,nf_info_associated) ///< returns true if the table has nf info associated
  RETURN_INT_MEMBER(pdf_table,nflo)               ///< return table's lowest nf value
  RETURN_INT_MEMBER(pdf_table,nfhi)               ///< return table's highest nf value
  RETURN_DBL_MEMBER(pdf_table,lnlnQ_min)          ///< return the minimum ln(ln(Q/Λ)) value of the table
  RETURN_DBL_MEMBER(pdf_table,lnlnQ_max)          ///< return the maximum ln(ln(Q/Λ)) value of the table
  RETURN_DBL_MEMBER(pdf_table,lambda_eff)         ///< return the effective Λ value of the table
  RETURN_DBL_MEMBER(pdf_table,dlnlnQ)             ///< return the spacing in ln(ln(Q/Λ)) of the table
  RETURN_OBJ_MEMBER(pdf_table,seginfo_no_nf,pdfseginfo)  ///< return the pdfseginfo object when there is no nf info
  RETURN_OBJ_MEMBER_I(pdf_table,seginfo,pdfseginfo) ///< return the pdfseginfo object for segment with nf == i
  RETURN_DBL_MEMBER_I(pdf_table,as2pi)            ///< return \f$\alpha_s(Q_i)/2\pi\f$ where \f$Q_i\f$ is Q_vals(i)
  RETURN_INT_MEMBER_I(pdf_table,nf_int)           ///< return the number of active flavours at Q index i
  RETURN_DBL_MEMBER_I(pdf_table,lnlnQ_vals)       ///< return the ln(ln(Q/Λ)) value at Q index i
  RETURN_DBL_MEMBER_I(pdf_table,Q_vals)           ///< return the Q value at Q index i

  ///
  ///@}


protected:
  /// @brief return the maximum flavour number for which the table has been allocate
  /// This is in Fortran numbering, which is why it is protected; use iflv_max() instead
  RETURN_INT_MEMBER(pdf_table,tab_iflv_max)

};

//-----------------------------------------------------------------------------
/// @brief  Holds a tabulation of a pdf. Most functions
///         are inherited from and documented in pdf_table_view
class pdf_table : public obj_owner<pdf_table_view> {
public:
  typedef obj_owner<pdf_table_view> base_t;
  using base_t::base_t; // ensures that constructors are inherited  

  /// @brief  default constructor
  pdf_table() {}

  /// @brief construct and allocate a pdf_table object
  pdf_table(const grid_def_view & grid,    ///< grid definition on which to build the table
            double Qmin,                   ///< minimum Q value of the table
            double Qmax,                   ///< maximum Q value of the table  
            double dlnlnQ,                 ///< spacing in ln(ln(Q/Λ)) of the table
            int lnlnQ_order = hoppet_pdf_table_def_lnlnQ_order, ///< order of ln(ln(Q)) interpolation
            bool freeze_at_Qmin = false,   ///< if true, access functions return pdf at Qmin for Q<Qmin; otherwise they return zero  
            int iflv_max_table = ncompmax  ///< maximum iflv index to allocate in the table (including info index)
          ) {
    _ptr = hoppet_cxx__pdf_table__new(grid.ptr(), Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table);
  }

  /// @brief  Associate nf information with this table, to allow for segmentation across different nf values
  /// @param coupling a variable-flavour-number coupling, whose flavour thresholds will be used to segment the table
  ///
  /// This function reallocates the table storage to allow for the segmentation
  /// in nf, so it should be called before filling the table.
  void add_nf_info(const running_coupling & coupling) {
    hoppet_cxx__pdf_table__add_nf_info(valid_ptr(), coupling.ptr());
  }
};


//-----------------------------------------------------------------------------
/// @brief  A view of an array of pdf_table objects
class pdf_table_array_view : public data_view<Empty, pdf_table_f> {
public:
  typedef data_view<Empty, pdf_table_f> base_t;
  using base_t::base_t; // ensures that constructors are inherited

  pdf_table_array_view() noexcept {}
  pdf_table_array_view(pdf_table_f * ptr, std::size_t sz): base_t(ptr, sz, Empty()) {}

  /// @brief return a view of the i-th pdf_table in the array (numbered from 0)
  pdf_table_view operator[](std::size_t i) {
    return pdf_table_view(hoppet_cxx__pdf_tables__table_i(this->data(), this->size(), i));
  }
};

//-----------------------------------------------------------------------------
/// @brief  An object that holds an array of pdf_table objects
///
/// The main advantage of this class (and its view) over a std::vector
/// of pdf_tables, is that it gives access to efficient operations on the
/// whole array of tables, notably interpolation in x,Q across all
/// tables in one go. These require an underlying Fortran array of
/// pdf_table objects, which is what this class effectively wraps.
///
/// To use this class, one first constructs it with the desired size
/// (number of pdf_table objects), and then allocates and fills each
/// table individually, e.g. using the operator[] to get a view of each
/// table in turn.
///
/// ```c++
/// hoppet::pdf_table_array tables(3); // array of 3 pdf_table objects
/// tables[0] = hoppet::pdf_table(...); // allocate table 0
/// tables[1] = hoppet::pdf_table(...); // allocate table 1
/// tables[2] = hoppet::pdf_table(...); // allocate table 2
/// ```
class pdf_table_array: public data_owner<pdf_table_array_view, pdf_table_array_f> {
public:
  typedef data_owner<pdf_table_array_view, pdf_table_array_f> base_t;
  using base_t::base_t; // ensures that constructors are inherited
  typedef pdf_table_array_view view_t;

  pdf_table_array(std::size_t sz) {alloc(sz);}

  void alloc_virtual(std::size_t sz, const Empty & ) override {alloc(sz);}
  void alloc(std::size_t sz) {
    _size = sz;
    hoppet_cxx__pdf_table_array__new(sz, &_ptr, &_data);
  }

};

} // end namespace hoppet -------------------------------------



/// objects globally defined in the streamlined interface
namespace hoppet {

/// Namespace for access to the globally defined objects from the streamlined interface  
///
/// The objects here are automatically set up when calling hoppetStart
/// (or hoppetStartExtended) followed by a call hoppetEvolve.
///
namespace sl {
  /// a view of the grid_def object being used in the streamlined interface
  extern grid_def_view grid;

  /// a view of the dglap_holder object being used in the streamlined interface
  extern dglap_holder_view dh;

  /// a view of the running_coupling object being used in the streamlined interface
  extern running_coupling_view coupling; 

  /// @brief  a view of the main pdf_table object being used in the streamlined interface (same as hoppet::sl::tables[0])
  extern pdf_table_view    table; 

  /// @brief  a view of the array of pdf_table objects in the streamlined interface
  extern pdf_table_array_view  tables; 

}
}

#undef DEFINE_RETURN_INT_MEMBER
#undef DEFINE_RETURN_DBL_MEMBER
#undef DEFINE_RETURN_OBJ_MEMBER
#undef DEFINE_RETURN_OBJ_MEMBER_I 
#undef DEFINE_RETURN_OBJ_MEMBER_IJ
#undef DEFINE_DELETE
#undef DEFINE_COPY

#undef RETURN_INT_MEMBER
#undef RETURN_DBL_MEMBER
#undef RETURN_OBJ_MEMBER
#undef RETURN_OBJ_MEMBER_I 
#undef RETURN_OBJ_MEMBER_IJ


#endif // __HOPPET_OO__