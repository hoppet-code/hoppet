#ifndef __HOPPET_OO__
#define __HOPPET_OO__
#include "hoppet.h"
#include "hoppet/base_types.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <tuple>
#include <functional>
#include <concepts>


// Elements to think about
// - [ ] do we separate things out into different files
// - [ ] do we provide "physics" aliases: e.g. grid_quant -> pdf_flav, grid_conv -> split_fn
// - [ ] make sure hoppet/qcd.h (etc) gets installed properly
// - [ ] think about hoppet__qcd v hoppet_cxx__qcd naming conventions
// - [ ] making sure we treat default epsilon correctly

// Next steps:
// - [x] add the running_coupling class
// - [x] add the dglap_holder    
// - [~] add the pdf_table
//       - [ ] grid_quant_3d
//       - [ ] evolution functions
//       - [ ] PDF access functions
// - [ ] add the evln_operator class
// - [x] add the mass thresholds 
// - [ ] grid_quant_3d?
// - [~] add streamlined interface functions
// - [ ] support for probes
// - [ ] add qed support, including pdf_table allocators
// - [x] add access to things like the beta function coefficients, qcd group constants, etc.
// - [ ] add structure function support
// - [ ] add documentation

// Things to perhaps add
// - [x] add take_view to the obj_view class
// - [x] basic checks of grid compatibility, so as to get C++ errors rather than Fortran errors
// - [ ] mvv interface for splitting functions?



#define DEFINE_RETURN_INT_MEMBER(classname, membername)           extern "C" int            hoppet_cxx__##classname##__##membername(const classname##_f * ptr);
#define DEFINE_RETURN_LOG_MEMBER(classname, membername)           extern "C" bool           hoppet_cxx__##classname##__##membername(const classname##_f * ptr);
#define DEFINE_RETURN_DBL_MEMBER(classname, membername)           extern "C" double         hoppet_cxx__##classname##__##membername(const classname##_f * ptr);
#define DEFINE_RETURN_INT_MEMBER_I(classname, membername)         extern "C" int            hoppet_cxx__##classname##__##membername(const classname##_f * ptr, int i);
#define DEFINE_RETURN_DBL_MEMBER_I(classname, membername)         extern "C" double         hoppet_cxx__##classname##__##membername(const classname##_f * ptr, int i);
#define DEFINE_RETURN_OBJ_MEMBER(classname, membername, typename) extern "C" typename##_f * hoppet_cxx__##classname##__##membername(const classname##_f * ptr);
#define DEFINE_RETURN_OBJ_MEMBER_I( classname, membername, typename) extern "C" typename##_f * hoppet_cxx__##classname##__##membername(const classname##_f * ptr, int i);
#define DEFINE_RETURN_OBJ_MEMBER_IJ(classname, membername, typename) extern "C" typename##_f * hoppet_cxx__##classname##__##membername(const classname##_f * ptr, int i, int j);
#define DEFINE_DELETE(classname) \
   extern "C" void hoppet_cxx__##classname##__delete(classname##_f ** ptr); \
   inline void generic_delete(classname##_f * ptr) {if (ptr) hoppet_cxx__##classname##__delete(&ptr);}
#define DEFINE_COPY(classname) \
   extern "C" classname##_f * hoppet_cxx__##classname##__copy(const classname##_f * ptr); \
   inline classname##_f * generic_copy(const classname##_f * ptr) {return hoppet_cxx__##classname##__copy(ptr);}

#define RETURN_INT_MEMBER(classname, membername)           inline int    membername() const {return hoppet_cxx__##classname##__##membername(valid_ptr());} 
#define RETURN_LOG_MEMBER(classname, membername)           inline bool   membername() const {return hoppet_cxx__##classname##__##membername(valid_ptr());} 
#define RETURN_DBL_MEMBER(classname, membername)           inline double membername() const {return hoppet_cxx__##classname##__##membername(valid_ptr());} 
#define RETURN_INT_MEMBER_I(classname, membername)         inline int    membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_DBL_MEMBER_I(classname, membername)         inline double membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_OBJ_MEMBER(classname, membername, typename) inline typename##_view membername() const {return typename##_view(hoppet_cxx__##classname##__##membername(valid_ptr()));} 
#define RETURN_OBJ_MEMBER_I( classname, membername, typename) inline typename##_view membername(int i) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i);} 
#define RETURN_OBJ_MEMBER_IJ(classname, membername, typename) inline typename##_view membername(int i, int j) const {return hoppet_cxx__##classname##__##membername(valid_ptr(), i, j);} 




/// "forward" declaration of the Fortran grid_def type;
/// note that we only ever use pointers to it, so we do not
/// need to know anything about its actual structure, which
/// is hidden in the Fortran code.
class grid_def_f;
class grid_quant_f;
class grid_quant_2d_f;
class grid_conv_f;
class split_mat_f;
class running_coupling_f;
class dglap_holder_f;
class mass_threshold_mat_f;
// things for pdf_table
class pdf_table_f;
class pdfseginfo_f;
class evln_operator_f;



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
DEFINE_RETURN_DBL_MEMBER(grid_def,dy)
DEFINE_RETURN_DBL_MEMBER(grid_def,ymax)
DEFINE_RETURN_DBL_MEMBER(grid_def,eps)
DEFINE_RETURN_INT_MEMBER(grid_def,nsub)
DEFINE_RETURN_INT_MEMBER(grid_def,order)
inline void generic_delete(grid_def_f * ptr) {if (ptr) hoppet_cxx__grid_def__delete(&ptr);}
inline grid_def_f * generic_copy(const grid_def_f * ptr) {  return hoppet_cxx__grid_def__copy(ptr);}

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
  grid_conv_f * hoppet_cxx__grid_conv__new_from_fn(const grid_def_f * grid_ptr, void * conv_ignd_c_fn_obj);
  grid_conv_f * hoppet_cxx__grid_conv__new_from_gc(const grid_conv_f * gc_other);
  void  hoppet_cxx__grid_conv__copy_contents(grid_conv_f * dest, const grid_conv_f * src); //< src copied into dest
  //void  hoppet_cxx__grid_conv__delete(grid_conv_f ** gridconv);
  void  hoppet_cxx__grid_conv__times_grid_quant(const grid_conv_f * conv, const double * q_in_data, double * q_out_data);

  void hoppet_cxx__grid_conv__add(grid_conv_f * conv, const grid_conv_f * other, const double * factor = nullptr);
  void hoppet_cxx__grid_conv__multiply(grid_conv_f * conv, const double factor);
  grid_conv_f * hoppet_cxx__grid_conv__alloc_and_conv(const grid_conv_f * conv1, const grid_conv_f * conv2);
}
DEFINE_DELETE(grid_conv)
DEFINE_RETURN_OBJ_MEMBER(grid_conv,grid,grid_def)
//inline void generic_delete(grid_conv_f * ptr) {if (ptr)hoppet_cxx__grid_conv__delete(&ptr);}
inline grid_conv_f * generic_copy(const grid_conv_f * ptr) {return hoppet_cxx__grid_conv__new_from_gc(ptr);}
//  if (ptr)  return hoppet_cxx__grid_conv__new_from_gc(ptr); else return nullptr;}


/// split_mat function wrappers
extern "C" {
  split_mat_f * hoppet_cxx__split_mat__new(int nf);
  split_mat_f * hoppet_cxx__split_mat__copy(const split_mat_f * splitmat);
  void hoppet_cxx__split_mat__copy_contents(split_mat_f * dest, const split_mat_f * src); //< src copied into dest
  void hoppet_cxx__split_mat__delete(split_mat_f ** splitmat);

  void hoppet_cxx__split_mat__add(split_mat_f * conv, const split_mat_f * other, const double * factor = nullptr);
  void hoppet_cxx__split_mat__multiply(split_mat_f * conv, const double factor);
  split_mat_f * hoppet_cxx__split_mat__times_grid_quant_2d (const split_mat_f * splitmat, const double * q_in_data, double * q_out_data);
  split_mat_f * hoppet_cxx__split_mat__alloc_and_conv(const split_mat_f * a, const split_mat_f * b);
  split_mat_f * hoppet_cxx__split_mat__alloc_and_commutate(const split_mat_f * a, const split_mat_f * b);

  int hoppet_cxx__split_mat__nf(const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__qq      (const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__qg      (const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__gq      (const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__gg      (const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__ns_plus (const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__ns_minus(const split_mat_f * splitmat);
  grid_conv_f * hoppet_cxx__split_mat__ns_v    (const split_mat_f * splitmat);
}
inline void generic_delete(split_mat_f * ptr) {if (ptr) hoppet_cxx__split_mat__delete(&ptr);}
inline split_mat_f * generic_copy(const split_mat_f * ptr) {return hoppet_cxx__split_mat__copy(ptr);}

/// mass_threshold_mat function wrappers
extern "C" {
  mass_threshold_mat_f * hoppet_cxx__mass_threshold_mat__new(int nf_heavy);
  void hoppet_cxx__mass_threshold_mat__set_nf(mass_threshold_mat_f * ptr, int nf_lcl);
  void hoppet_cxx__mass_threshold_mat__multiply(mass_threshold_mat_f * ptr, const double factor);
  void hoppet_cxx__mass_threshold_mat__add     (mass_threshold_mat_f * ptr, const mass_threshold_mat_f * other, double * factor = nullptr);
  mass_threshold_mat_f * hoppet_cxx__mass_threshold_mat__times_grid_quant_2d(const mass_threshold_mat_f * mtm, const double * q_in_data, double * q_out_data);
  void hoppet_cxx__mass_threshold_mat__copy_contents(mass_threshold_mat_f * dest, const mass_threshold_mat_f * src); //< src copied into dest
}
DEFINE_COPY(mass_threshold_mat)
DEFINE_DELETE(mass_threshold_mat)
DEFINE_RETURN_INT_MEMBER(mass_threshold_mat,nf_int)

#define MTM_REF(NAME)  DEFINE_RETURN_OBJ_MEMBER(mass_threshold_mat,NAME,grid_conv)
MTM_REF(pshq      ) //< A^PS_Qq    Q+Qbar from singlet(nflight)
MTM_REF(pshg      ) //< A^PS_Qg    Q+Qbar from gluon  (nflight)
MTM_REF(nsqq_h    ) //< A^NS_qq,Q  ΔNS(nfheavy) from NS(nflight)
MTM_REF(sgg_h     ) //< A^S_gg,Q   Δg(nfheavy) from g(nflight)
MTM_REF(sgq_H     ) //< A^S_gq,Q   Δg(nfheavy) from singlet(nflight)
MTM_REF(psqq_h    ) //< A^PS_qq,Q  Δsinglet(nfheavy) from singlet(nflight)
MTM_REF(sqg_h     ) //< A^S_qg,Q   Δsinglet(nfheavy) from gluon(nflight)
MTM_REF(nsmqq_h   ) //< A^{NSm}_qq,Q ΔNSminus(1:nflight) from NSminus(1:nflight)
MTM_REF(pshg_msbar) //< replaces PShg when masses are MSbar (not yet supported at N3LO)
#undef MTM_REF


/// running_coupling function wrappers
extern "C" {
  running_coupling_f * hoppet_cxx__running_coupling__new_fixnf(double alphaS, double Q, int nloop, int fixnf);
  running_coupling_f * hoppet_cxx__running_coupling__new_varnf(double alphaS, double Q, int nloop,
                                                              double mc, double mb, double mt,
                                                              bool masses_are_MSbar, double muMatch_mQuark);
  void hoppet_cxx__running_coupling__delete(running_coupling_f ** rc);
  double hoppet_cxx__running_coupling__value(const running_coupling_f * rc, double Q, int * fixnf = 0);
  int    hoppet_cxx__running_coupling__num_loops(const running_coupling_f * rc);
  void   hoppet_cxx__running_coupling__nf_range(const running_coupling_f * rc, int * nflo, int * nfhi);
  int    hoppet_cxx__running_coupling__nf_at_q(const running_coupling_f * rc, double Q, double * Qlo=0, double * Qhi=0, const double * muM_mQ=0);
  double hoppet_cxx__running_coupling__quark_mass(const running_coupling_f * rc, int iflv);
  bool   hoppet_cxx__running_coupling__quark_masses_are_msbar(const running_coupling_f * rc);
  void   hoppet_cxx__running_coupling__q_range_at_nf(const running_coupling_f * rc, int nflcl, double * Qlo, double * Qhi);
}
inline void generic_delete(running_coupling_f * ptr) {if (ptr) hoppet_cxx__running_coupling__delete(&ptr);}

// things for the dglap_holder
extern "C" {
  dglap_holder_f * hoppet_cxx__dglap_holder__new(const grid_def_f * gridptr, const int * factscheme, const int * nloop, const int * nflo, const int * nfhi);
  void hoppet_cxx__dglap_holder__set_nf(dglap_holder_f * ptr, int nf);
}
DEFINE_DELETE(dglap_holder)
DEFINE_RETURN_INT_MEMBER(dglap_holder,nloop)
DEFINE_RETURN_INT_MEMBER(dglap_holder,nf)
DEFINE_RETURN_INT_MEMBER(dglap_holder,factscheme)
DEFINE_RETURN_OBJ_MEMBER_IJ(dglap_holder,allp,split_mat)
DEFINE_RETURN_OBJ_MEMBER_IJ(dglap_holder,allmtm,mass_threshold_mat)

DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_lo,split_mat)
DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_nlo,split_mat)
DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_nnlo,split_mat)
DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_n3lo,split_mat)

//----- things for pdfseginfo ----------
DEFINE_RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_lo)
DEFINE_RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_hi)
DEFINE_RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_lo)
DEFINE_RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_hi)
DEFINE_RETURN_DBL_MEMBER(pdfseginfo,dlnlnQ)
DEFINE_RETURN_DBL_MEMBER(pdfseginfo,inv_dlnlnQ)

//----- things for pdf_table ----------
extern "C" {
  extern int hoppet_pdf_table_def_lnlnQ_order; //< default lnlnQ order for pdf_table
  pdf_table_f * hoppet_cxx__pdf_table__new(const grid_def_f * grid, double Qmin, double Qmax, 
                    double dlnlnQ, int lnlnQ_order, bool freeze_at_Qmin, int iflv_max_table); 
  void hoppet_cxx__pdf_table__copy_contents(pdf_table_f * dest, const pdf_table_f * src); //< src copied into dest                  
  void hoppet_cxx__pdf_table__add_nf_info(pdf_table_f * tab, const running_coupling_f * coupling);
  double * hoppet_cxx__pdf_table__tab_ptr(const pdf_table_f * tab);
  int      hoppet_cxx__pdf_table__size_flv(const pdf_table_f * tab);
}
DEFINE_COPY(pdf_table)
DEFINE_DELETE(pdf_table)

DEFINE_RETURN_OBJ_MEMBER(pdf_table,grid,grid_def)
DEFINE_RETURN_INT_MEMBER(pdf_table,nQ)
DEFINE_RETURN_INT_MEMBER(pdf_table,tab_iflv_max)
DEFINE_RETURN_INT_MEMBER(pdf_table,lnlnQ_order)
DEFINE_RETURN_LOG_MEMBER(pdf_table,freeze_at_Qmin)
DEFINE_RETURN_LOG_MEMBER(pdf_table,nf_info_associated)
DEFINE_RETURN_INT_MEMBER(pdf_table,nflo)
DEFINE_RETURN_INT_MEMBER(pdf_table,nfhi)
DEFINE_RETURN_DBL_MEMBER(pdf_table,lnlnQ_min)
DEFINE_RETURN_DBL_MEMBER(pdf_table,lnlnQ_max)
DEFINE_RETURN_DBL_MEMBER(pdf_table,lambda_eff)
DEFINE_RETURN_DBL_MEMBER(pdf_table,dlnlnQ)
DEFINE_RETURN_OBJ_MEMBER(pdf_table,seginfo_no_nf,pdfseginfo)

DEFINE_RETURN_OBJ_MEMBER_I(pdf_table,seginfo,pdfseginfo)
DEFINE_RETURN_DBL_MEMBER_I(pdf_table,as2pi)
DEFINE_RETURN_INT_MEMBER_I(pdf_table,nf_int)
DEFINE_RETURN_DBL_MEMBER_I(pdf_table,lnlnQ_vals)
DEFINE_RETURN_DBL_MEMBER_I(pdf_table,Q_vals)

// DEFINE_RETURN_OBJ_MEMBER(pdf_table,grid,grid_def)
// DEFINE_RETURN_INT_MEMBER(pdf_table,nq)
// DEFINE_RETURN_INT_MEMBER(pdf_table,tab_iflv_max)
// DEFINE_RETURN_INT_MEMBER(pdf_table,lnlnq_order)
// DEFINE_RETURN_LOG_MEMBER(pdf_table,freeze_at_qmin)
// DEFINE_RETURN_LOG_MEMBER(pdf_table,nf_info_associated)
// DEFINE_RETURN_INT_MEMBER(pdf_table,nflo)
// DEFINE_RETURN_INT_MEMBER(pdf_table,nfhi)
// DEFINE_RETURN_DBL_MEMBER(pdf_table,lnlnq_min)
// DEFINE_RETURN_DBL_MEMBER(pdf_table,lnlnq_max)
// DEFINE_RETURN_DBL_MEMBER(pdf_table,lambda_eff)
// DEFINE_RETURN_DBL_MEMBER(pdf_table,dlnlnq)
// DEFINE_RETURN_OBJ_MEMBER(pdf_table,seginfo_no_nf,pdfseginfo)
// DEFINE_RETURN_OBJ_MEMBER_I(pdf_table,seginfo,pdfseginfo)
// DEFINE_RETURN_DBL_MEMBER_I(pdf_table,as2pi)
// DEFINE_RETURN_INT_MEMBER_I(pdf_table,nf_int)
// DEFINE_RETURN_INT_MEMBER_I(pdf_table,lnlnq_vals)
// DEFINE_RETURN_INT_MEMBER_I(pdf_table,q_vals)


/// namespace hoppet contains the object-oriented C++ interface to HOPPET
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
template <typename F>
concept DoubleFnDouble =
    std::invocable<F, double> &&
    std::same_as<std::invoke_result_t<F, double>, double>;

template <typename F>
concept DoubleFnDoubleInt =
    std::invocable<F, double, int> &&
    std::same_as<std::invoke_result_t<F, double, int>, double>;
template <typename F>
concept VoidFnDoubleDoubleDoubleptr =
    std::invocable<F, double, double, double *> &&
    std::same_as<std::invoke_result_t<F, double, double, double *>, void >;



//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_def Fortran type, non-owning
///
/// This version provides a "view" onto an existing Fortran grid_def object,
/// without taking ownership of it.
class grid_def_view : public obj_view<grid_def_f> {
public:

  typedef obj_view<grid_def_f> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  int ny() const {return hoppet_cxx__grid_def__ny(valid_ptr()); }

  std::vector<double> y_values() const {
    std::vector<double> yvals(ny()+1);
    hoppet_cxx__grid_def__y_values(valid_ptr(), yvals.data());
    return yvals;
  }
  std::vector<double> x_values() const {
    std::vector<double> xvals(ny()+1);
    hoppet_cxx__grid_def__x_values(valid_ptr(), xvals.data());
    return xvals;
  }

  inline bool operator==(const grid_def_view & other) const noexcept {
    // equivalence is only true if both pointers are non-null
    // they are either the same pointer or equivalent grids
    return ptr() && other.ptr() && (
      ptr() == other.ptr() || hoppet_cxx__grid_def__equiv(ptr(), other.ptr())
    );
  }
  inline bool operator!=(const grid_def_view & other) const noexcept {
    return !(*this == other);
  }

  void ensure_compatible(const grid_def_view & other) const {
    if (*this != other) {
      throw std::runtime_error(
        "hoppet::grid_def_view::ensure_compatible: grids are not equivalent");
    }
  }

  RETURN_DBL_MEMBER(grid_def,dy)
  RETURN_DBL_MEMBER(grid_def,ymax)
  RETURN_DBL_MEMBER(grid_def,eps)
  RETURN_INT_MEMBER(grid_def,nsub)
  RETURN_INT_MEMBER(grid_def,order)

};


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_def Fortran type, with ownership
///
/// This version takes ownership of the underlying Fortran grid_def object
class grid_def : public obj_owner<grid_def_view> {
public:

  typedef obj_owner<grid_def_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  /// @brief construct and allocate a new grid_def object
  ///
  /// @param dy      grid spacing
  /// @param ymax    maximum y value
  /// @param order   usual interpolation order parameter
  /// @param eps     accuracy for the adaptive integration
  //
  explicit grid_def(double dy, double ymax, int order=-5, double eps=1e-7)
    : grid_def(hoppet_cxx__grid_def__new(dy, ymax, order, eps)) {}

  /// @brief construct a grid_def from multiple grid_def objects
  ///
  /// @param grids   vector of grid_def objects
  /// @param locked  whether the new grid should be "locked"
  explicit grid_def(const std::vector<grid_def> & grids, bool locked=false) {
    int ngrids = static_cast<int>(grids.size());
    std::vector<const grid_def_f *> grid_ptrs(ngrids);
    for (int i=0; i<ngrids; ++i) {
      grid_ptrs[i] = grids[i].ptr();
    }
    _ptr = hoppet_cxx__grid_def__new_from_grids(grid_ptrs.data(), ngrids, locked);
  }

};  

/// @brief  create a grid_def object with the default choice of nested, locked grids
/// @param dy      spacing of the coarsest grid
/// @param ymax    maximum y value for the coarsest grid
/// @param order   usual interpolation order parameter
/// @return 
inline grid_def grid_def_default(double dy, double ymax, int order) {
  return grid_def(hoppet_cxx__grid_def__new_default(dy, ymax, order));
}


//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_quant Fortran type, non-owning
///
/// This version provides a "view" onto an existing Fortran grid_quant object,
/// without taking ownership of it.
///
/// Note hoppet doesn't have grid_quant objects by default, instead it just
/// uses `real(dp) :: gq(0:ny)` arrays. The c++ wrapper instead explicitly
/// goes via a fortran grid_quant object, which manages the data array internally.
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
  //template<typename F, IsDoubleFunction<F> = 0>
  void assign(const DoubleFnDouble auto & fn) {
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
  void assign(double val) {ensure_valid(); std::fill(data(), data()+size(), val);}

  /// @brief  ensure that the grid_quant_view is valid (i.e. associated with data), otherwise throw an exception
  void ensure_valid() const {
    if (!data()) {throw std::runtime_error("hoppet::grid_quant_view::ensure_valid: data pointer is null");}
  }

  /// @brief  assign from a function, similar to assign()
  //template<typename F, IsDoubleFunction<F> = 0>
  grid_quant_view & operator=(const DoubleFnDouble auto & fn) {assign(fn); return *this;}

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
  grid_quant(const grid_quant & other)      {copy(other);}
  grid_quant(const grid_quant_view & other) {copy(other);}
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
  //template<class F>
  grid_quant(const grid_def_view & grid, const DoubleFnDouble auto & fn) : grid_quant(grid) {
    assign(fn);
  }

  //template<typename F>
  grid_quant & operator=(const DoubleFnDouble auto & fn) {assign(fn); return *this;}
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
inline grid_quant operator*(const grid_def_view & grid, DoubleFnDouble auto && fn) {
  // create an output grid_quant
  grid_quant result(grid);
  result = fn;
  return result;
}

///@}

typedef grid_quant_view pdf_flav_view;
typedef grid_quant      pdf_flav;


//-----------------------------------------------------------------------------
struct gq2d_extras {
  grid_def_view grid;
  std::size_t   size_dim1;
  std::size_t   size_dim0;
  gq2d_extras() : grid(), size_dim1(0), size_dim0(0) {}
  gq2d_extras(const gq2d_extras & other) : grid(other.grid), size_dim1(other.size_dim1), size_dim0(other.size_dim0) {}
  gq2d_extras(const grid_def_view & grid, std::size_t size_dim0) : grid(grid), size_dim1(grid.ny() + 1), size_dim0(size_dim0) {}
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
class grid_quant_2d_view : public data_view<gq2d_extras> {
public:  
  typedef data_view<gq2d_extras> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  grid_quant_2d_view() {}
  grid_def_view grid() const {return extras().grid;}
  grid_quant_view operator[](std::size_t i) {
    grid_quant_view result(this->data() + i * extras().size_dim1, extras().size_dim1, extras().grid);
    return result;
  }
  grid_quant_view operator()(std::size_t i) {return (*this)[i];}
  double & operator()(std::size_t i, std::size_t j) {return this->data()[i * extras().size_dim1 + j];}
  const double & operator()(std::size_t i, std::size_t j) const {return this->data()[i * extras().size_dim1 + j];}

  void assign(const VoidFnDoubleDoubleDoubleptr auto & fn, double Q) {
  //void assign(void (*fn)(double,double,double*), double Q) {
    if (!data()) {
      throw std::runtime_error("grid_quant_2d_view::assign(fn): grid_quant_2d_view object not associated");
    }
    if (extras().size_dim0 < iflv_max+1) {
      throw std::runtime_error("grid_quant_2d_view::assign(fn): grid_quant_2d_view dim1_sz = " 
                  + std::to_string(extras().size_dim0) + " < iflv_max+1 = " + std::to_string(iflv_max+1));
    }
    std::vector<double> xvals = grid().x_values();
    std::vector<double> xpdf (extras().size_dim0);
    for (std::size_t iy=0; iy < grid().ny()+1; ++iy) {
      fn(xvals[iy], Q, xpdf.data());
      for (std::size_t i=0; i < extras().size_dim0; ++i) {
        this->data()[i*extras().size_dim1 + iy] = i <= iflv_max ? xpdf[i] : 0.0;
      }
    }
  }

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
class grid_quant_2d : public data_owner<grid_quant_2d_view, grid_quant_2d_f> {

public:

  typedef grid_quant_2d_view view_type;

  grid_quant_2d() {}
  grid_quant_2d(const grid_def_view & grid, std::size_t dim1_size) {
    alloc(gq2d_extras(grid, dim1_size));
  }
  void alloc_virtual(const gq2d_extras & extras_in) override {
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
  typedef obj_view<grid_conv_f> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  grid_conv_view() = default;
  //grid_conv_view(const grid_def_view & grid) :  base_type(nullptr, grid) {}
  grid_conv_view(grid_conv_f * ptr) : base_type(ptr) {}

  const grid_def_view grid() const {return grid_def_view(hoppet_cxx__grid_conv__grid(ptr()));}

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
};

//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the grid_conv Fortran type, with ownership
class grid_conv : public obj_owner<grid_conv_view> {
public:

  typedef obj_owner<grid_conv_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  grid_conv() {}

  /// construct a grid_conv object and initialise it with the given function
  ///
  /// @param grid          the grid definition
  /// @param conv_ignd_fn  the convolution integrand function, double(double y, int piece)
  ///
  grid_conv(const grid_def_view & grid, DoubleFnDoubleInt auto && conv_ignd_fn) : base_type(nullptr) {

    //std::cout << "grid_conv: constructing from function object, grid = " << grid.ptr() <<" " << this->grid().ptr() << "\n";
    using FuncType = decltype(conv_ignd_fn);
    // under the hood, the fortran calls hoppet_grid_conv_f__wrapper
    // (defined in hoppet_oo.cc), with a pointer to the function object
    std::function<double(double,int)> fn_ptr = std::forward<FuncType>(conv_ignd_fn);
    _ptr = hoppet_cxx__grid_conv__new_from_fn(grid.ptr(), &fn_ptr);
  }
};

inline grid_quant operator*(const grid_conv_view & conv, const grid_quant_view & q) {
  // create an output grid_quant
  conv.grid().ensure_compatible(q.grid());
  grid_quant result(q.grid());

  hoppet_cxx__grid_conv__times_grid_quant(conv.ptr(), q.data(), result.data());
  return result;
}

inline grid_conv operator*(const grid_def_view & grid, DoubleFnDoubleInt auto && conv_ignd_fn) {
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

  typedef obj_view<split_mat_f> base_type;
  typedef grid_def_view extra_type;
  using base_type::base_type; // ensures that constructors are inherited

  split_mat_view & operator=(const split_mat_view & other) {
    if (!ptr()) throw std::runtime_error("split_mat_view::operator=: split_mat_view object not associated");
    hoppet_cxx__split_mat__copy_contents(ptr(), other.ptr());
    return *this;
  }

  const grid_def_view grid() const { return qq().grid(); }
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

  typedef obj_owner<split_mat_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited

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

  typedef obj_view<mass_threshold_mat_f> base_type;
  typedef grid_def_view extra_type;
  using base_type::base_type; // ensures that constructors are inherited

  mass_threshold_mat_view & operator=(const mass_threshold_mat_view & other) {
    if (!ptr()) throw std::runtime_error("mass_threshold_mat_view::operator=: mass_threshold_mat_view object not associated");
    hoppet_cxx__mass_threshold_mat__copy_contents(ptr(), other.ptr());
    return *this;
  }

  const grid_def_view grid() const { return pshq().grid(); }
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

  typedef obj_owner<mass_threshold_mat_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited

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
  typedef obj_view<running_coupling_f> base_type;
  using base_type::base_type; // ensures that constructors are inherited

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
  /// bypass fast interpolation and might need solve the evolution beyond the 
  /// cached range
  double operator()(double Q, int fixnf) const {return hoppet_cxx__running_coupling__value(valid_ptr(), Q, &fixnf);}

  /// @brief  Return the number of loops (nloop) for this running coupling
  /// @return nloop
  int nloop() const { return hoppet_cxx__running_coupling__num_loops(valid_ptr()); }

  /// @brief  Returns the range of active flavours (nflcl) for this running coupling
  /// @return  a tuple (nflcl_lo, nflcl_hi)
  ///
  /// use this as `auto [nflcl_lo, nflcl_hi] = rc.nf_range();`
  /// or (if nflcl_lo, nflcl_hi are already defined) `std::tie(nflcl_lo, nflcl_hi) = rc.nf_range();`
  std::tuple<int,int> nf_range() const {
    int lo=0, hi=0;
    hoppet_cxx__running_coupling__nf_range(valid_ptr(), &lo, &hi);
    return {lo, hi};
  }

  /// @brief  Returns the number of active flavours at a given scale Q  
  /// @param Q 
  /// @return 
  int nf_at_Q(double Q) const {
    return hoppet_cxx__running_coupling__nf_at_q(valid_ptr(), Q);
  }

  /// @brief  Returns the number of active flavours at a given scale Q, 
  ///         along with the range of Q for which this number of flavours is valid
  ///
  /// @param Q 
  /// @param Qlo is set to the lower edge of the Q range
  /// @param Qhi is set to the upper edge of the Q range
  /// @return number of active flavours at scale Q
  int nf_at_Q(double Q, double &Qlo, double & Qhi) const {
    return hoppet_cxx__running_coupling__nf_at_q(valid_ptr(), Q, &Qlo, &Qhi);
  }

  /// @brief  Returns the number of active flavours at a given scale Q, 
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
  typedef obj_owner<running_coupling_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  running_coupling() {}

  running_coupling(double alphas_at_Q, double Q, int nloop, int fixnf) {
    _ptr = hoppet_cxx__running_coupling__new_fixnf(alphas_at_Q, Q, nloop, fixnf);
  }

  running_coupling(double alphas_at_Q, double Q, int nloop,
                   double mc, double mb, double mt,
                   bool masses_are_MSbar = false, double muMatch_mQuark = 1.0) {
    _ptr = hoppet_cxx__running_coupling__new_varnf(alphas_at_Q, Q, nloop,
                                                  mc, mb, mt, masses_are_MSbar, muMatch_mQuark);
  }

};



//-----------------------------------------------------------------------------
/// @brief Object-oriented wrapper around the dglap_holder Fortran type, non-owning
class dglap_holder_view : public obj_view<dglap_holder_f> {
public:
  typedef obj_view<dglap_holder_f> base_type;
  using base_type::base_type; // ensures that constructors are inherited

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

  grid_def_view grid() const {return p_lo().grid();}  

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
  typedef obj_owner<dglap_holder_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  dglap_holder() {}
  dglap_holder(const hoppet::grid_def_view & grid, int factscheme, int nloop, int nflo, int nfhi) {
    _ptr = hoppet_cxx__dglap_holder__new(grid.ptr(), &factscheme, &nloop, &nflo, &nfhi);
  }
};

class pdfseginfo_view : public obj_view<pdfseginfo_f> {
public:
  typedef obj_view<pdfseginfo_f> base_type;
  using base_type::base_type; // ensures that constructors are inherited
  RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_lo)
  RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_hi)
  RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_lo)
  RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_hi)
  RETURN_DBL_MEMBER(pdfseginfo,dlnlnQ)
  RETURN_DBL_MEMBER(pdfseginfo,inv_dlnlnQ)
};

class pdf_table_view : public obj_view<pdf_table_f> {
public:
  typedef obj_view<pdf_table_f> base_type;
  using base_type::base_type; // ensures that constructors are inherited

  pdf_table_view & operator=(const pdf_table_view & other) {
    if (!ptr()) throw std::runtime_error("pdf_table_view::operator=: pdf_table_view object not associated");
    hoppet_cxx__pdf_table__copy_contents(ptr(), other.ptr());
    return *this;
  }

  grid_quant_view at_iQf(size_t iQ, size_t iflv) const {
    double * tab_ptr = hoppet_cxx__pdf_table__tab_ptr(valid_ptr());
    // these are the sizes of the two dimensions of result (shifted by 1 index wrt table)
    size_t size_dim0 = size_flv();
    size_t size_dim1 = static_cast<size_t>(grid().ny()+1);
    size_t iQ_size = size_dim0 * size_dim1;
    double * iQf_ptr  = tab_ptr + iQ * iQ_size + iflv * size_dim1;
    return grid_quant_view(iQf_ptr, size_dim1, grid());
  }

  grid_quant_2d_view at_iQ(size_t iQ) const {
    double * tab_ptr = hoppet_cxx__pdf_table__tab_ptr(valid_ptr());
    // these are the sizes of the two dimensions of result (shifted by 1 index wrt table)
    size_t size_dim0 = size_flv();
    size_t size_dim1 = static_cast<size_t>(grid().ny()+1);
    size_t iQ_size = size_dim0 * size_dim1;
    double * iQ_ptr  = tab_ptr + iQ * iQ_size;
    return grid_quant_2d_view(iQ_ptr, iQ_size, gq2d_extras(grid(), size_dim0));
  }

  /// return the maximum valid iflv index, in C++ numbering; 
  size_t iflv_max() const {return static_cast<size_t>(tab_iflv_max() - iflv_min_fortran);  }
  /// @brief return the size of the flavour (dim1=2nd) dimension of the table 
  size_t size_flv() const {return static_cast<size_t>(hoppet_cxx__pdf_table__size_flv(valid_ptr())); }

  // think carefully which of these interfaces should be public, which perhaps renamed
  RETURN_OBJ_MEMBER(pdf_table,grid,grid_def)
  RETURN_INT_MEMBER(pdf_table,nQ)
  RETURN_INT_MEMBER(pdf_table,lnlnQ_order)
  RETURN_LOG_MEMBER(pdf_table,freeze_at_Qmin)
  RETURN_LOG_MEMBER(pdf_table,nf_info_associated)
  RETURN_INT_MEMBER(pdf_table,nflo)
  RETURN_INT_MEMBER(pdf_table,nfhi)
  RETURN_DBL_MEMBER(pdf_table,lnlnQ_min)
  RETURN_DBL_MEMBER(pdf_table,lnlnQ_max)
  RETURN_DBL_MEMBER(pdf_table,lambda_eff)
  RETURN_DBL_MEMBER(pdf_table,dlnlnQ)
  RETURN_OBJ_MEMBER(pdf_table,seginfo_no_nf,pdfseginfo)
  RETURN_OBJ_MEMBER_I(pdf_table,seginfo,pdfseginfo)
  RETURN_DBL_MEMBER_I(pdf_table,as2pi)
  RETURN_INT_MEMBER_I(pdf_table,nf_int)
  RETURN_DBL_MEMBER_I(pdf_table,lnlnQ_vals)
  RETURN_DBL_MEMBER_I(pdf_table,Q_vals)


protected:
  /// @brief return the maximum flavour number for which the table has been allocate
  /// This is in Fortran numbering, which is why it is protected; use iflv_max() instead
  RETURN_INT_MEMBER(pdf_table,tab_iflv_max)

//  RETURN_OBJ_MEMBER(pdf_table,grid,grid_def)
//  RETURN_INT_MEMBER(pdf_table,nq)
//  RETURN_INT_MEMBER(pdf_table,tab_iflv_max)
//  RETURN_INT_MEMBER(pdf_table,lnlnq_order)
//  RETURN_LOG_MEMBER(pdf_table,freeze_at_qmin)
//  RETURN_LOG_MEMBER(pdf_table,nf_info_associated)
//  RETURN_INT_MEMBER(pdf_table,nflo)
//  RETURN_INT_MEMBER(pdf_table,nfhi)
//  RETURN_DBL_MEMBER(pdf_table,lnlnq_min)
//  RETURN_DBL_MEMBER(pdf_table,lnlnq_max)
//  RETURN_DBL_MEMBER(pdf_table,lambda_eff)
//  RETURN_DBL_MEMBER(pdf_table,dlnlnq)
//  RETURN_OBJ_MEMBER(pdf_table,seginfo_no_nf,pdfseginfo)
//  RETURN_OBJ_MEMBER_I(pdf_table,seginfo,pdfseginfo)
//  RETURN_DBL_MEMBER_I(pdf_table,as2pi)
//  RETURN_INT_MEMBER_I(pdf_table,nf_int)
//  RETURN_INT_MEMBER_I(pdf_table,lnlnq_vals)
//  RETURN_INT_MEMBER_I(pdf_table,q_vals)

};

class pdf_table : public obj_owner<pdf_table_view> {
public:
  typedef obj_owner<pdf_table_view> base_type;
  using base_type::base_type; // ensures that constructors are inherited  
  pdf_table() {}
  pdf_table(const grid_def_view & grid, double Qmin, double Qmax, 
            double dlnlnQ, int lnlnQ_order = hoppet_pdf_table_def_lnlnQ_order, 
            bool freeze_at_Qmin = false, 
            int iflv_max_table = ncompmax) {
    _ptr = hoppet_cxx__pdf_table__new(grid.ptr(), Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table);
  }

  void add_nf_info(const running_coupling & coupling) {
    hoppet_cxx__pdf_table__add_nf_info(valid_ptr(), coupling.ptr());
  }
};

} // end namespace hoppet -------------------------------------



/// objects globally defined in the streamlined interface
namespace hoppet {
namespace sl {
  /// a view of the grid_def object being used in the streamlined interface
  extern grid_def_view grid;
  extern dglap_holder_view dh;
  extern pdf_table_view    table; //< this should evolve to become an array of tables
  extern running_coupling_view coupling; 
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