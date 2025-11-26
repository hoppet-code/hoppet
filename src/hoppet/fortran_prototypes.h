/// @file hoppet/fortran_prototypes.h
/// @brief C++ prototypes for the Fortran routines used in the object-oriented interface
///
/// See hoppet_cxx_oo.f90 for the Fortran implementations of these routines.

#ifndef __HOPPET_FORTRAN_PROTOTYPES_H__
#define __HOPPET_FORTRAN_PROTOTYPES_H__

#include "hoppet.h"
#include "hoppet/base_types.h"

#ifdef __cplusplus
  #define _Bool bool
#endif

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


/// "forward" declarations of classes corresponding to Fortran types.
///
/// note that we only ever use pointers to to these class, so we do not
/// need to know anything about their actual structure, which is hidden in
/// the Fortran code.
///
/// @{
/// "forward" declaration of the Fortran grid_def type;
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
class pdf_table_array_f;
/// @}


/// grid_def function wrappers
extern "C" {
  extern const double hoppet_default_conv_eps;
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
DEFINE_RETURN_LOG_MEMBER(grid_def,locked)
DEFINE_RETURN_OBJ_MEMBER_I(grid_def,subgd,grid_def)
DEFINE_RETURN_INT_MEMBER_I(grid_def,subiy)
inline void generic_delete(grid_def_f * ptr) {if (ptr) hoppet_cxx__grid_def__delete(&ptr);}
inline grid_def_f * generic_copy(const grid_def_f * ptr) {  return hoppet_cxx__grid_def__copy(ptr);}

/// grid_quant function wrappers
extern "C" {
  grid_quant_f * hoppet_cxx__grid_quant__new(const grid_def_f * griddef);
  void   hoppet_cxx__grid_quant__delete(grid_quant_f ** gridquant);
  double hoppet_cxx__grid_quant__at_y(const grid_def_f * griddef, const double * gq_data, double y);
  double   hoppet_cxx__grid_quant__trunc_mom(const grid_def_f * griddef, const double * data, double n, const double * ymax = 0);
  double * hoppet_cxx__grid_quant__data_ptr(grid_quant_f * gridquant);

  void hoppet_cxx__grid_quant__luminosity(const grid_def_f * griddef, const double * gq1_data, const double * gq2_data, double * result_data);
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
  grid_conv_f * hoppet_cxx__grid_conv__new_from_fn(const grid_def_f * grid_ptr, void * conv_ignd_c_fn_obj, 
                                                   const double * split_array_ptr = nullptr, int * nsplit = nullptr);
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
  running_coupling_f * hoppet_cxx__running_coupling__new_fixnf(double alphaS, double Q, int nloop, int fixnf, const double * Qmax = nullptr);
  running_coupling_f * hoppet_cxx__running_coupling__new_varnf(double alphaS, double Q, int nloop,
                                                              double mc, double mb, double mt,
                                                              bool masses_are_MSbar, double muMatch_mQuark, const double * Qmax = nullptr);
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

DEFINE_RETURN_OBJ_MEMBER(dglap_holder,grid,grid_def)
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
  double   hoppet_cxx__pdf_table__at_yqf(const pdf_table_f * tab, double y, double Q, int iflv);
  void     hoppet_cxx__pdf_table__at_yq_into(const pdf_table_f * tab, double y, double Q, double * res);
  void     hoppet_cxx__pdf_table__evolve(pdf_table_f *tab, double q0, const double * pdf_at_q0,
                                         const dglap_holder_f *dh, const running_coupling_f *coupling, 
                                         double mur_over_q, int nloop, _Bool untie_nf);
  void     hoppet_cxx__pdf_table__pre_evolve(pdf_table_f *tab, double q0, 
                                         const dglap_holder_f *dh, const running_coupling_f *coupling, 
                                         double mur_over_q, int nloop, _Bool untie_nf);
  void     hoppet_cxx__pdf_table__evolve_frompre(pdf_table_f *tab, const double * pdf_at_q0);
  void hoppet_cxx__pdf_table__write_lhapdf(const pdf_table_f * tab, const running_coupling_f * coupling, 
                                           const char * filename_cstr, int pdf_index, const int * iy_increment = nullptr, 
                                           const int * n_flav = nullptr, const int * flav_indices = nullptr, 
                                           const int * flav_pdg_ids = nullptr, const double * flav_rescale = nullptr);
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

//----- things for pdf_table arrays ----------
extern "C" {
  void hoppet_cxx__pdf_table_array__new(int sz, pdf_table_array_f **ptr, pdf_table_f **table0_ptr);
  
  pdf_table_f * hoppet_cxx__pdf_tables__table_i(const pdf_table_f * arr_start, int sz, int i);
}
DEFINE_DELETE(pdf_table_array)
#endif // __HOPPET_FORTRAN_PROTOTYPES_H__