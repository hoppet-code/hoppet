#pragma once

extern "C" {
  extern double hoppet__qcd__ca;
  extern double hoppet__qcd__cf;
  extern double hoppet__qcd__tr;

  extern int    hoppet__qcd__nf_int;
  extern int    hoppet__qcd__nf_u  ;
  extern int    hoppet__qcd__nf_d  ;
  extern double hoppet__qcd__nf ;
  extern double hoppet__qcd__tf ;

  extern double hoppet__qcd__beta0;
  extern double hoppet__qcd__twopi_beta0;
  extern double hoppet__qcd__beta1;
  extern double hoppet__qcd__beta2;
  extern double hoppet__qcd__beta3;

  extern double hoppet__qcd__alphastep11;
  extern double hoppet__qcd__alphastep22;
  extern double hoppet__qcd__alphastep21;
  extern double hoppet__qcd__alphastep20_msbar;
  extern double hoppet__qcd__alphastep20_pole;
  extern double hoppet__qcd__alphastep33;
  extern double hoppet__qcd__alphastep32;
  extern double hoppet__qcd__alphastep32_inv;
  extern double hoppet__qcd__alphastep31_msbar;
  extern double hoppet__qcd__alphastep31_msbar_inv;
  extern double hoppet__qcd__alphastep31_msbar_nl;
  extern double hoppet__qcd__alphastep30_msbar;
  extern double hoppet__qcd__alphastep30_msbar_nl;
  extern double hoppet__qcd__alphastep31_pole;
  extern double hoppet__qcd__alphastep31_pole_inv;
  extern double hoppet__qcd__alphastep31_pole_nl;
  extern double hoppet__qcd__alphastep30_pole;
  extern double hoppet__qcd__alphastep30_pole_nl;

  extern double hoppet__qcd__cmw_K;
  extern double hoppet__qcd__cmw_K2;
  extern double hoppet__qcd__mvv_A3;
  extern double hoppet__qcd__mvv_A3G;


  void hoppet__qcd__set_nf(const int & nf_in);
  void hoppet__qcd__set_group(const double & ca_in, const double & cf_in, const double & tr_in);
}

namespace hoppet {
/// @brief namespace containing QCD-related global constants and functions  
namespace qcd {

  inline const double & ca     = hoppet__qcd__ca    ;  ///< \f$C_A\f$ colour factor
  inline const double & cf     = hoppet__qcd__cf    ;  ///< \f$C_F\f$ colour factor
  inline const double & tr     = hoppet__qcd__tr    ;  ///< \f$T_R\f$ colour factor
  inline const int    & nf_int = hoppet__qcd__nf_int;  ///< number of active flavours as an integer
  inline const int    & nf_u   = hoppet__qcd__nf_u  ;  ///< number of active up-type flavours
  inline const int    & nf_d   = hoppet__qcd__nf_d  ;  ///< number of active down-type flavours
  inline const double & nf     = hoppet__qcd__nf    ;  ///< number of active flavours as a double
  inline const double & tf     = hoppet__qcd__tf    ;  ///< \f$T_F = n_f T_R\f$ factor

  inline const double & beta0  = hoppet__qcd__beta0;
  inline const double & beta1  = hoppet__qcd__beta1;
  inline const double & beta2  = hoppet__qcd__beta2;
  inline const double & beta3  = hoppet__qcd__beta3;

  inline const double & twopi_beta0 = hoppet__qcd__twopi_beta0;

  inline const double & alphastep11 = hoppet__qcd__alphastep11;
  inline const double & alphastep22 = hoppet__qcd__alphastep22;
  inline const double & alphastep21 = hoppet__qcd__alphastep21;
  inline const double & alphastep20_msbar = hoppet__qcd__alphastep20_msbar;
  inline const double & alphastep20_pole  = hoppet__qcd__alphastep20_pole;
  inline const double & alphastep33   = hoppet__qcd__alphastep33;
  inline const double & alphastep32   = hoppet__qcd__alphastep32;
  inline const double & alphastep32_inv = hoppet__qcd__alphastep32_inv;
  inline const double & alphastep31_msbar = hoppet__qcd__alphastep31_msbar;
  inline const double & alphastep31_msbar_inv = hoppet__qcd__alphastep31_msbar_inv;
  inline const double & alphastep31_msbar_nl  = hoppet__qcd__alphastep31_msbar_nl;
  inline const double & alphastep30_msbar = hoppet__qcd__alphastep30_msbar;
  inline const double & alphastep30_msbar_nl = hoppet__qcd__alphastep30_msbar_nl;
  inline const double & alphastep31_pole = hoppet__qcd__alphastep31_pole;
  inline const double & alphastep31_pole_inv = hoppet__qcd__alphastep31_pole_inv;
  inline const double & alphastep31_pole_nl = hoppet__qcd__alphastep31_pole_nl;
  inline const double & alphastep30_pole = hoppet__qcd__alphastep30_pole;
  inline const double & alphastep30_pole_nl = hoppet__qcd__alphastep30_pole_nl;

  inline const double & cmw_K  = hoppet__qcd__cmw_K;
  inline const double & cmw_K2 = hoppet__qcd__cmw_K2;
  inline const double & mvv_A3 = hoppet__qcd__mvv_A3;
  inline const double & mvv_A3G= hoppet__qcd__mvv_A3G;

  /// @brief  set the number of active flavours (nf) globally
  ///
  /// @param nf_in  the new number of active flavours (integer)
  ///
  /// This function updates the values of all the QCD constants defined
  /// in this namespace, including e.g. beta-function coefficients
  inline void set_nf(const int & nf_in) {hoppet__qcd__set_nf(nf_in);}

  /// @brief  set the QCD group constants globally
  ///
  /// @param ca_in  the new value of \f$C_A\f$
  /// @param cf_in  the new value of \f$C_F\f$
  /// @param tr_in  the new value of \f$T_R\f$
  ///
  /// This function updates the values of all the QCD constants defined
  /// in this namespace, including e.g. beta-function coefficients
  inline void set_group(const double & ca_in, const double & cf_in, const double & tr_in) {
    hoppet__qcd__set_group(ca_in, cf_in, tr_in);
  }
}
}