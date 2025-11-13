
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

  const double & ca     = hoppet__qcd__ca    ;  ///< $C_A$ colour factor
  const double & cf     = hoppet__qcd__cf    ;  ///< $C_F$ colour factor
  const double & tr     = hoppet__qcd__tr    ;  ///< $T_R$ colour factor
  const int    & nf_int = hoppet__qcd__nf_int;  ///< number of active flavours as an integer
  const int    & nf_u   = hoppet__qcd__nf_u  ;  ///< number of active up-type flavours
  const int    & nf_d   = hoppet__qcd__nf_d  ;  ///< number of active down-type flavours
  const double & nf     = hoppet__qcd__nf    ;  ///< number of active flavours as a double
  const double & tf     = hoppet__qcd__tf    ;  ///< $T_F = n_f T_R$ factor

  const double & beta0  = hoppet__qcd__beta0;
  const double & beta1  = hoppet__qcd__beta1;
  const double & beta2  = hoppet__qcd__beta2;
  const double & beta3  = hoppet__qcd__beta3;

  const double & twopi_beta0 = hoppet__qcd__twopi_beta0;

  const double & alphastep11 = hoppet__qcd__alphastep11;
  const double & alphastep22 = hoppet__qcd__alphastep22;
  const double & alphastep21 = hoppet__qcd__alphastep21;
  const double & alphastep20_msbar = hoppet__qcd__alphastep20_msbar;
  const double & alphastep20_pole  = hoppet__qcd__alphastep20_pole;
  const double & alphastep33   = hoppet__qcd__alphastep33;
  const double & alphastep32   = hoppet__qcd__alphastep32;
  const double & alphastep32_inv = hoppet__qcd__alphastep32_inv;
  const double & alphastep31_msbar = hoppet__qcd__alphastep31_msbar;
  const double & alphastep31_msbar_inv = hoppet__qcd__alphastep31_msbar_inv;
  const double & alphastep31_msbar_nl  = hoppet__qcd__alphastep31_msbar_nl;
  const double & alphastep30_msbar = hoppet__qcd__alphastep30_msbar;
  const double & alphastep30_msbar_nl = hoppet__qcd__alphastep30_msbar_nl;
  const double & alphastep31_pole = hoppet__qcd__alphastep31_pole;
  const double & alphastep31_pole_inv = hoppet__qcd__alphastep31_pole_inv;
  const double & alphastep31_pole_nl = hoppet__qcd__alphastep31_pole_nl;
  const double & alphastep30_pole = hoppet__qcd__alphastep30_pole;
  const double & alphastep30_pole_nl = hoppet__qcd__alphastep30_pole_nl;

  const double & cmw_K  = hoppet__qcd__cmw_K;
  const double & cmw_K2 = hoppet__qcd__cmw_K2;
  const double & mvv_A3 = hoppet__qcd__mvv_A3;
  const double & mvv_A3G= hoppet__qcd__mvv_A3G;

  /// @brief  set the number of active flavours (nf) globally
  /// @param nf_in  the new number of active flavours (integer)
  inline void set_nf(const int & nf_in) {hoppet__qcd__set_nf(nf_in);}

  /// @brief  set the QCD group constants globally
  /// @param ca_in  the new value of C_A
  /// @param cf_in  the new value of C_F
  /// @param tr_in  the new value of T_R
  auto & set_group = hoppet__qcd__set_group;
}
}