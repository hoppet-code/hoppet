
extern "C" {
  extern double hoppet_global_qcd_ca;
  extern double hoppet_global_qcd_cf;
  extern double hoppet_global_qcd_tr;

  extern int    hoppet_global_qcd_nf_int;
  extern int    hoppet_global_qcd_nf_u  ;
  extern int    hoppet_global_qcd_nf_d  ;
  extern double hoppet_global_qcd_nf ;
  extern double hoppet_global_qcd_tf ;


  void hoppet_qcd_set_nf(const int & nf_in);
}

namespace hoppet {
namespace qcd {
  const double & ca     = hoppet_global_qcd_ca    ;
  const double & cf     = hoppet_global_qcd_cf    ;
  const double & tr     = hoppet_global_qcd_tr    ;
  const int    & nf_int = hoppet_global_qcd_nf_int;
  const int    & nf_u   = hoppet_global_qcd_nf_u  ;
  const int    & nf_d   = hoppet_global_qcd_nf_d  ;
  const double & nf     = hoppet_global_qcd_nf    ;
  const double & tf     = hoppet_global_qcd_tf    ;

  inline void set_nf(int nf_in) {hoppet_qcd_set_nf(nf_in);}
}
}