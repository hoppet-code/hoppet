
extern "C" {
  extern double hoppet_module_qcd__ca;
  extern double hoppet_module_qcd__cf;
  extern double hoppet_module_qcd__tr;

  extern int    hoppet_module_qcd__nf_int;
  extern int    hoppet_module_qcd__nf_u  ;
  extern int    hoppet_module_qcd__nf_d  ;
  extern double hoppet_module_qcd__nf ;
  extern double hoppet_module_qcd__tf ;


  void hoppet_module_qcd__set_nf(const int & nf_in);
  void hoppet_module_qcd__set_group(const double & ca_in, const double & cf_in, const double & tr_in);
}

namespace hoppet {
namespace qcd {
  const double & ca     = hoppet_module_qcd__ca    ;
  const double & cf     = hoppet_module_qcd__cf    ;
  const double & tr     = hoppet_module_qcd__tr    ;
  const int    & nf_int = hoppet_module_qcd__nf_int;
  const int    & nf_u   = hoppet_module_qcd__nf_u  ;
  const int    & nf_d   = hoppet_module_qcd__nf_d  ;
  const double & nf     = hoppet_module_qcd__nf    ;
  const double & tf     = hoppet_module_qcd__tf    ;

  auto & set_nf    = hoppet_module_qcd__set_nf;
  auto & set_group = hoppet_module_qcd__set_group;
}
}