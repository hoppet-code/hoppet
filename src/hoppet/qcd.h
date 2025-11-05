
extern "C" {
  extern double hoppet__qcd__ca;
  extern double hoppet__qcd__cf;
  extern double hoppet__qcd__tr;

  extern int    hoppet__qcd__nf_int;
  extern int    hoppet__qcd__nf_u  ;
  extern int    hoppet__qcd__nf_d  ;
  extern double hoppet__qcd__nf ;
  extern double hoppet__qcd__tf ;


  void hoppet__qcd__set_nf(const int & nf_in);
  void hoppet__qcd__set_group(const double & ca_in, const double & cf_in, const double & tr_in);
}

namespace hoppet {
namespace qcd {
  const double & ca     = hoppet__qcd__ca    ;
  const double & cf     = hoppet__qcd__cf    ;
  const double & tr     = hoppet__qcd__tr    ;
  const int    & nf_int = hoppet__qcd__nf_int;
  const int    & nf_u   = hoppet__qcd__nf_u  ;
  const int    & nf_d   = hoppet__qcd__nf_d  ;
  const double & nf     = hoppet__qcd__nf    ;
  const double & tf     = hoppet__qcd__tf    ;

  auto & set_nf    = hoppet__qcd__set_nf;
  auto & set_group = hoppet__qcd__set_group;
}
}