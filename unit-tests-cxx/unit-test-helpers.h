#ifndef UNIT_TEST_HELPERS_H
#define UNIT_TEST_HELPERS_H

#include "hoppet_oo.h"

constexpr double cf = 4.0/3.0;
inline double pqq_fn(double y, int piece) {
  double x = exp(-y);
  double offset = -cf; // the result should be independent of the offset value

//  double res = 0.0;
//  if      (piece == hoppet::cc_REAL    ) res =   cf * (1+x*x)/(1-x);
//  else if (piece == hoppet::cc_VIRT    ) res = - cf * (1+x*x)/(1-x) + offset*(1+x);
//  else if (piece == hoppet::cc_REALVIRT) res = offset*(1+x);
//  else if (piece == hoppet::cc_DELTA   ) res = -offset*3.0/2.0;
//
//  if (piece != hoppet::cc_DELTA) res *= x;
//  //return res;
//  
//  double res1 = 0.0;
//  if (piece == hoppet::cc_REAL)     res1   = cf * (1.0 + x*x)/ (1.0-x);
//  if (piece == hoppet::cc_VIRT)     res1   = -2*cf/(1.0-x);
//  if (piece == hoppet::cc_REALVIRT) res1   = -cf*(1.0+x);
//  if (piece == hoppet::cc_DELTA)    res1 = cf * (3.0/2.0);
//  if (piece != hoppet::cc_DELTA) res1 *= x;


//  if (abs(res-res1) > 1e-8) cout << "pqq_fn: piece=" << piece << ", res=" << res << ", res1=" << res1 << ", x=" << x << endl;
  //return res1;


  double res2 = 0.0;
  if (piece == hoppet::cc_REAL
   || piece == hoppet::cc_REALVIRT) res2 = cf * (1.0 + x*x)/ (1.0-x);

  if (piece == hoppet::cc_VIRT 
   || piece == hoppet::cc_REALVIRT) res2 -= cf * 2.0/(1.0-x);

  if (piece == hoppet::cc_DELTA) res2 = cf * (3.0/2.0);

  if (piece != hoppet::cc_DELTA) res2 *= x;
//  if (abs(res-res2) > 1e-8) cout << "pqq_fn: piece=" << piece << ", res=" << res << ", res2=" << res2 << ", x=" << x << endl;

  return res2;
}

inline double pgq_fn(double y, int piece) {
  double x = exp(-y);
  double res;
  if      (piece == hoppet::cc_REAL
        || piece == hoppet::cc_REALVIRT) res = cf * (1+pow(1-x,2))/x;
  else                               res = 0.0;
  return x*res;
}

inline double pgg_fn(double y, int piece, double nf) {
  double x = exp(-y);
  constexpr double ca = 3.0;
  double res = 0.0;
  if (piece == hoppet::cc_REAL
   || piece == hoppet::cc_REALVIRT) res = 2.0*ca*(x/(1.0-x) + (1.0-x)/x + x*(1.0-x));

  if (piece == hoppet::cc_VIRT 
   || piece == hoppet::cc_REALVIRT) res -= 2.0*ca/(1.0-x);

  if (piece == hoppet::cc_DELTA) res = (11.0*ca - 2.0*nf)/6.0;

  if (piece != hoppet::cc_DELTA) res *= x;
  return res;

}

inline double pqg_fn(double y, int piece, double nf) {
  double x = exp(-y);
  double res = 0.0;
  if      (piece == hoppet::cc_REAL
        || piece == hoppet::cc_REALVIRT) res = nf*(x*x + (1-x)*(1-x));
  return x*res;
}

#endif // UNIT_TEST_HELPERS_H