#include <iostream>
#include <cmath>
#include "hoppet.h"
using namespace std;

extern "C" {
  double pqq_c(const double x) {
    double result = (1.0 + x*x)/(1.0 - x);
    return result;
  }

  double xpqqy_c(const double y) {

    double x = exp(-y);
    double result = 0.0;
    constexpr double CF = 4.0/3.0; // CF for QCD
    switch (hoppet_global_cc_piece) {
    case hoppet::cc_REAL:
    case hoppet::cc_REALVIRT:
      result = CF * (1.0+x*x)/(1.0-x);
    }

    switch (hoppet_global_cc_piece) {
    case hoppet::cc_VIRT:
    case hoppet::cc_REALVIRT:
      result -= CF * 2.0 / (1.0 - x);
      break;
    case hoppet::cc_DELTA:
      cout << "In pqqy_c with cc_piece=DELTA" << endl;
      result = CF * 1.5;
    }

    if (hoppet_global_cc_piece != hoppet::cc_DELTA) result *= x;
    return result;
  }
}

//class A{};
//
//A a;
//A b;
//template<FA, FB> 
//class U {
//public:
//  U() 
//  A * FA, *FB;
//}
//
//U<a,a> u1;