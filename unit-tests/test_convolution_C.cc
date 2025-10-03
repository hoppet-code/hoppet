#include <iostream>
using namespace std;

extern "C" {
  double pqq_c(const double & x) {
    //cout << "pqq_c called with x = " << x << endl;
    return (1.0 + x*x)/(1.0 - x);
  }
}