#include "apfel/splittingfunctionsunp_tl.h"
#include "apfel/splittingfunctionsunp_sl.h"
#include <iostream>

using namespace std;

int main(int argc, char ** argv) {
  int i, n=100;
  int nf = 5;
  auto split_obj = apfel::P1Tqg(nf);
  cout.precision(15);
  for (i = 1; i < n; i++) {
    double x = (i-0.9)*1.0/n;
    cout << x << " " << 0.25*(split_obj.Regular(x)+split_obj.Singular(x)) << endl;
  }
}
