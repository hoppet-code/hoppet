#include "apfel/splittingfunctionsunp_tl.h"
#include "apfel/splittingfunctionsunp_sl.h"
#include <iostream>
#include <vector>

using namespace std;

double val(const apfel::Expression & split, double x) {
  return split.Regular(x)+split.Singular(x);
}

int i, n=100;
int nf=5;
vector<double> xvals;

void test_qg() {
  auto split_obj = apfel::P1Tqg(nf);
  for (double x: xvals)  {
    cout << x << " " << 0.25*(val(split_obj,x)) << endl;
  }
}
void test_gq() {
  auto split_obj = apfel::P1Tgq(nf);
  for (double x: xvals)  {
    cout << x << " " << 0.25*(val(split_obj,x)) << endl;
  }
}
void test_qq() {
  auto nsp = apfel::P1Tnsp(nf);
  auto ps =  apfel::P1Tps(nf);
  for (double x: xvals)  {
    cout << x << " " << 0.25*(val(nsp,x) + val(ps,x)) << endl;
  }
}

int main(int argc, char ** argv) {
  cout.precision(15);
  xvals.resize(n-1);
  for (i = 1; i < n; i++) {
    xvals[i-1] = (i-0.9)*1.0/n;
  }
  //test_qq();
  test_gq();
}
