// code for testing thread safety of hoppet PDF and coupling evaluations.
//
// It runs multiple threads that all evaluate the PDF at a range of
// (x,Q) points and the coupling at a range of Q points and checks that
// the results are identical across threads



#include "hoppet.h"
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

constexpr char red[]{"\033[0;31m"};
constexpr char green[]{"\033[0;32m"};
constexpr char reset[]{"\033[0m"};

// definition of the initial condition function
void heralhc_init(const double & x,
                   const double & Q,
                   double * pdf);

/// base class with limits and number of x,Q points
class TestBase {
public:
  constexpr static double xmin = 1e-5, xmax = 0.99;
  constexpr static double Qmin = 1.5, Qmax = 10000.0;
  constexpr static int nxQ = 1000;
  vector<double> results;
};

/// class for a PDF (and coupling) evaluation task
class PDFThreadTask : public TestBase {
public:
  PDFThreadTask() {}
  void run() {
    for (unsigned i = 0; i <= nxQ; ++i) {
      double x = xmin * exp(i * log(xmax/xmin) / nxQ);
      double Q = Qmin * exp(i * log(Qmax/Qmin) / nxQ);
      double pdf[13];
      hoppetEval(x, Q, pdf);
      results.push_back(pdf[6+2] - pdf[6-2]); // u - ubar
      results.push_back(pdf[6+1] - pdf[6-1]); // d - dbar
      results.push_back(2 * (pdf[6-1] + pdf[6-2])); // 2(ubar + dbar)
      results.push_back(pdf[6+4] + pdf[6-4]); // c + cbar
      results.push_back(pdf[6+0]); // gluon
      results.push_back(hoppetAlphaS(Q)); // alpha_s
    }
  }

};                  

//----------------------------------------------------------------------
int main(int argc, char** argv) {

  double dy    = 0.1;
  int    nloop = 3;

  // initialise with NNLO, VFN
  hoppetStart(dy, nloop);
  // hoppetSetPoleMassVFN(1.414213563, 4.5, 175.0);
  
  // evolve the initial condition
  double asQ0 = 0.35, Q0=sqrt(2.0);
  hoppetEvolve(asQ0, Q0, nloop, 1.0, heralhc_init, Q0);

  int nrep = 20;
  bool global_fail = false;
  for (int irep = 0; irep < nrep; ++irep) {
    constexpr int nthreads = 8;

    // the tasks are defined above in PDFThreadTask
    PDFThreadTask tasks[nthreads];

    // create the threads, each of which runs a task
    vector<thread> threads;
    for (int i = 0; i < nthreads; ++i) {
      threads.emplace_back(&PDFThreadTask::run, &tasks[i]);
    }

    // wait for all threads to finish
    for (auto& t : threads) {
      t.join();
    }

    // check the results are the same across all tasks
    bool local_fail = false;
    for (int i = 1; i < nthreads; ++i) {
      if (tasks[i].results != tasks[0].results) {
        cout << red << "irep=" << irep << ": mismatch in thread " << i << reset << endl;
        local_fail = true;
        break;
      }
    }
    if (!local_fail) cout << green << "irep=" << irep << ": all threads produced identical results." << endl;
    global_fail = global_fail || local_fail;
  }

  if (global_fail) {
    cout << red << "Thread safety test failed!" << reset << endl;
    return 1;
  } else {
    cout << green << "Thread safety test passed!" << reset << endl;
    return 0;
  }
}


//----------------------------------------------------------------------
// the initial condition
void  heralhc_init(const double & x,
                   const double & Q,
                   double * pdf) {
  double uv, dv;
  double ubar, dbar;
  double N_g=1.7, N_ls=0.387975;
  double N_uv=5.107200, N_dv=3.064320;
  double N_db=N_ls/2;

  uv = N_uv * pow(x,0.8) * pow((1-x),3);
  dv = N_dv * pow(x,0.8) * pow((1-x),4);
  dbar = N_db * pow(x,-0.1) * pow(1-x,6);
  ubar = dbar * (1-x);

  pdf[ 0+6] = N_g * pow(x,-0.1) * pow(1-x,5);
  pdf[-3+6] = 0.2*(dbar + ubar);
  pdf[ 3+6] = pdf[-3+6];
  pdf[ 2+6] = uv + ubar;
  pdf[-2+6] = ubar;
  pdf[ 1+6] = dv + dbar;
  pdf[-1+6] = dbar;

  pdf[ 4+6] = 0;
  pdf[ 5+6] = 0;
  pdf[ 6+6] = 0;
  pdf[-4+6] = 0;
  pdf[-5+6] = 0;
  pdf[-6+6] = 0;
}

