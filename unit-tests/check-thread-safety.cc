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
  double Qval(int i) const { 
    return Qmin * exp(i * log(Qmax/Qmin) / nxQ); 
  }
};

/// class for a multi-flavour PDF evaluation task
class PDFThreadTask : public TestBase {
public:
  PDFThreadTask() {}
  void run() {
    for (unsigned i = 0; i <= nxQ; ++i) {
      double x = xmin * exp(i * log(xmax/xmin) / nxQ);
      double Q = Qval(i);
      double pdf[13];
      hoppetEval(x, Q, pdf);
      results.push_back(pdf[6+2] - pdf[6-2]); // u - ubar
      results.push_back(pdf[6+1] - pdf[6-1]); // d - dbar
      results.push_back(2 * (pdf[6-1] + pdf[6-2])); // 2(ubar + dbar)
      results.push_back(pdf[6+4] + pdf[6-4]); // c + cbar
      results.push_back(pdf[6+0]); // gluon
    }
  }
  string name() const { return "PDFThreadTask"; }
};

/// class for a multi-flavour PDF evaluation task
class PDFOneFlavThreadTask : public TestBase {
public:
  PDFOneFlavThreadTask() {}
  void run() {
    for (unsigned i = 0; i <= nxQ; ++i) {
      double x = xmin * exp(i * log(xmax/xmin) / nxQ);
      double Q = Qval(i);
      double pdf[13];
      results.push_back(hoppetEvalPID(x, Q, i%5));
    }
  }
  string name() const { return "PDFOneFlavThreadTask"; }
};

/// class for a AlphaS evaluation task
class AlphaSThreadTask : public TestBase {
public:
  AlphaSThreadTask() {}
  void run() {
    for (unsigned i = 0; i <= nxQ; ++i) {
      double Q = Qval(i);
      results.push_back(hoppetAlphaS(Q)); // alpha_s
    }
  }
  string name() const { return "AlphaSThreadTask"; }
};


/// Code to run the thread safety check for a class of type T
/// It calls T::run() in multiple threads and checks that the T::results
/// vector is identical across all the results
template<class T>
bool check_thread_safety(int nrep = 20) {
  cout << "Checking thread safety of " << T().name() << " with " << nrep << " repetitions..." << endl;
  bool fail = false;
  int nfail = 0, nfail_max = 5;
  for (int irep = 0; irep < nrep; ++irep) {
    constexpr int nthreads = 8;

    // the tasks are defined above in PDFThreadTask
    T tasks[nthreads];

    // create the threads, each of which runs a task
    vector<thread> threads;
    for (int i = 0; i < nthreads; ++i) {
      threads.emplace_back(&T::run, &tasks[i]);
    }

    // wait for all threads to finish
    for (auto& t : threads) {
      t.join();
    }

    // check the results are the same across all tasks
    for (int i = 1; i < nthreads; ++i) {
      //if (tasks[i].results != tasks[0].results) {
      //  cout << red 
      //       << " ↳ mismatch in thread " << i 
      //       << ", repetition " << irep;
      if (tasks[i].results.size() != tasks[0].results.size()) {
        cout << red << " ↳ Failure for repetition " << irep << ", thread " << i << ": different results sizes "
              << tasks[i].results.size() << " vs " 
              << tasks[0].results.size() << reset << endl;
        fail = true;
        nfail++;
      } else {
        //cout << ": same size " << tasks[i].results.size() << ", first few mismatches:" << endl;
        for (size_t j = 0; j < tasks[i].results.size(); ++j) {
          if (tasks[i].results[j] != tasks[0].results[j]) {
            fail = true;
            cout << red << " ↳ Failure for repetition " << irep 
                  << "   index " << j 
                  << ", thread 0: " << tasks[0].results[j]
                  << ", thread " << i << ": " << tasks[i].results[j]
                  << ", diff: " << fabs(tasks[0].results[j] - tasks[i].results[j])
                  << reset << endl;
            nfail++;
            if (nfail >= nfail_max) break;
          }
        }
      }
      if (nfail >= nfail_max) break;
    }
    if (nfail >= nfail_max) {
      cout << red << " ↳ ... too many failures, stopping output." << reset << endl;
      break;
    }
  }
  if (!fail) cout << green << " ↳ all threads produced identical results." << reset << endl;
  return !fail;
}

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

  int nrep = 0;
  bool global_pass = false;


  global_pass |= check_thread_safety<PDFThreadTask>();
  global_pass |= check_thread_safety<PDFOneFlavThreadTask>();
  global_pass |= check_thread_safety<AlphaSThreadTask>();

  hoppetSetYLnlnQInterpOrders(4,3); // to try a case with warnings
  global_pass |= check_thread_safety<PDFThreadTask>();

  if (!global_pass) {
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

