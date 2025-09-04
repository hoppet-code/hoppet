/// An example in C++ showing how to replace LHAPDF interpolation with
/// (faster) hoppet interpolation
///
#include "../src/hoppet.h"
#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>

using namespace std;
using namespace hoppet; // To access factscheme_MSbar
using namespace LHAPDF;

// Global PDF pointer
PDF *pdf = nullptr;

/// Interface to LHAPDF as needed by hoppetAssign
void lhapdf_interface(const double & x, const double & Q, double * res)  {
  vector<double> data(13, 0.0); // Pre-allocate and zero;
  pdf->xfxQ2(x, Q*Q, data);
  copy(data.begin(), data.end(), res); // Fast copy to the output array
};

/// Routine that loads an LHAPDF set, extracts some information from it
/// and transfers the PDF to hoppet. It needs the lhapdf_interface
/// defined above.
void load_lhapdf_assign_hoppet(const string & pdfname, int imem=0) {
  // Start by loading a PDF set from LHAPDF
  pdf = mkPDF(pdfname, imem);

  // Now let's access some basic information about the PDF to set up hoppet
  double mc = pdf->quarkMass(4);
  double mb = pdf->quarkMass(5);
  double mt = pdf->quarkMass(6);
  int nloop = pdf->orderQCD() + 1; // LHAPDF is zero indexed
  double xmin = pdf->xMin();
  double xmax = pdf->xMax();
  double Qmin = pdf->qMin();
  double Qmax = pdf->qMax();

  cout << "LHAPDF set: " << pdfname << " loaded successfully" << endl;

  // Now let us define some hoppet specific parameters. These are
  // typical values, and should guarantee similar accuracy as can be
  // expected from LHAPDF
  double ymax = static_cast<float>(std::ceil(log(1.0/xmin))); // To get a nice value of ymax that can contain the full LHAPDF grid
  double dy = 0.05;
  double dlnlnQ = dy/4.0;
  if(ymax > 15.0){
    dlnlnQ = dy/8.0; // for large ymax we need a finer grid in Q
  }
  int order = -6; // Default
  int yorder = 2; // Quadratic interpolation in y
  int lnlnQorder = 2; // Quadratic interpolation in lnlnQ

  // Print all the relevant parameters that we have passed to hoppet
  cout << "Hoppet starting with:" << endl;
  cout << " ymax:       " << ymax << endl;
  cout << " dy:         " << dy << endl;
  cout << " Qmin:       " << Qmin << endl;
  cout << " Qmax:       " << Qmax << endl;
  cout << " dlnlnQ:     " << dlnlnQ << endl;
  cout << " nloop:      " << nloop << endl;
  cout << " order:      " << order << endl;
  cout << " yorder:     " << yorder << endl;
  cout << " lnlnQorder: " << lnlnQorder << endl;

  hoppetSetPoleMassVFN(mc,mb,mt); // set the pole masses
  hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder); // Set the interpolation orders
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme_MSbar); // Start hoppet

  // Now we fill the hoppet grid using the LHAPDF grid directly,
  // rather than evolving ourselves
  double Q0 = Qmin;
  hoppetAssignWithCoupling(lhapdf_interface, pdf->alphasQ(Q0), Q0, nloop);

  // If instead we want to evolve the PDF with hoppet starting from
  // some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
  // of hoppetAssign
  //hoppetEvolve(pdf->alphasQ(Q0), Q0, nloop, 1.0, lhapdf_interface, Q0);

}

//----------------------------------------------------------------------
int main () {
  string pdfname = "PDF4LHC21_40";
  int imem = 0;
  load_lhapdf_assign_hoppet(pdfname, imem);

  double x = 0.01;
  double Q = 13.0;
 
  // Standard LHAPDF call
  vector<double> lhapdf;
  pdf->xfxQ2(x, Q*Q, lhapdf); // Get the PDF from LHAPDF

  // Equivalent hoppet call
  double hoppetpdf[13];
  hoppetEval(x, Q, hoppetpdf);

  // Print the two PDFs
  cout << "LHAPDF: ";
  for (const auto & val : lhapdf) {
    cout << val << " ";
  }
  cout << endl;

  cout << "Hoppet: ";
  for (const auto & val : hoppetpdf) {
    cout << val << " ";
  }
  cout << endl;


  int npoints = 3000; // or your desired value
  double xmin = 1e-5, xmax = 0.9;
  double qmin = 1.5, qmax = 100000.0;

  // Create log-spaced grids
  vector<double> xvals(npoints), qvals(npoints);
  for (int i = 0; i < npoints; ++i) {
      xvals[i] = exp(log(xmin) + i * (log(xmax) - log(xmin)) / (npoints - 1));
      qvals[i] = exp(log(qmin) + i * (log(qmax) - log(qmin)) / (npoints - 1));
  }

  // Benchmark hoppet_eval
  auto t1 = chrono::high_resolution_clock::now();
  for (double x : xvals)
      for (double q : qvals)
          hoppetEval(x, q, hoppetpdf); 
  auto t2 = chrono::high_resolution_clock::now();
  cout << "hoppetEval time: " << chrono::duration<double>(t2 - t1).count()/npoints/npoints*1e9 << " ns\n";

  // Benchmark xfxq2
  t1 = chrono::high_resolution_clock::now();
  for (double x : xvals)
      for (double q : qvals)
        pdf->xfxQ2(x, q*q, lhapdf);
  t2 = chrono::high_resolution_clock::now();
  cout << "LHAPDF time: " << chrono::duration<double>(t2 - t1).count()/npoints/npoints*1e9 << " ns\n";
}



