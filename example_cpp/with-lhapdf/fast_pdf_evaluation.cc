/// An example in C++ showing how to replace LHAPDF interpolation with
/// (faster) hoppet interpolation
///
/// Usage (from CMake build directory):
///
///   example_cpp/fast_pdf_evaluation [PDFsetName]
///
#include "hoppet.h"
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
  auto t1 = chrono::high_resolution_clock::now();
  pdf = mkPDF(pdfname, imem);
  auto t2 = chrono::high_resolution_clock::now();
  cout << "Time to load LHAPDF set: " << chrono::duration<double,milli>(t2 - t1).count() << " ms" << endl;

  // Now let's access some basic information about the PDF to set up hoppet
  int nloop = pdf->orderQCD() + 1; // LHAPDF is zero indexed
  double xmin = pdf->xMin();
  double xmax = pdf->xMax();
  double Qmin = pdf->qMin();
  double Qmax = pdf->qMax();
  double mc = pdf->quarkMass(4);
  double mb = pdf->quarkMass(5);
  double mt = pdf->quarkMass(6);
  if(!pdf->hasFlavor(6)) mt = 2.0*Qmax; // If no top mass defined, use 173 GeV

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
  t1 = chrono::high_resolution_clock::now();
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme_MSbar); // Start hoppet
  t2 = chrono::high_resolution_clock::now();
  cout << "Time to start hoppet: " << chrono::duration<double,milli>(t2 - t1).count() << " ms" << endl;

  // Now we fill the hoppet grid using the LHAPDF grid directly,
  // rather than evolving ourselves
  double Q0 = Qmin;
  hoppetSetCoupling(pdf->alphasQ(Q0), Q0, nloop);
  t1 = chrono::high_resolution_clock::now();
  hoppetAssign(lhapdf_interface);
  t2 = chrono::high_resolution_clock::now();
  cout << "Time to fill hoppet grid from LHAPDF: " << chrono::duration<double,milli>(t2 - t1).count() << " ms" << endl;

  // If instead we want to evolve the PDF with hoppet starting from
  // some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
  // of hoppetAssign
  //hoppetEvolve(pdf->alphasQ(Q0), Q0, nloop, 1.0, lhapdf_interface, Q0);

}

//----------------------------------------------------------------------
int main (int argc, char *argv[]) {
  string pdfname = "PDF4LHC21_40";
  if (argc > 1) {
    pdfname = argv[1];
  }
  cout << "Using PDF set: " << pdfname << endl;
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
  cout << "PDF Evaluation at x = " << x << ", Q = " << Q << " GeV" << endl;
  cout << "LHAPDF(-6,...,6): ";
  for (const auto & val : lhapdf) cout << val << " ";
  cout << endl;

  cout << "Hoppet(-6,...,6): ";
  for (const auto & val : hoppetpdf) cout << val << " ";
  cout << endl;

  cout << "Difference (H-L): ";
  for (size_t i = 0; i < lhapdf.size(); ++i) {
      cout << hoppetpdf[i] - lhapdf[i] << " ";
  }
  cout << endl;

  int npoints = 3000; // or your desired value
  double xmin = 1e-5, xmax = 0.9;
  double qmin = 1.5, qmax = 7000.0;

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
  cout << "hoppetEval time: " << chrono::duration<double,nano>(t2 - t1).count()/npoints/npoints << " ns\n";
  t1 = chrono::high_resolution_clock::now();
  for (double x : xvals)
      for (double q : qvals)
          hoppetpdf[0] = hoppetEvalPID(x, q, 1); // Just one flavour
  t2 = chrono::high_resolution_clock::now();
  cout << "hoppetEvalPID time: " << chrono::duration<double,nano>(t2 - t1).count()/npoints/npoints << " ns\n";

  // Benchmark xfxq2
  t1 = chrono::high_resolution_clock::now();
  for (double x : xvals)
      for (double q : qvals)
        pdf->xfxQ2(x, q*q, lhapdf);
  t2 = chrono::high_resolution_clock::now();
  cout << "LHAPDF pdf->xfxQ2(x,Q2,vec) time: " << chrono::duration<double,nano>(t2 - t1).count()/npoints/npoints << " ns\n";
  // Benchmark xfxq2
  t1 = chrono::high_resolution_clock::now();
  for (double x : xvals)
      for (double q : qvals)
        lhapdf[0] = pdf->xfxQ2(1, x, q*q);
  t2 = chrono::high_resolution_clock::now();
  cout << "LHAPDF pdf->xfxQ2(id,x,Q2) time: " << chrono::duration<double,nano>(t2 - t1).count()/npoints/npoints << " ns\n";

  // Now check the evaluation time of hoppetAlphaS vs pdf->alphasQ
  // over the range of Q values
  t1 = chrono::high_resolution_clock::now();
  for (double q : qvals)
      hoppetAlphaS(q);
  t2 = chrono::high_resolution_clock::now();
  cout << "hoppetAlphaS time: " << chrono::duration<double,nano>(t2 - t1).count()/npoints << " ns\n";

  t1 = chrono::high_resolution_clock::now();
  for (double q : qvals)
      pdf->alphasQ(q);
  t2 = chrono::high_resolution_clock::now();
  cout << "LHAPDF alphasQ time: " << chrono::duration<double,nano>(t2 - t1).count()/npoints << " ns\n";

  hoppetDeleteAll();

}



