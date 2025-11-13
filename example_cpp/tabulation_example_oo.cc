#include "hoppet.h"
#include "hoppet_oo.h"
#include<iostream>
#include<cmath>
#include<cstdio>
#include <string>

using namespace std;


// definition of the initial condition function
void  heralhc_init(const double & x,
                   const double & Q,
                   vector<double> & pdf);


using namespace hoppet;                   
//----------------------------------------------------------------------
int main () {

  double ymax  = 12.0;  // the maximum y=ln(1/x) value
  double dy    = 0.1;   // the y=ln(1/x) grid spacing
  double dlnlnQ = dy/4; // the table lnlnQ step size
  int    nloop = 3;     // number of loops (LO=1, etc.)

  // set up the grid and the dglap splitting functions etc.
  grid_def grid = grid_def_default(dy, ymax, -6);
  dglap_holder dh(grid, factscheme_MSbar, nloop, 3, 6);

  // set up a running coupling with variable nf
  double mc = 1.414213563;
  double mb = 4.5;
  double mt = 175.0;
  double asQ0 = 0.35;  
  double Q0=sqrt(2.0);
  running_coupling alphas(asQ0, Q0, nloop, mc, mb, mt);

  // set up an initial condition
  pdf pdf0 = pdf_qcd(grid);        // allocate space for a a pdf with QCD flavours
  pdf0.assign_xQ_into(heralhc_init, Q0);   // fill it from the heralhc_init function (below)

  // set up a table (in y=ln1/x and Q)
  double Qmin = 1.2, Qmax = 1e4;
  pdf_table table(grid, Qmin, Qmax, dlnlnQ); // create the table
  table.add_nf_info(alphas);                 // make sure it has the nf boundaries (from the coupling)
  
  // then fill it by evolution from the initial condition
  bool pre_evolve = false;
  if (!pre_evolve) {
    cout << "Using direct evolution, faster with just a single initial condition\n";
    table.evolve(Q0, pdf0, dh, alphas);  // take the initial condition and evolve it to fill the table
  } else {
    cout << "Using pre-evolution, faster for multiple initial conditions\n";
    table.pre_evolve(Q0, dh, alphas);    // "pre-evolve" to set up the evolution operators
    table.evolve(pdf0);                  // apply to initial condition 
  }

  // output the results
  double results[13];
  vector<double> xvals = {1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9};
  double Q = 100;
  printf("           Evaluating PDFs at Q = %8.3f GeV\n",Q);
  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon\n");
  for (double x: xvals) {
    // put the x*pdf(x,Q) values for all flavours into results[]
    table.at_xQ_into(x, Q, results); 
    printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E\n",x,
           results[iflv_u] - results[iflv_ubar], 
           results[iflv_d] - results[iflv_dbar], 
           2*(results[iflv_dbar] + results[iflv_ubar]),
           (results[iflv_c] + results[iflv_cbar]),
           results[iflv_g]);
  }

  // write out the PDF in LHAPDF format (basename/ directory must exist)
  //const std::string basename = "test_cpp";
  //const int pdf_index = 0;
  //hoppetWriteLHAPDFGrid(basename, pdf_index);

}


//----------------------------------------------------------------------
// the initial condition
void  heralhc_init(const double & x,
                   const double & Q,
                   vector<double>& pdf) {

  pdf.resize(13); // ensure pdf has enough space

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

