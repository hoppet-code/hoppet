#include "../src/hoppet_v1.h"
#include<iostream>
#include<cmath>

using namespace std;


// definition of the initial condition function
void  heralhc_init(const double & x,
                   const double & Q,
                   double * pdf);


//----------------------------------------------------------------------
int main () {

  double dy    = 0.1;
  int    nloop = 3;

  // initialise with NNLO, VFN
  hoppetStart(dy, nloop);
  
  // evolve the initial condition
  double asQ0 = 0.35, Q0=sqrt(2.0);
  hoppetEvolve(asQ0, Q0, nloop, 1.0, heralhc_init, Q0);
  // alternatively preprepare an evolution and then use its cached version.x
  //hoppetPreEvolve(asQ0, Q0, nloop, 1.0, Q0);
  //hoppetCachedEvolve(heralhc_init);

  // output the results
  double pdf[13];
  double xvals[9]={1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9};
  double Q = 100;
  printf("           Evaluating PDFs at Q = %8.3f GeV\n",Q);
  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon\n");
  for (int ix = 0; ix < 9; ix++) {
    hoppetEval(xvals[ix], Q, pdf);
    printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E\n",xvals[ix],
           pdf[6+2]-pdf[6-2], 
           pdf[6+1]-pdf[6-1], 
           2*(pdf[6-1]+pdf[6-2]),
           (pdf[6-4]+pdf[6+4]),
           pdf[6+0]);
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

