#include "hoppet_oo.h"
#include "hoppet/qcd.h"
#include "unit-test-helpers.h"
#include <iostream>
#include <cmath>
#include <string>

//#include "catch_amalgamated.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
// define a shorthands for Catch::Matchers::WithinAbs etc.
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// need these to be able to provide our own main
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>

using namespace std;

string red = "\033[1;31m";
string grn = "\033[1;32m";
string yel = "\033[1;33m";
string blu = "\033[1;34m";
string mag = "\033[1;35m";
string cyn = "\033[1;36m";
string bold = "\033[1m";
string endc= "\033[0m";

hoppet::grid_def grid100;
hoppet::grid_def big_grid;

extern "C" void hoppet_cxx__print_c_str(const char * cstr_ptr); // declaration of fortran routine

//-----------------------------------------------------------------------------
/// We supply the main routine, to make sure that we can 
/// do any global-object initialisation that's needed
int main(int argc, char* argv[]) {

  cout << bold << "Starting C++ unit tests with hoppet " << hoppetVersion() << endc << endl;
  string hello="Hello from Fortran to C++";
  hoppet_cxx__print_c_str(hello.c_str()); // test calling Fortran from C++

  // Global setup
  // make a simple grid for use in the tests
  grid100  = hoppet::grid_def(0.1, 10.0);
  big_grid = hoppet::grid_def_default(0.1, 18.0, -6);

  // and also set up the objects in the hoppet streamlined interface
  double ymax = 12.0, dy = 0.2, Qmin = 1.0, Qmax = 1e4;
  double dlnlnQ = dy/4.0;
  int nloop = 3, order = -6;
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order);
  // evolve the initial condition
  double asQ0 = 0.35, Q0=sqrt(2.0), muR_Q = 1.0;
  hoppetEvolve(asQ0, Q0, nloop, muR_Q, hoppetBenchmarkPDFunpol, Q0);

  hoppetEnableCxxExceptionsFromFortran(); // enable C++ exceptions from Fortran

  // then run the tests
  int result = Catch::Session().run(argc, argv);

  // Global teardown (if any)
  hoppetDeleteAll();
  return result;
}


//-----------------------------------------------------------------------------
TEST_CASE( "grid_def", "[hoppet]" ) {  
  REQUIRE(grid100.ptr() != nullptr);
  REQUIRE(grid100.ny() == 100);

  //std::cout << "grid100.ny(): "  << grid100.ny() << std::endl;
  //std::cout << "grid100.ptr(): " << grid100.ptr() << std::endl;

  std::vector<double> xvals = grid100.x_values();
  std::vector<double> yvals = grid100.y_values();

  int iy = grid100.ny() / 2;
  REQUIRE_THAT( yvals[iy], WithinRel(5.0, 1e-12) );
  REQUIRE_THAT( xvals[iy], WithinRel(exp(-5.0), 1e-12) );

  hoppet::split_fn sf(grid100, pqq_fn);
  REQUIRE(sf.grid().ptr() != grid100.ptr());
  REQUIRE(sf.grid() == grid100);
  //hoppet::grid_def grid2 = grid100; // copy constructor
  //REQUIRE( grid100.ptr() != grid2.ptr() ); // different underlying pointers
  //REQUIRE( grid100 == grid2);              // copy should still be equivalent
  //REQUIRE_THAT( grid100.y_values()[iy], WithinRel( grid2.y_values()[iy], 1e-12) );

  // check construction of grid with subgrids & various query operations 
  hoppet::grid_def grid3({grid100, hoppet::grid_def(0.03,0.3)}, true);
  REQUIRE(grid3.nsub()   == 2);
  REQUIRE(grid3.locked() == true);
  REQUIRE(grid3.subgd(2) == grid100);
  REQUIRE(grid3.subiy(1) == 0);
  REQUIRE(grid3.subiy(2) == grid3.subgd(1).ny()+1);

  // now locked=false grid with subgrids, which should preserve order of input grids
  hoppet::grid_def grid4({grid100, hoppet::grid_def(0.03,0.3)}, false);
  REQUIRE(grid4.subgd(1) == grid100);

  //// assignment and equality of grids with subgrids
  //hoppet::grid_def grid5;
  //grid5 = grid3;
  //REQUIRE(grid3.ptr() != grid5.ptr());
  //REQUIRE(grid3 == grid5);

  // further tests of (in)equality
  REQUIRE(!(grid100 == big_grid));
  REQUIRE(grid100 != big_grid);
  hoppet::grid_def null_grid;
  REQUIRE(null_grid != grid100);
  REQUIRE(null_grid != null_grid);
  REQUIRE(grid100 != null_grid);
  // enforcement tests
  REQUIRE_THROWS_AS(null_grid.ensure_valid(), std::runtime_error);
  REQUIRE_THROWS_AS(grid100.ensure_compatible(big_grid), std::runtime_error);
}


//-----------------------------------------------------------------------------
TEST_CASE( "grid_quant", "[hoppet]" ) {

  REQUIRE( grid100.ptr() != nullptr );
  cout << "Defining grid_quant on grid_def, with grid100.ptr(): " << grid100.ptr() << endl;
  hoppet::grid_quant q(grid100);

  q = 4.0;
  REQUIRE_THAT( q.at_y(5.0), WithinAbs(4.0, 1e-12) );

  std::cout << "grid_quant.ptr(): " << q.ptr() << std::endl;
  q = [](double y) { return y*y; };

  // test the self-assignment operator
  q = q;

  // check value is OK (approx, in case we change the grid100 definition)
  REQUIRE_THAT( q.at_y(5.0), WithinAbs(25.0, 1e-6));

  // make sure we can copy things OK and that we get a genuine copy
  int iy = grid100.ny() / 2;
  hoppet::grid_quant q2 = q; // copy constructor
  REQUIRE ( q.ptr() != q2.ptr() ); // different underlying pointers
  REQUIRE ( q[iy] == q2[iy] );     // same values

  hoppet::grid_quant q3(grid100);
  q3.copy(q); // copy 
  REQUIRE ( q.ptr()  != q3.ptr() ); // different underlying pointers
  REQUIRE ( q.data() != q3.data() ); // different underlying data
  REQUIRE ( q[iy] == q3[iy] );
 
  // test move via explicit std::move
  q2 = q;
  auto q2_ptr = q2.ptr();
  q3 = std::move(q2); // assignment operator, should go via a move
  REQUIRE ( q3.ptr() == q2_ptr ); // q3 has taken ownership of q2's pointer
  REQUIRE ( q2.ptr() == nullptr ); // q2 has been nulled

  //------- compound operations ------------------------
  q3 += q;    REQUIRE ( q3[iy] == 2.0 * q[iy] );
  q3 *= 2.0;  REQUIRE ( q3[iy] == 4.0 * q[iy] );
  q3 /= 2.0;  REQUIRE ( q3[iy] == 2.0 * q[iy] );
  q3 -= q;    REQUIRE ( q3[iy] == q[iy] );


  //------- binary operations ------------------------
  hoppet::grid_quant q4;
  q4 = q+q;   REQUIRE ( q4[iy] == 2.0 * q[iy] );
  q4 = q4-q;  REQUIRE ( q4[iy] ==       q[iy] );
  q4 = q+q4;  REQUIRE ( q4[iy] == 2.0 * q[iy] );
  q4 = q*2.0; REQUIRE ( q4[iy] == 2.0 * q[iy] );
  q4 = 2.0*q; REQUIRE ( q4[iy] == 2.0 * q[iy] );
  q4 = q/2.0; REQUIRE ( q4[iy] == 0.5 * q[iy] );

  //-- test operations with views ------------------------
  hoppet::grid_quant_view qv(q); // "view" copy constructor (effectively a reference)
  REQUIRE( q.data() == qv.data() ); 
  q4 += qv;   REQUIRE ( q4[iy] == 1.5 * q[iy] );
  q2 = q;
  hoppet::grid_quant_view qv2(q2); // "view" copy constructor (effectively a reference)
  qv2 += q2;   REQUIRE ( qv2[iy] == 2.0 * q[iy] );
  qv2 += 2*q;  REQUIRE ( qv2[iy] == 4.0 * q[iy] );
  q3 = qv2 + q; REQUIRE ( q3[iy] == 5.0 * q[iy] );
  //qv = q + q4; REQUIRE ( qv[iy] == 5.0 * q[iy] );

  /// this should assign the contents of qv2 into qv
  qv = qv2;
  REQUIRE (qv.data() != qv2.data() ); // different underlying data
  REQUIRE (qv[iy] == qv2[iy] );
  qv.take_view(qv2);
  REQUIRE (qv.data() == qv2.data() ); // same underlying data
  hoppet::grid_quant_view qv3 (q);

  //------- moments ------------------------
  auto qdist = big_grid * [](double y) { double x = exp(-y); return 5*pow(1-x,4)*x;};
  REQUIRE_THAT(qdist.truncated_moment(0.0), WithinAbs(1.0, 1e-5));
  REQUIRE_THAT(qdist.truncated_moment(1.0), WithinAbs(1/6.0, 1e-7));
  double ytrunc = 5.0; REQUIRE(qdist.truncated_moment(0.0, ytrunc) < qdist.truncated_moment(0.0));

  //------- luminosity ------------------------
  qdist = big_grid * [](double y) { double x = exp(-y); return pow(1-x,4);};
  auto lumi = qdist.luminosity(qdist);
  auto expected_lumi = [](double y) { 
    double x = exp(-y);
    // function specfic to (1-x)^4 with (1-x)^4 luminosity
    return -(log(x)*(1 + 16*x + 36*pow(x,2) + 16*pow(x,3) + pow(x,4))) + 
            (5*(-5 - 32*x + 32*pow(x,3) + 5*pow(x,4)))/6.;
  };
  for (double y : {10.0, 8.0, 5.0, 2.0}) {
    REQUIRE_THAT( lumi.at_y(y), WithinRel( expected_lumi(y), 1e-5) );
  }

}


//-----------------------------------------------------------------------------
TEST_CASE( "grid_quant_2d", "[hoppet]" ) {
  hoppet::grid_quant_2d pdf(grid100, 14); // 14 is size of 2nd dimension

  pdf = 4.0; // set all values to 4.0
  REQUIRE_THAT( pdf[hoppet::iflv_g   ].at_y(5.0), WithinAbs(4.0, 1e-12) );

  hoppet::grid_quant_2d_view pdf_view(pdf);
  pdf_view = 3.0;
  REQUIRE_THAT( pdf[hoppet::iflv_g   ].at_y(5.0), WithinAbs(3.0, 1e-12) );
  pdf = 0.0;


  pdf[hoppet::iflv_g   ] =  grid100 * [](double y) { return y*y; };
  pdf[hoppet::iflv_dbar] =            [](double y) { return pow(y,3); };
  REQUIRE_THAT( pdf[hoppet::iflv_g   ].at_y(5.0), WithinAbs(25.0, 1e-6));
  REQUIRE_THAT( pdf[hoppet::iflv_dbar].at_y(5.0), WithinAbs(125.0, 1e-6));


  hoppet::grid_quant_2d pdf2 = pdf; // copy constructor
  auto pdf3 = pdf + pdf2;

  hoppet::grid_quant_2d pdf4(grid100, 14);
  //pdf4.assign_xQ_into(hoppetBenchmarkPDFunpol, 10.0);
  //cout << "g    " << pdf4[hoppet::iflv_g   ].at_y(5.0) << endl;
  //cout << "ubar " << pdf4[hoppet::iflv_ubar].at_y(5.0) << endl;
}

//-----------------------------------------------------------------------------
TEST_CASE( "pdf_qcd", "[hoppet]" ) {
  hoppet::pdf pdf = pdf_qcd(big_grid);
  double dummy_Q = 0.0;
  pdf.assign_xQ_into(hoppetBenchmarkPDFunpol, dummy_Q);
  // hard-coded values from the benchmark function
  auto fn_g    = [](double x) { return 1.7*pow(x,-0.1) * pow(1-x,5); };
  auto fn_ubar = [](double x) { return 0.387975*0.5*pow(x,-0.1)*pow(1-x,7); };
  auto xvals = {1e-4, 1e-2, 0.1, 0.5};
  for (double x : xvals) {
    REQUIRE_THAT( pdf[hoppet::iflv_g   ].at_x(x), WithinAbs(fn_g(x), 1e-5));
    REQUIRE_THAT( pdf[hoppet::iflv_ubar].at_x(x), WithinAbs(fn_ubar(x), 1e-6));
  }
  auto control = pdf[hoppet::ncompmax];
  bool allzero = true;
  for (int i = 0; i <= pdf.grid().ny(); ++i) {
    if (control[i] != 0.0) allzero = false;
  }
  REQUIRE(allzero);

  hoppet::pdf pdf2 = pdf + pdf;
  for (double x : xvals) {
    REQUIRE_THAT( pdf2[hoppet::iflv_g   ].at_x(x), WithinAbs(2*fn_g(x), 1e-5));
  }
  //cout << pdf2.ptr() << "=pdf2.ptr()\n";
}



//-----------------------------------------------------------------------------
TEST_CASE( "grid_conv", "[hoppet]" ) {
  hoppet::grid_conv pqq(big_grid, pqq_fn);

  hoppet::grid_conv pgq(big_grid, [](double y, int piece) {
      double x = exp(-y);
      double res;
      if      (piece == hoppet::cc_REAL
            || piece == hoppet::cc_REALVIRT) res = cf * (1+pow(1-x,2))/x;
      else                               res = 0.0;
      return x*res;
    });

  auto power = [](double powval) { return grid100 * [powval](double y) { return exp(-powval*y); }; };

  auto q = big_grid * [](double y) { double x = exp(-y); return 5*pow(1-x,4)*x;};
  auto pqq_q = pqq * q;
  auto pgq_q = pgq * q;
  double pqq_q_mom1 = pqq_q.truncated_moment(1.0);
  double pgq_q_mom1 = pgq_q.truncated_moment(1.0);
  double q_mom1 = 1.0/6.0;
  double pqq_mom1 = -cf * 4.0/3.0;
  REQUIRE_THAT(pqq_q_mom1/q_mom1, WithinAbs(pqq_mom1, 1e-6)); //< check momentum moment (1/6 * (-4/3))
  REQUIRE_THAT(          pqq_q.truncated_moment(0.0), WithinAbs(0.0, 1e-5));    //< check quark number conserved
  REQUIRE_THAT((pqq_q + pgq_q).truncated_moment(1.0), WithinAbs(0.0, 1e-6));    //< check momentum conserved

  //------- compound operations ------------------------
  hoppet::grid_conv p2 = pqq;
  REQUIRE(p2.ptr() != nullptr);
  REQUIRE(p2.ptr() != pqq.ptr());
  REQUIRE(p2.grid() == pqq.grid()); // grids should be equivalent
  p2 += pgq;
  REQUIRE_THAT( (p2 * q) .truncated_moment(1.0), WithinAbs(0.0, 1e-6));    //< check momentum conserved
  p2 -= pgq;
  REQUIRE_THAT( (p2 * q) .truncated_moment(1.0), WithinAbs( pqq_q_mom1, 1e-6)); //< check back to pqq alone
  p2 *= 2.0;
  REQUIRE_THAT( (p2 * q) .truncated_moment(1.0), WithinAbs( 2.0 * pqq_q_mom1, 1e-6)); //< check back to pqq alone
  p2 /= 2.0;
  REQUIRE_THAT( (p2 * q) .truncated_moment(1.0), WithinAbs( pqq_q_mom1, 1e-6)); //< check back to pqq alone

  //------- binary operations ------------------------
  hoppet::grid_conv p3 = pqq + pgq;
  REQUIRE_THAT( (p3 * q) .truncated_moment(1.0), WithinAbs(0.0, 1e-6));    //< check momentum conserved
  p3 = p3 - pgq;
  REQUIRE_THAT( (p3 * q) .truncated_moment(1.0), WithinAbs( pqq_q_mom1, 1e-6)); //< check back to pqq alone
  p3 = p3 * 2.0;
  REQUIRE_THAT( (p3 * q) .truncated_moment(1.0), WithinAbs( 2.0 * pqq_q_mom1, 1e-6)); //< check back to pqq alone
  p3 = 2.0 * p3;
  REQUIRE_THAT( (p3 * q) .truncated_moment(1.0), WithinAbs( 4.0 * pqq_q_mom1, 1e-6)); //< check back to pqq alone
  p3 = p3 / 4.0;
  REQUIRE_THAT( (p3 * q) .truncated_moment(1.0), WithinAbs( pqq_q_mom1, 1e-6)); //< check back to pqq alone

  auto pqq_pqq = pqq * pqq;
  auto pqq_pqq_q = pqq_pqq * q;
  REQUIRE_THAT( (pqq_pqq_q).truncated_moment(1.0), WithinAbs(pqq_mom1 * pqq_mom1 * q_mom1, 1e-6));
  REQUIRE_THAT( (pqq_pqq_q).at_x(0.3), WithinAbs( (pqq * (pqq * q)).at_x(0.3), 1e-6) );



  //------- checking throw on incompatible grids ------------------------
  hoppet::grid_conv pqq_100(grid100, pqq_fn);
  REQUIRE_THROWS(pqq_100 += pqq); // self-add, should give an error because of incompatible grid


  //-------- view tests ------------------------
  hoppet::grid_conv pgen = pqq; // general copy
  REQUIRE_THAT( (pgen  * q) .truncated_moment(1.0), WithinAbs(pqq_q_mom1, 1e-6)); // since pview was a view of pqq, pqq should also now be pgq
  hoppet::grid_conv_view pview(pgen); // view of general copy
  auto pgen_ptr = pgen.ptr();
  pview = pgq;
  CHECK( pview.ptr() == pgen_ptr ); // pview should still point to pgen's data
  REQUIRE_THAT( (pview * q) .truncated_moment(1.0), WithinAbs(pgq_q_mom1, 1e-6)); // pview should be equal to pgq
  REQUIRE_THAT( (pgen  * q) .truncated_moment(1.0), WithinAbs(pgq_q_mom1, 1e-6)); // since pview was a view of pgen, pgen should also now be pgq
  pview += pgq;
  REQUIRE_THAT( (pgen * q) .truncated_moment(1.0), WithinAbs(2*pgq_q_mom1, 1e-6)); // pgen should be equal to 2*pgq
}

//-----------------------------------------------------------------------------
TEST_CASE( "split_mat", "[hoppet]" ) {
  int nf = 5;
  hoppet::split_mat p_lo(nf);
  REQUIRE(p_lo.nf() == nf);

  using hoppet::iflv_g;

  // assign the various components
  p_lo.qq() = hoppet::grid_conv(big_grid, pqq_fn);
  p_lo.gq() = hoppet::grid_conv(big_grid, pgq_fn);
  p_lo.gg() = hoppet::grid_conv(big_grid, [nf](double y, int piece) { return pgg_fn(y, piece, nf); });
  p_lo.qg() = hoppet::grid_conv(big_grid, [nf](double y, int piece) { return pqg_fn(y, piece, nf); });
  p_lo.ns_minus() = p_lo.qq(); 
  p_lo.ns_plus()  = p_lo.qq();
  p_lo.ns_v()     = p_lo.qq();
  REQUIRE( p_lo.qq().ptr() != p_lo.ns_v().ptr() ); // different underlying pointers

  auto p_copy = p_lo; // copy constructor

  hoppet::pdf pdf = pdf_qcd(big_grid);  
  // check assignment with the classic LHAPDF interface
  double dummy_Q = 0.0;
  pdf.assign_xQ_into(hoppetBenchmarkPDFunpol, dummy_Q);

  // check assignment with the LHAPDF PDF::xfxQ(x,Q,v) interface
  hoppet::pdf pdf2 = pdf_qcd(big_grid);
  auto benchmark_via_vec = [](double x, double Q, vector<double> & v){
    v.resize(13);
    hoppetBenchmarkPDFunpol(x, Q, v.data());
  };
  pdf2.assign_xQ_into(benchmark_via_vec, dummy_Q);
  for (double x : {1e-4, 1e-2, 0.1, 0.5}) {
    for (int iflv = hoppet::iflv_min; iflv <= hoppet::iflv_max; ++iflv)
      REQUIRE( pdf[iflv].at_x(x) == pdf2[iflv].at_x(x) );
  }

  //  pdf[hoppet::iflv_g] *= 0.0; // zero gluon to make life simpler
  //pdf *= 0.0; // zero everything to make life simpler
  //pdf[hoppet::iflv_g ] = big_grid * [](double y) { double x = exp(-y); return 2*pow(1-x,4)*x; };
  //pdf[hoppet::iflv_ubar] = big_grid * [](double y)


  hoppet::pdf dpdf_lo = p_lo * pdf;

  auto mom0 = [](const auto & qv) { return qv.truncated_moment(0.0); };
  auto mom1 = [](const auto & qv) { return qv.truncated_moment(1.0); };

  // check u-ubar number sum rule
  REQUIRE_THAT( mom0(dpdf_lo[hoppet::iflv_u ]-dpdf_lo[hoppet::iflv_ubar]), WithinAbs(0.0, 1e-4) );

  // check momentum sum rule
  hoppet::grid_quant d_all_quarks = dpdf_lo[hoppet::iflv_tbar];
  for (int i = hoppet::iflv_tbar+1; i < hoppet::ncompmax; ++i) {
    if (i != hoppet::iflv_g) d_all_quarks += dpdf_lo[i];
    //cout << "i=" << i << " mom1=" << mom1(dpdf_lo[i]) << endl;
  }
  REQUIRE_THAT( mom1(d_all_quarks + dpdf_lo[hoppet::iflv_g]), WithinAbs(0.0, 1e-4) );

  //--- test addition of two split_mats
  auto pp = p_lo;
  pp += p_lo;             // compound add
  auto pp2 = pp + 2*p_lo; // binary add
  REQUIRE_THAT( mom1((pp *pdf)[iflv_g]), WithinAbs(2*mom1(dpdf_lo[iflv_g]), 1e-6) );
  REQUIRE_THAT( mom1((pp2*pdf)[iflv_g]), WithinAbs(4*mom1(dpdf_lo[iflv_g]), 1e-6) );
  pp2 = -pp2 - p_lo/2;      // binary subtract
  REQUIRE_THAT( mom1((pp2*pdf)[iflv_g]), WithinAbs(-4.5*mom1(dpdf_lo[iflv_g]), 1e-6) );

  // test multiplication of two split_mats
  auto ppp = p_lo * p_lo;
  auto ppp_pdf = ppp * pdf;
  REQUIRE_THAT( mom1(ppp_pdf[iflv_g]), WithinAbs(mom1((p_lo * dpdf_lo)[iflv_g]), 3e-5) );

  //------- test commutator of two split_mats
  hoppet::split_mat p_lo_swapped = p_lo;
  // change some components to make it different
  p_lo_swapped.qg() = p_lo.gq();
  p_lo_swapped.gq() = p_lo.qg();
  // get comm
  auto p_comm = hoppet::commutator(p_lo, p_lo_swapped);
  auto p_comm_pdf = p_comm * pdf;
  auto p_comm_pdf_direct = (p_lo * (p_lo_swapped * pdf)) - (p_lo_swapped * (p_lo * pdf));
  REQUIRE_THAT( mom1(p_comm_pdf[iflv_g]), WithinAbs(mom1(p_comm_pdf_direct[iflv_g]), 3e-5) );

  hoppet::grid_quant_2d pdf_big(big_grid, 21); // 21 is size of 2nd dimension (more than standard 14)
  pdf_big.assign(0.01);
  pdf_big[hoppet::ncompmax] = 0.0;
  auto p_lo_pdf_big = p_lo * pdf_big; // test multiplication with larger pdf_2d
  REQUIRE( p_lo_pdf_big.size()  == pdf_big.size() );
  REQUIRE( p_lo_pdf_big[15].at_y(5.0) == 0.0); // should be set to zero
  REQUIRE( p_lo_pdf_big[hoppet::iflv_g].at_y(5.0) != 0.0); // should not be zero

  hoppet::grid_quant_2d pdf_sml(big_grid, 7); // 
  pdf_sml.assign(0.0);
  REQUIRE_THROWS_AS( p_lo * pdf_sml, std::runtime_error ); // should throw because too small

  // test copying from the streamlined interface
  int nloop = 2;
  hoppet::split_mat p_gen = hoppet::sl::dh.p(nloop, nf);
}

//-----------------------------------------------------------------------------
TEST_CASE( "mass_threshold_mat", "[hoppet]" ) {
  typedef hoppet::mass_threshold_mat MTM;
  MTM::view_t mtm_view;
  int nf_heavy = 4;
  int nloop    = 3;
  MTM mtm(nf_heavy);
  auto & grid = hoppet::sl::grid;
  auto & dh   = hoppet::sl::dh;
  mtm_view.take_view(dh.mtm(nloop, nf_heavy));
  mtm = mtm_view; // test assignment from view

  // now do some basic logic checks, namely that heavy-quark entries are zero/non-zero as appropriate
  hoppet::pdf pdf = pdf_qcd(grid);
  pdf = 0.0;
  pdf[hoppet::iflv_g] = grid * [](double y) { double x = exp(-y); return 2*pow(1-x,6)*x; };
  auto dpdf = mtm * pdf;
  REQUIRE(dpdf[hoppet::iflv_c].at_y(5.0) != 0.0);
  REQUIRE(dpdf[hoppet::iflv_b].at_y(5.0) == 0.0);

  // try replacing a one of the entries in the mtm and check that it propagates
  mtm.sgg_h() = dh.p_lo().gg(); // set to LO gg->gg (physically wrong, but serves the purpose)
  auto ddpdf = mtm * pdf;
  REQUIRE(ddpdf[hoppet::iflv_g].at_y(5.0) != dpdf[hoppet::iflv_g].at_y(5.0));
  REQUIRE(ddpdf[hoppet::iflv_c].at_y(5.0) == dpdf[hoppet::iflv_c].at_y(5.0));
}

//-----------------------------------------------------------------------------
TEST_CASE( "qcd", "[hoppet]" ) {
  int nf_store = hoppet::qcd::nf_int;

  REQUIRE(hoppet::qcd::ca == 3.0);
  REQUIRE(hoppet::qcd::tr == 0.5);
  REQUIRE_THAT(hoppet::qcd::cf, WithinAbs(4.0/3.0, 1e-14));

  hoppet::qcd::set_nf(3);
  REQUIRE(hoppet::qcd::nf_int == 3);
  REQUIRE(hoppet::qcd::nf_u   == 1);
  REQUIRE(hoppet::qcd::nf_d   == 2);
  REQUIRE(hoppet::qcd::nf     == 3.0);

  hoppet::qcd::set_group(4.0, 1.0, 0.25);
  REQUIRE(hoppet::qcd::ca == 4.0);
  REQUIRE(hoppet::qcd::cf == 1.0);
  REQUIRE(hoppet::qcd::tr == 0.25);

  // restore standard QCD values
  hoppet::qcd::set_group(3.0, 4.0/3.0, 0.5);
  hoppet::qcd::set_nf(nf_store);
}

//-----------------------------------------------------------------------------
TEST_CASE( "running_coupling", "[hoppet]" ) {
  // first with fixnf
  double asmz = 0.118;
  double mz = 91.1880;
  int nloop = 3;
  int nf = 5;
  hoppet::running_coupling alphas(asmz, mz, nloop, nf);

  // deferred assignment
  hoppet::running_coupling alphas2;
  REQUIRE_THROWS_AS(alphas2(mz), std::runtime_error);
  alphas2 = hoppet::running_coupling(asmz, mz, nloop, nf); // move assignment?

  REQUIRE_THAT(alphas(mz), WithinAbs(asmz, 1e-6));
  REQUIRE(alphas2(mz) == alphas(mz));

  //alphas2 = alphas; // test copy assignment operator
  //auto alphas3(alphas); // test copy constructor

  hoppet::running_coupling alphas3 = std::move(alphas2); // test move via std::move
  REQUIRE(alphas3.ptr() != nullptr);
  REQUIRE(alphas2.ptr() == nullptr);

  hoppet::running_coupling::view_t alphas_view;
  alphas_view.take_view(alphas);
  REQUIRE(alphas_view(mz) == alphas(mz));
  REQUIRE(alphas_view.ptr() == alphas.ptr());

  // now with variable nf ---------------
  double mc = 1.67, mb = 4.78, mt = 172.57;
  nloop = 4;
  // reference results generated by unit-tests/gen-alphas-ref.py (as used in test_alphas.F90)
  vector<double> Qvals = {1.000000, 3.000000, 9.000000, 10000.000000};
  vector<double> ref_4loop_vfn = {4.830605117006e-01, 2.536789691195e-01, 1.827498038423e-01, 7.181752397345e-02};
  alphas = hoppet::running_coupling(asmz, mz, nloop, mc, mb, mt); // variable nf
  REQUIRE_THAT(alphas(mz), WithinAbs(asmz, 1e-6));
  REQUIRE(alphas(1.0) > alphas3(1.0)); // should be larger at low Q, because of fewer flavours

  for (unsigned i = 0; i < 4; ++i) {
    REQUIRE_THAT(alphas(Qvals[i]), WithinAbs(ref_4loop_vfn[i], 1e-9));
  }

  REQUIRE(alphas.nloop() == nloop);

  REQUIRE_THAT(alphas.quark_mass(4)             , WithinAbs(mc, 1e-12));
  REQUIRE_THAT(alphas.quark_mass(hoppet::iflv_c), WithinAbs(mc, 1e-12));

  auto [nflo, nfhi] = alphas.nf_range();
  REQUIRE(nflo == 3);
  REQUIRE(nfhi == 6);

  double nf_at_1 = alphas.nf_at_Q(1.0);
  REQUIRE(nf_at_1 == 3);

  double nf_at_3 = alphas.nf_at_Q(3.0);
  REQUIRE(nf_at_3 == 4);

  double Qlo, Qhi;
  double nf_at_5 = alphas.nf_at_Q(5.0, Qlo, Qhi);
  REQUIRE(nf_at_5 == 5);
  REQUIRE_THAT(Qlo, WithinAbs(mb, 1e-12));
  REQUIRE_THAT(Qhi, WithinAbs(mt, 1e-12));
  std::tie(Qlo, Qhi) = alphas.Q_range_at_nf(4);
  REQUIRE_THAT(Qlo, WithinAbs(mc, 1e-12));
  REQUIRE_THAT(Qhi, WithinAbs(mb, 1e-12));

  bool msbar = false;
  double muMatch_mQuark = 2.0;
  alphas2 = hoppet::running_coupling(asmz, mz, nloop, mc, mb, mt, msbar, muMatch_mQuark); // variable nf
  // check that the matching point is at 2*mQuark
  std::tie(Qlo, Qhi) = alphas2.Q_range_at_nf(4);
  int nf_at_12_3 = alphas.nf_at_Q(12.0, Qlo, Qhi, 3.0);
  REQUIRE(nf_at_12_3 == 4);
  REQUIRE_THAT(Qlo, WithinAbs(3*mc, 1e-12));
  REQUIRE_THAT(Qhi, WithinAbs(3*mb, 1e-12));
  // sanity check: MSbar masses would be lower, but we've kept them the same numerical values as the pole
  // which means the switch is effectively higher, giving more running 

}

//-----------------------------------------------------------------------------
TEST_CASE("dglap_holder", "[hoppet]") {
  using namespace hoppet;
  int nloop = 1;
  int nflo = 4;
  int nfhi = 5;
  dglap_holder dh(big_grid, hoppet::factscheme_MSbar, nloop, nflo, nfhi);
  dh.set_nf(5);

  REQUIRE(dh.nloop() == nloop);

  pdf_flav q = dh.grid() * [](double y) { double x = exp(-y); return 5*pow(1-x,4)*x;};
  pdf_flav plo_qq_q1 = dh.p_lo().qq() * q;
  
  split_fn pqq = dh.grid() * pqq_fn;
  pdf_flav plo_qq_q2 = pqq * q;

  REQUIRE_THAT(plo_qq_q1.at_y(5.05), WithinAbs(plo_qq_q2.at_y(5.05), 1e-6));
}


//-----------------------------------------------------------------------------
TEST_CASE("pdf_table", "[hoppet]") {

  hoppet::pdf_table tab = setup_table(big_grid);

  auto sl_tab = hoppet::sl::table;
  size_t iQ = sl_tab.nQ() / 2;
  double Q = sl_tab.Q_vals(iQ);
  int nf = sl_tab.nf_int(iQ);
  int nloop = 1;
  double x = 0.09;

  // try various evaluations to check equivalence between the streamlined interface
  // and explicit local evaluations here
  double gluon_sl = hoppetEvalIFlv(x,Q,hoppet::iflv_g);
  double gluon_oo = sl_tab.at_iQ(iQ)[hoppet::iflv_g].at_x(x);
  hoppet::pdf_flav_view gluon_flav = sl_tab.at_iQ(iQ)[hoppet::iflv_g];
  hoppet::pdf_flav_view gluon_flav_alt = sl_tab.at_iQf(iQ,hoppet::iflv_g);
  REQUIRE_THAT(gluon_sl, WithinAbs(gluon_flav.at_x(x), 1e-6));
  REQUIRE((gluon_flav.data() == gluon_flav_alt.data()));

  // access via xQf
  gluon_oo = sl_tab.at_xQf(x, Q, hoppet::iflv_g);
  REQUIRE_THAT(gluon_sl, WithinAbs(gluon_oo, 1e-6));

  // access via at_xQ_into
  double f_into[hoppet::ncompmax];
  sl_tab.at_xQ_into(x, Q, f_into);
  REQUIRE_THAT(gluon_sl, WithinAbs(f_into[hoppet::iflv_g], 1e-6));


  // explore applying splitting functions, locally versus streamlined
  hoppet::pdf p_pdf = hoppet::sl::dh.p(nloop,nf) * sl_tab.at_iQ(iQ);
  vector<double> p_pdf_raw(sl_tab.size_flv());
  hoppetEvalSplit(x, Q, nloop, nf, p_pdf_raw.data());
  for (int i = 0; i <= sl_tab.iflv_max(); ++i) {
    REQUIRE_THAT(p_pdf[i].at_x(x), WithinAbs(p_pdf_raw[i], 1e-6));
  }

  // now try to evolve a table  
  hoppet::pdf_table tab_copy = hoppet::sl::table; // copy constructor
  hoppet::pdf_table_view tab_view(tab_copy);
  REQUIRE(tab_view.ptr() != hoppet::sl::table.ptr()); // different underlying pointers
  hoppet::pdf pdf_Q0 = hoppet::pdf_qcd(tab_view.grid());
  double Q0 = sqrt(2.0);
  pdf_Q0.assign_xQ_into(hoppetBenchmarkPDFunpol, 0.0);
  tab_view.evolve(Q0, pdf_Q0, hoppet::sl::dh, hoppet::sl::coupling);
  REQUIRE_THAT( tab_view.at_xQf(x, Q, hoppet::iflv_g), WithinAbs( hoppetEvalIFlv(x,Q,hoppet::iflv_g), 1e-6) );
  //tab_view.write_lhapdf(hoppet::sl::coupling, "/tmp/base");

  // next, try pre-evolution
  hoppet::pdf_table tab_pre = hoppet::sl::table; // copy constructor
  tab_pre.pre_evolve(Q0, hoppet::sl::dh, hoppet::sl::coupling);
  tab_pre.evolve(pdf_Q0);
  REQUIRE_THAT( tab_pre.at_xQf(x, Q, hoppet::iflv_g), WithinAbs( hoppetEvalIFlv(x,Q,hoppet::iflv_g), 1e-6) );

  cout << hoppet::sl::table.nQ() << "=hoppet::sl::table.nQ()\n";
}

//-----------------------------------------------------------------------------
TEST_CASE("pdf_table_assignment", "[hoppet]") {

  auto silly_pdf_dblptr = [&](double x, double Q, double * f){
    for (int iflv = hoppet::iflv_min; iflv <= hoppet::iflv_max; ++iflv) {
      f[iflv] = silly_pdf(x, Q, iflv);
    }
  };
  auto silly_pdf_dblvec = [&](double x, double Q, vector<double> & f){
    f.resize(hoppet::iflv_max+1);
    for (int iflv = hoppet::iflv_min; iflv <= hoppet::iflv_max; ++iflv) {
      f[iflv] = silly_pdf(x, Q, iflv);
    }
  };

  auto tab1 = setup_table(big_grid);
  tab1.assign_xQ_into(silly_pdf_dblptr);
  auto tab2 = setup_table(big_grid);
  tab2.assign_xQ_into(silly_pdf_dblvec);

  vector<double> xvals = {1e-4, 1e-2, 0.1, 0.5};
  vector<double> Qvals = {2.0, 10.0, 100.0, 1000.0};
  for (double x : xvals) {
    for (double Q : Qvals) {
      for (int iflv = hoppet::iflv_min; iflv <= tab1.iflv_max(); ++iflv) {
        REQUIRE_THAT( tab1.at_xQf(x, Q, iflv), WithinAbs(silly_pdf(x,Q,iflv), 1e-6) );
        REQUIRE_THAT( tab2.at_xQf(x, Q, iflv), WithinAbs(silly_pdf(x,Q,iflv), 1e-6) );
      }
    }
  }

  // check the control component is zero
  int iy_test = 1;
  int iQ = tab1.nQ() / 2;
  for (auto & tab : {tab1, tab2}) {
    auto control = tab.at_iQ(iQ)[hoppet::ncompmax];
    //for (unsigned i = 0; i < control.size(); ++i) {cout << i << " " << control[i] << "\n";} 
    REQUIRE(control[iy_test] == 0);
  }
}

TEST_CASE("pdf_table_array", "[hoppet]") {

  int ntab = 3;
  hoppet::pdf_table_array tab_array(ntab);
  tab_array[0] = hoppet::sl::tables[0];
  tab_array[1] = hoppet::sl::tables[0];
  tab_array[2] = hoppet::sl::tables[0];
  // check we've copied things
  REQUIRE(tab_array[0].ptr() != hoppet::sl::tables[0].ptr());
  REQUIRE(tab_array[0].ptr() != tab_array[1].ptr());
  // and that structure is the same
  REQUIRE(tab_array[2].nQ() == hoppet::sl::tables[0].nQ());

  // modify entries in different tables
  for (int iQ = 0; iQ <= tab_array[0].nQ(); ++iQ) {
    tab_array[1](iQ,hoppet::iflv_g) *= 2.0; // double gluon in tab 1
    tab_array[2](iQ,hoppet::iflv_g) *= 3.0; // triple gluon in tab 2
  }
  double x = 0.1;
  double Q = 100.0;
  REQUIRE_THAT(tab_array[2].at_xQf(x,Q,hoppet::iflv_g), WithinAbs(3.0*hoppetEvalIFlv(x,Q,hoppet::iflv_g), 1e-6));
  REQUIRE_THAT(tab_array[1].at_xQf(x,Q,hoppet::iflv_g), WithinAbs(2.0*hoppetEvalIFlv(x,Q,hoppet::iflv_g), 1e-6));
  REQUIRE_THAT(tab_array[0].at_xQf(x,Q,hoppet::iflv_g), WithinAbs(1.0*hoppetEvalIFlv(x,Q,hoppet::iflv_g), 1e-6));
}

//-----------------------------------------------------------------------------
TEST_CASE( "streamlined-objects", "[hoppet]" ) {
  hoppet::grid_quant q(hoppet::sl::grid);
  hoppet::dglap_holder_view dh = hoppet::sl::dh;
  //hoppet::grid_quant q(hoppet_sl_grid_ptr);
  q = [](double y) { return y*y; };
  REQUIRE_THAT( q.at_y(5.0), WithinAbs(25.0, 1e-6));


  auto qq = q;
  //q.take_view(qq);

  dh.set_nf(4);
  REQUIRE( dh.nf() == 4 );
  REQUIRE( dh.nloop() == 3 );

  // check pdf tables
  REQUIRE( hoppet::sl::tables[0].ptr() == hoppet::sl::table.ptr() );
}
