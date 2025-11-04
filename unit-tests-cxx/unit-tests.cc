#include "hoppet_oo.h"
#include <iostream>
#include <cmath>

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

hoppet::grid_def grid100;
hoppet::grid_def big_grid;

/// We supply the main routine, to make sure that we can 
/// do any global-object initialisation that's needed
int main(int argc, char* argv[]) {
  // Global setup
  // make a simple grid for use in the tests
  grid100  = hoppet::grid_def(0.1, 10.0);
  big_grid = hoppet::grid_def_default(0.1, 16.0, -6);

  // and also set up the objects in the hoppet streamlined interface
  double ymax = 12.0, dy = 0.2, Qmin = 1.0, Qmax = 1e4;
  double dlnlnQ = dy/4.0;
  int nloop = 1, order = -6;
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order);

  // then run the tests
  int result = Catch::Session().run(argc, argv);

  // Global teardown (if any)
  hoppetDeleteAll();
  return result;
}


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

  hoppet::grid_def grid2 = grid100; // copy constructor
  REQUIRE( grid100.ptr() != grid2.ptr() ); // different underlying pointers
  REQUIRE( grid100 == grid2);              // copy should still be equivalent
  REQUIRE_THAT( grid100.y_values()[iy], WithinRel( grid2.y_values()[iy], 1e-12) );

  hoppet::grid_def grid3({grid100, hoppet::grid_def(0.03,0.3)}, true);
  hoppet::grid_def grid4;
  grid4 = grid3;

  // tests of (in)equality
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

}


//-----------------------------------------------------------------------------
TEST_CASE( "grid_quant_2d", "[hoppet]" ) {
  hoppet::grid_quant_2d pdf(grid100, 14); // 14 is size of 2nd dimension
  pdf[hoppet::iflv_g   ] =  grid100 * [](double y) { return y*y; };
  pdf[hoppet::iflv_dbar] =            [](double y) { return pow(y,3); };
  REQUIRE_THAT( pdf[hoppet::iflv_g   ].at_y(5.0), WithinAbs(25.0, 1e-6));
  REQUIRE_THAT( pdf[hoppet::iflv_dbar].at_y(5.0), WithinAbs(125.0, 1e-6));

  hoppet::grid_quant_2d pdf2 = pdf; // copy constructor
  auto pdf3 = pdf + pdf2;

  hoppet::grid_quant_2d pdf4(grid100, 14);
  //pdf4.assign(hoppetBenchmarkPDFunpol, 10.0);
  //cout << "g    " << pdf4[hoppet::iflv_g   ].at_y(5.0) << endl;
  //cout << "ubar " << pdf4[hoppet::iflv_ubar].at_y(5.0) << endl;
}

TEST_CASE( "pdf_qcd", "[hoppet]" ) {
  hoppet::pdf pdf = pdf_qcd(big_grid);
  double dummy_Q = 0.0;
  pdf.assign(hoppetBenchmarkPDFunpol, dummy_Q);
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
TEST_CASE( "streamlined-objects", "[hoppet]" ) {
  hoppet::grid_quant q(hoppet::sl::grid);
  //hoppet::grid_quant q(hoppet_sl_grid_ptr);
  q = [](double y) { return y*y; };
  REQUIRE_THAT( q.at_y(5.0), WithinAbs(25.0, 1e-6));
}

//-----------------------------------------------------------------------------
TEST_CASE( "grid_conv", "[hoppet]" ) {
  auto pqq_fn = [](double y, int piece) {
      double x = exp(-y);
      double res;
      if      (piece == hoppet::cc_REAL) res =  (1+x*x)/(1-x);
      else if (piece == hoppet::cc_VIRT) res = -(1+x*x)/(1-x);
      else                               res = 0.0;
      return x*res;
    };  
  hoppet::grid_conv pqq(big_grid, pqq_fn);

  hoppet::grid_conv pgq(big_grid, [](double y, int piece) {
      double x = exp(-y);
      double res;
      if      (piece == hoppet::cc_REAL
            || piece == hoppet::cc_REALVIRT) res = (1+pow(1-x,2))/x;
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
  double pqq_mom1 = -4.0/3.0;
  REQUIRE_THAT(pqq_q_mom1, WithinAbs(pqq_mom1 * q_mom1, 1e-6)); //< check momentum moment (1/6 * (-4/3))
  REQUIRE_THAT(          pqq_q.truncated_moment(0.0), WithinAbs(0.0, 1e-5));    //< check quark number conserved
  REQUIRE_THAT((pqq_q + pgq_q).truncated_moment(1.0), WithinAbs(0.0, 1e-6));    //< check momentum conserved

  //------- compound operations ------------------------
  hoppet::grid_conv p2 = pqq;
  REQUIRE(p2.ptr() != nullptr);
  REQUIRE(p2.ptr() != pqq.ptr());
  REQUIRE(p2.grid().ptr() == pqq.grid().ptr());
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
  REQUIRE_THAT( (pgen * q) .truncated_moment(1.0), WithinAbs(2*pgq_q_mom1, 1e-6)); // pview should be
}