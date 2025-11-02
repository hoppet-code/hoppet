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

/// We supply the main routine, to make sure that we can 
/// do any global-object initialisation that's needed
int main(int argc, char* argv[]) {
  // Global setup
  // make a simple grid for use in the tests
  grid100 = hoppet::grid_def(0.1, 10.0);

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


TEST_CASE( "hoppet grid_def and grid_quant basic functionality", "[hoppet]" ) {
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
  REQUIRE_THAT( grid100.y_values()[iy], WithinRel( grid2.y_values()[iy], 1e-12) );

  hoppet::grid_def grid3({grid100, hoppet::grid_def(0.03,0.3)}, true);
  hoppet::grid_def grid4;
  grid4 = grid3;

  auto grid5 = hoppet::grid_def_default(0.2, 12.0, -5);
}


TEST_CASE( "hoppet grid_quant basic functionality", "[hoppet]" ) {

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
}

TEST_CASE( "hoppet streamlined objects", "[hoppet]" ) {
  hoppet::grid_quant q(hoppet::sl::grid);
  //hoppet::grid_quant q(hoppet_sl_grid_ptr);
  q = [](double y) { return y*y; };
  REQUIRE_THAT( q.at_y(5.0), WithinAbs(25.0, 1e-6));
}

TEST_CASE( "hoppet grid_conv basic functionality", "[hoppet]" ) {
  hoppet::grid_conv conv(grid100, 
    [](double y, int piece) {
      double x = exp(-y);
      double res;
      if      (piece == hoppet::cc_REAL) res =  (1+x*x)/(1-x);
      else if (piece == hoppet::cc_VIRT) res = -(1+x*x)/(1-x);
      else                               res = 0.0;
      return x*res;
    });

  auto power = [](double powval) { return grid100 * [powval](double y) { return exp(-powval*y); }; };

  //hoppet::grid_quant q = grid100 * power(1.0);
  hoppet::grid_quant result = conv * power(1.0);
  cout << "After convolution, result.at_y(10.0): " << result.at_y(10.0) << endl;
  //q = power(0.0);
  result = conv * power(0.0);
  cout << "After convolution, result.at_y(10.0): " << result.at_y(10.0) << endl;
  result = conv * power(4.0);
  cout << "After convolution, result.at_y(10.0): " << result.at_y(10.0) << endl;
}