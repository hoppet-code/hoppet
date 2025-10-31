#include "hoppet_oo.h"
#include <iostream>
#include <cmath>

//#include "catch_amalgamated.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
// define a shorthands for Catch::Matchers::WithinAbs etc.
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


using namespace std;

hoppet::grid_def grid;

TEST_CASE( "hoppet grid_def and grid_quant basic functionality", "[hoppet]" ) {
  grid = hoppet::grid_def(0.1, 10.0);
  std::cout << "grid.ny(): "  << grid.ny() << std::endl;
  std::cout << "grid.ptr(): " << grid.ptr() << std::endl;
  std::vector<double> xvals = grid.x_values();
  std::vector<double> yvals = grid.y_values();

  int iy = grid.ny() / 2;
  REQUIRE_THAT( yvals[iy], WithinRel(5.0, 1e-12) );
  REQUIRE_THAT( xvals[iy], WithinRel(exp(-5.0), 1e-12) );

  hoppet::grid_def grid2 = grid; // copy constructor
  REQUIRE( grid.ptr() != grid2.ptr() ); // different underlying pointers
  REQUIRE_THAT( grid.y_values()[iy], WithinRel( grid2.y_values()[iy], 1e-12) );

  hoppet::grid_def grid3({grid, hoppet::grid_def(0.03,0.3)}, true);
  hoppet::grid_def grid4;
  grid4 = grid3;

  auto grid5 = hoppet::grid_def_default(0.2, 12.0, -5);
}


TEST_CASE( "hoppet grid_quant basic functionality", "[hoppet]" ) {

  hoppet::grid_quant q(grid);

  std::cout << "grid_quant.ptr(): " << q.ptr() << std::endl;
  q = [](double y) { return y*y; };

  // test the self-assignment operator
  q = q;


  // check value is OK (approx, in case we change the grid definition)
  REQUIRE_THAT( q.at_y(5.0), WithinAbs(25.0, 1e-6));

  // make sure we can copy things OK and that we get a genuine copy
  int iy = grid.ny() / 2;
  hoppet::grid_quant q2 = q; // copy constructor
  REQUIRE ( q.ptr() != q2.ptr() ); // different underlying pointers
  REQUIRE ( q[iy] == q2[iy] );     // same values

  hoppet::grid_quant q3(grid);
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
}

