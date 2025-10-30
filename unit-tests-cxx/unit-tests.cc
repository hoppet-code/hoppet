#include "hoppet_oo.h"
#include <iostream>

using namespace std;

void example() {
  hoppet::grid_def grid(0.1, 10.0);
  std::cout << "grid.ny(): "  << grid.ny() << std::endl;
  std::cout << "grid.ptr(): " << grid.ptr() << std::endl;
  std::vector<double> xvals = grid.x_values();
  std::vector<double> yvals = grid.y_values();
  std::cout << "x values: ";
  for (double x : xvals) {
    std::cout << x << " ";
  }
  std::cout << std::endl;

  hoppet::grid_quant q(grid);
  std::cout << "grid_quant.ptr(): " << q.ptr() << std::endl;
  q = [](double y) { return y*y; };
  //q = [](double y) { double x = exp(-y);return pow(1-x,4); };

  hoppet::grid_quant q2 = q; // copy constructor
  cout << "q.at_y(5.0): " << q.at_y(5.0) << endl;

  int iy = grid.ny() / 2;
  cout << "q[" << iy << "] = " << q[iy] << " " << pow(yvals[iy],2) <<  endl;

}

int main() {
  example();
  return 0;
}