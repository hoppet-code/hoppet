#include "hoppet_oo.h"
#include <iostream>

void example() {
  hoppet::grid_def grid(0.1, 10.0);
  std::cout << "grid.ny(): "  << grid.ny() << std::endl;
  std::cout << "grid.ptr(): " << grid.ptr() << std::endl;
  std::vector<double> xvals = grid.x_values();
  std::cout << "x values: ";
  for (double x : xvals) {
    std::cout << x << " ";
  }
  std::cout << std::endl;
}

int main() {
  example();
  return 0;
}