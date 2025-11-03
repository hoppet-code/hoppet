#include "hoppet_oo.h"



extern "C" double hoppet_grid_conv_f__wrapper(double y, int piece, void*ctx) {
  auto func = static_cast<std::function<double(const double, const int)>*>(ctx);
  return (*func)(y, piece);
}

namespace hoppet {
namespace sl {
  grid_def_view grid;
}
}

extern "C" {
void hoppetStartCXX(grid_def_f * grid_ptr) {
  hoppet::sl::grid = hoppet::grid_def_view(grid_ptr);
}
} // extern "C"