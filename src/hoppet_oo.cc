#include "hoppet_oo.h"



extern "C" double hoppet_grid_conv_f__wrapper(double y, int piece, void*ctx) {
  auto func = static_cast<std::function<double(const double, const int)>*>(ctx);
  return (*func)(y, piece);
}

namespace hoppet {
namespace sl {
  grid_def_view grid;
  dglap_holder_view dh;
}
}

extern "C" {
void hoppetStartCXX(grid_def_f * grid_ptr, dglap_holder_f * dh_ptr) {
  hoppet::sl::grid = hoppet::grid_def_view(grid_ptr);
  hoppet::sl::dh = hoppet::dglap_holder_view(dh_ptr);
}
} // extern "C"