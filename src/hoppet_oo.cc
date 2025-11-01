#include "hoppet_oo.h"


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