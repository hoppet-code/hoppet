#include "hoppet_oo.h"


namespace hoppet {
namespace sl {
  grid_def_view grid;
}
}

extern "C" {
void hoppetStartCXX() {
  hoppet::sl::grid = hoppet::grid_def_view(hoppet_sl_grid_ptr);
}
} // extern "C"