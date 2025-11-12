#include "hoppet_oo.h"


extern "C" void hoppet_throw_runtime_error() {throw std::runtime_error("hoppet wae_error");}

extern "C" double hoppet_grid_conv_f__wrapper(double y, int piece, void*ctx) {
  auto func = static_cast<std::function<double(const double, const int)>*>(ctx);
  return (*func)(y, piece);
}

namespace hoppet {
namespace sl {
  grid_def_view     grid;
  dglap_holder_view dh;
  pdf_table_view    table;
  running_coupling_view coupling; 
}
}

extern "C" {
/// @brief set up the C++ streamlined interface objects to view the Fortran ones  
void hoppetStartCXX(grid_def_f * grid_ptr, dglap_holder_f * dh_ptr) {
  hoppet::sl::grid = hoppet::grid_def_view(grid_ptr);
  hoppet::sl::dh.take_view(dh_ptr);
}
// maybe fix up the name here
void hoppetSetEvolvedPointers(running_coupling_f * coupling_ptr, pdf_table_f * tab_ptr) {
  hoppet::sl::coupling.take_view(coupling_ptr);
  hoppet::sl::table.take_view(tab_ptr);
}
} // extern "C"