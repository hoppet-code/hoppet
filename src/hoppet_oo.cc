#include "hoppet_oo.h"
#include <cstring>


extern "C" void hoppet_throw_runtime_error(char * msg) {
  std::string msg_str(msg ? msg : "hoppet_throw_runtime_error: unknown error");
  // we have ownership of msg, so now that we've copied it into a string, we must delete it
  if (msg) delete[] msg;
  //std::cout << "hoppet_throw_runtime_error: throwing runtime_error with message: " << msg_str << std::endl;
  throw std::runtime_error(msg_str);
}

extern "C" void hoppet_print_cloc(void * ptr) {
  std::cout << "hoppet_print_cloc: pointer value = " << ptr << std::endl;
}

/// @brief  A wrapper to call a C++ convolution function object from Fortran
/// @param y      the y-value at which to evaluate the function
/// @param piece  the piece index (cc_REAL, etc.)
/// @param ctx    a pointer to the std::function<double(double,int)> object
/// @return the value of the function
///
/// This is used by the Fortran convolution routines to call back into C++ code
/// so that we can use arbitrary C++ function objects as convolution kernels
///
extern "C" double hoppet_grid_conv_f__wrapper(double y, int piece, void*ctx) {
  auto func = static_cast<std::function<double(const double, const int)>*>(ctx);
  return (*func)(y, piece);
}

/// @brief Allocate a C-style string of given size (for use in fortran -> C++ string conversion)
/// @param size 
/// @return a pointer to the allocated char array
extern "C" char * hoppet_allocate_cstr(int size) {
  return new char[size];
}

/// @brief Return the length of a C-style string (up to but not including the null terminator)
/// @param cstr_ptr 
/// @return 
extern "C" int hoppet_cstr_len(char * cstr_ptr) {return strlen(cstr_ptr);}

namespace hoppet {
namespace sl {
  grid_def_view        grid;
  dglap_holder_view    dh;
  pdf_table_view       table;
  pdf_table_array_view tables;
  running_coupling_view coupling; 
}
}

/// @brief functions to enable the C++ streamlined interface to view the underlying Fortran objects
extern "C" {

  /// @brief set up the C++ streamlined interface grid and dglap_holder views  
  void hoppet_sl_register_objects(grid_def_f * grid_ptr, dglap_holder_f * dh_ptr) {
    hoppet::sl::grid = hoppet::grid_def_view(grid_ptr);
    hoppet::sl::dh.take_view(dh_ptr);
  }

  void hoppet_sl_register_coupling(running_coupling_f * coupling_ptr) {
    hoppet::sl::coupling.take_view(coupling_ptr);
  }

  void hoppet_sl_register_tables(pdf_table_f * tab_ptr, int n_tables) {
    hoppet::sl::table.take_view(tab_ptr);
    hoppet::sl::tables.take_view(hoppet::pdf_table_array_view(tab_ptr, n_tables));
  }

//  /// @brief set the pointers to the objects that contain the evolved information
//  void hoppetSetEvolvedPointers(running_coupling_f * coupling_ptr, pdf_table_f * tab_ptr, int n_tables) {
//    if (coupling_ptr) hoppet::sl::coupling.take_view(coupling_ptr);
//    if (tab_ptr) {
//      hoppet::sl::table.take_view(tab_ptr);
//      hoppet::sl::tables.take_view(hoppet::pdf_table_array_view(tab_ptr, n_tables));
//    }
//  }
} // extern "C"