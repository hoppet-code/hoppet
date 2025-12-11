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


//--  implementations that don't need to be in the .hh files ----------------------------

namespace hoppet {

std::vector<double> ostream_y_values;
void set_ostream_y_values(const std::vector<double> & y_values) {
  ostream_y_values = y_values;
}

/// @brief  stream output operator for grid_quant_view objects
///
/// Outputs the monotonically increasing y values and the corresponding
/// grid_quant_view values, one pair per line.
///
std::ostream & operator<<(std::ostream & os, const grid_quant_view & gq) {
  using namespace std;
  if (ostream_y_values.size() != 0) 
    throw std::runtime_error("operator<<(ostream&, grid_quant_view&): non-empty ostream_y_values not yet implemented");

  vector<int> unique_indices = gq.grid().monotonic_indices();
  vector<double> yvals = gq.grid().y_values();
  for (size_t i = 0; i < unique_indices.size(); ++i) {
    int idx = unique_indices[i];
    os << yvals[idx] << " " << gq[idx] << "\n";
  }
  return os;  
}

/// @brief  stream output operator for grid_quant_2d_view objects
///
/// Outputs the monotonically increasing y=ln(1/x) values and the corresponding
/// grid_quant_2d_view values, xf(x,:) across all flavours, one line per y point.
///
std::ostream & operator<<(std::ostream & os, const grid_quant_2d_view & gq) {
  using namespace std;
  if (ostream_y_values.size() == 0) {
    vector<int> unique_indices = gq.grid().monotonic_indices();
    vector<double> yvals = gq.grid().y_values();
    for (size_t i = 0; i < unique_indices.size(); ++i) {
      int idx = unique_indices[i];
      os << yvals[idx];
      for (size_t flav = 0; flav < gq.size_flv(); ++flav) {
        os << " " << gq(flav, idx);
      }
      os << "\n";
    }
  } else {
    for (double y : ostream_y_values) {
      os << y;
      // not an efficient way of writing things...
      for (std::size_t flav = 0; flav < gq.size_flv(); ++flav) {
        os << " " << gq[flav].at_y(y);
      }
      os << "\n";
    }
  }
  return os;  
}


}