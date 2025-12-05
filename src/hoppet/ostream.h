#ifndef __HOPPET_OSTREAM__
#define __HOPPET_OSTREAM__

#include <ostream>
#include "hoppet_oo.h"

namespace hoppet {
    
/// @brief  stream output operator for grid_quant_view objects
///
/// Outputs the monotonically increasing y=ln(1/x) values and the corresponding
/// grid_quant_view values, xf(x), one pair per line.
///
std::ostream & operator<<(std::ostream & os, const grid_quant_view & gq);    


/// @brief  stream output operator for grid_quant_2d_view objects
///
/// Outputs the monotonically increasing y=ln(1/x) values and the corresponding
/// grid_quant_2d_view values, xf(x,:) across all flavours, one line per y point.
///
/// Note that QCD PDFs typically have an extra flavour index that
/// encodes the representation (human v. evolution basis); this gets
/// included in the output.
std::ostream & operator<<(std::ostream & os, const grid_quant_2d_view & gq);


} // namespace hoppet

#endif // __HOPPET_OSTREAM__