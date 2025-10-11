
#include  "ome/AggQ.h"
#include  "ome/AgqQ.h"
#include  "ome/AQg.h"
#include  "ome/AqgQ.h"
#include  "ome/AQqPS.h"
#include  "ome/AqqQNSEven.h"
#include  "ome/AqqQNSOdd.h"
#include  "ome/AqqQPS.h"

#include "../hoppet.h"

namespace hoppet {
namespace ome {
  
}} // namespace hoppet::ome

namespace ome {
  typedef const rpd_distribution<ome_as_view<double>, ome_as_plus_view<double>, ome_as_const_view<double>> ome_rpd;  
}


extern "C" {


//const void * ome_AggQ       = & ome::AggQ       ;  
//const void * ome_AgqQ       = & ome::AgqQ       ; 
//const void * ome_AQg        = & ome::AQg        ; 
//const void * ome_AqgQ       = & ome::AqgQ       ; 
//const void * ome_AQqPS      = & ome::AQqPS      ;  
//const void * ome_AqqQNSEven = & ome::AqqQNSEven ; 
//const void * ome_AqqQNSOdd  = & ome::AqqQNSOdd  ;  
//const void * ome_AqqQPS     = & ome::AqqQPS     ; 

const ome::ome_rpd * ome_AggQ       = & ome::AggQ       ;  
const ome::ome_rpd * ome_AgqQ       = & ome::AgqQ       ; 
const ome::ome_rpd * ome_AQg        = & ome::AQg        ; 
const ome::ome_rpd * ome_AqgQ       = & ome::AqgQ       ; 
const ome::ome_rpd * ome_AQqPS      = & ome::AQqPS      ;  
const ome::ome_rpd * ome_AqqQNSEven = & ome::AqqQNSEven ; 
const ome::ome_rpd * ome_AqqQNSOdd  = & ome::AqqQNSOdd  ;  
const ome::ome_rpd * ome_AqqQPS     = & ome::AqqQPS     ; 


/// @brief Evaluate a piece of an OME in a way suitable for HOPPET's standard integrator
/// @param ptr    pointer to the OME object 
/// @param y      value of y = ln(1/x)
/// @param piece  piece of the OME to evaluate (hoppet::cc_REAL, hoppet::cc_VIRT, hoppet::cc_REALVIRT, hoppet::cc_DELTA)
/// @param order  order of the OME to evaluate
/// @param LM     value of LM = log(m^2/mu^2) (double)
/// @param NF     value of NF (double)
/// @return       result of the evaluation
//double ome_piece_hoppet(void * ptr, 
double ome_piece_hoppet(const ome::ome_rpd & rpd, 
                        const double & y, 
                        const int & piece, 
                        const int & order, 
                        const double & LM, 
                        const double & NF) {
  //std::cout << "ome_piece_hoppet called with ptr=" << *ptr << ", y=" << y << ", piece=" << piece << ", order=" << order << ", LM=" << LM << ", NF=" << NF << std::endl;
  //const auto & rpd = *static_cast<const ome::rpd_distribution<ome::ome_as_view<double>, ome::ome_as_plus_view<double>, ome::ome_as_const_view<double>>*>(ptr);
  //const auto & rpd = ome::AqqQNSEven;
  //double NF = iNF;

  double x = exp(-y);

  //auto regular = [&]()->double {return rpd.has_regular() ? rpd.get_regular().value()(order, LM, NF, x) : 0.0;};
  //auto plus    = [&]()->double {return rpd.has_plus()    ? rpd.get_plus().value()(order, LM, NF, x) : 0.0;};
  //auto delta   = [&]()->double {return rpd.has_delta()   ? rpd.get_delta().value()(order, LM, NF) : 0.0;};

  auto regular = [&]()->double {return rpd.has_regular() ? rpd.get_regular().value()[order](LM, NF, x) : 0.0;};
  auto plus    = [&]()->double {return rpd.has_plus()    ? rpd.get_plus()   .value()[order](LM, NF, x) : 0.0;};
  auto delta   = [&]()->double {return rpd.has_delta()   ? rpd.get_delta()  .value()[order](LM, NF) : 0.0;};



  //auto regular = [&]()->double {return rpd.has_regular() ? rpd.get_regular().value()[order][0](NF, x) : 0.0;};
  //auto plus    = [&]()->double {return rpd.has_plus()    ? rpd.get_plus()   .value()[order][0](NF, x) : 0.0;};
  //auto delta   = [&]()->double {return rpd.has_delta()   ? rpd.get_delta()  .value()[order][0](NF) : 0.0;};

  // print out the arguments for debugging
  //std::cout << "ome_piece called with piece=" << piece << ", order=" << order << ", LM=" << LM << ", NF=" << NF << ", x=" << x << std::endl;
  //std::cout << "  rpd has parts: regular=" << rpd.has_regular() << ", plus=" << rpd.has_plus() << ", delta=" << rpd.has_delta() << std::endl;
  //std::cout << "  rpd regular value: " << regular() << std::endl;
  //std::cout << "  rpd plus value: " << plus() << std::endl;
  //std::cout << "  rpd delta value: " << delta() << std::endl;
  //std::cout << "  ome_AqqQNSEven_reg_coeff_as: " << ome_AqqQNSEven_reg_coeff_as(order, LM, NF, x) << std::endl;

  //double sum = 0.0;
  //for (int trunc_order = 0; trunc_order <= order; ++trunc_order) {
  //  double as = 0.1;
  //  //sum += pow(as, trunc_order) * ome_AqqQNSEven_reg_coeff_as_LM_NF(order, LM, NF, x);
  //  sum += pow(as, trunc_order) * rpd.get_regular().value()[trunc_order](LM, NF, x);
  //  //std::cout << " reg,  trunc_order=" << trunc_order 
  //  //          << ", sum = " << sum 
  //  //          << ", ome_trunc = " << ome_AqqQNSEven_reg_trunc_as(trunc_order, as, LM, NF, x)
  //  //          << ", coeff = " << ome_AqqQNSEven_reg_coeff_as(trunc_order, LM, NF, x)
  //  //          << std::endl;
  //  
  //}

  double result = 0.0;
  if        (piece == hoppet::cc_REAL) {
    if (rpd.has_regular()) result  = regular();
    if (rpd.has_plus())    result += plus();

  } else if (piece == hoppet::cc_VIRT) {
    if (rpd.has_regular()) result  = -plus();

  } else if (piece == hoppet::cc_REALVIRT) {
    if (rpd.has_regular()) result  = regular();

  } else if (piece == hoppet::cc_DELTA) {
    if (rpd.has_regular()) result  = delta();

  } else {
    // invalid piece
    return 0.0;
  }

  //std::cout << "  result = " << result << std::endl;
  if (piece != hoppet::cc_DELTA) {
    result *= x;
  }

  // hoppet wants a normalisation to (alpha_s/2pi)^order
  result /= pow(2,order);
  //std::cout << "  result (hoppet norm) = " << result << std::endl;
  return result;
}

} // extern "C"