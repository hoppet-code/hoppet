
#include  "ome/AggQ.h"
#include  "ome/AgqQ.h"
#include  "ome/AQg.h"
#include  "ome/AqgQ.h"
#include  "ome/AQqPS.h"
#include  "ome/AqqQNSEven.h"
#include  "ome/AqqQNSOdd.h"
#include  "ome/AqqQPS.h"

#include "../hoppet.h"

namespace ome {
  typedef const rpd_distribution<ome_as_view<double>, ome_as_plus_view<double>, ome_as_const_view<double>> ome_rpd;  
}


extern "C" {

  // C-style pointers to the ome:A[...] objects, intended for 
  // accessibility from the Fortran interface
  const ome::ome_rpd * ome_AggQ       = & ome::AggQ       ;  
  const ome::ome_rpd * ome_AgqQ       = & ome::AgqQ       ; 
  const ome::ome_rpd * ome_AQg        = & ome::AQg        ; 
  const ome::ome_rpd * ome_AqgQ       = & ome::AqgQ       ; 
  const ome::ome_rpd * ome_AQqPS      = & ome::AQqPS      ;  
  const ome::ome_rpd * ome_AqqQNSEven = & ome::AqqQNSEven ; 
  const ome::ome_rpd * ome_AqqQNSOdd  = & ome::AqqQNSOdd  ;  
  const ome::ome_rpd * ome_AqqQPS     = & ome::AqqQPS     ; 
}

/// @brief Get the name of an OME object
/// @param rpd    reference to the OME object
/// @return       name of the OME object
std::string ome_name(const ome::ome_rpd & rpd) {
  if      (&rpd == ome_AggQ      ) return "AggQ";
  else if (&rpd == ome_AgqQ      ) return "AgqQ";
  else if (&rpd == ome_AQg       ) return "AQg";
  else if (&rpd == ome_AqgQ      ) return "AqgQ";
  else if (&rpd == ome_AQqPS     ) return "AQqPS";
  else if (&rpd == ome_AqqQNSEven) return "AqqQNSEven";
  else if (&rpd == ome_AqqQNSOdd ) return "AqqQNSOdd";
  else if (&rpd == ome_AqqQPS    ) return "AqqQPS";
  else return "unknown ome_rpd pointer";
}

extern "C" {

/// @brief Evaluate a piece of an OME in a way suitable for HOPPET's standard integrator
/// @param rpd    reference to the OME object 
/// @param y      value of y = ln(1/x)
/// @param piece  piece of the OME to evaluate (hoppet::cc_REAL, hoppet::cc_VIRT, hoppet::cc_REALVIRT, hoppet::cc_DELTA)
/// @param order  order of the OME to evaluate
/// @param LM     value of LM = log(m^2/mu^2) (double)
/// @param NF     value of NF (double)
/// @return       result of the evaluation
double ome_piece_hoppet(const ome::ome_rpd & rpd, 
                        const double & y, 
                        const int & piece, 
                        const int & order, 
                        const double & LM, 
                        const double & NF) {
  //std::cout << "ome_piece_hoppet called with ptr=" << ome_name(rpd) << ", y=" << y << ", piece=" << piece << ", order=" << order << ", LM=" << LM << ", NF=" << NF << std::endl;

  double x = exp(-y);

  double result = 0.0;
  if (LM == 0.0) {
    auto regular = [&]()->double {return rpd.has_regular() ? rpd.get_regular().value()[order][LM](NF, x) : 0.0;};
    auto plus    = [&]()->double {return rpd.has_plus()    ? rpd.get_plus()   .value()[order][LM](NF, x) : 0.0;};
    auto delta   = [&]()->double {return rpd.has_delta()   ? rpd.get_delta()  .value()[order][LM](NF   ) : 0.0;};
  
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
  } else {
    auto regular = [&]()->double {return rpd.has_regular() ? rpd.get_regular().value()[order](LM, NF, x) : 0.0;};
    auto plus    = [&]()->double {return rpd.has_plus()    ? rpd.get_plus()   .value()[order](LM, NF, x) : 0.0;};
    auto delta   = [&]()->double {return rpd.has_delta()   ? rpd.get_delta()  .value()[order](LM, NF   ) : 0.0;};

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