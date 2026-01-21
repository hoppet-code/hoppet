
#include "inc/ftlMacros.inc"
#include "inc/extraMacros.inc"


! a macro to generate the splitting function bind(C) interface functions
#define DEFINE_SPLIT_FN(NAME) \
  real(c_double) function CAT(hoppet_cxx__splitfn__,NAME)(y, piece) result(res) bind(C) ;\
    use splitting_functions; use convolution_communicator ;\
    real(c_double), intent(in), value :: y ;\
    integer(c_int), intent(in), value :: piece ;\
    cc_piece = piece ;\
    res = CAT(sf_,NAME)(y) ;\
  end function


!! Module that contains bind(C) interfaces for the splitting function routines.
!!
!! They will all have a y,piece signature, and the piece argument will be used to
!! set the global cc_piece variable used internally by the splitting function routines.
module hoppet_cpp_split_fns
  use iso_c_binding
  implicit none
contains

  DEFINE_SPLIT_FN(pqq)
  DEFINE_SPLIT_FN(pqg)
  DEFINE_SPLIT_FN(pgq)
  DEFINE_SPLIT_FN(pgg)

  DEFINE_SPLIT_FN(p1qqv)
  DEFINE_SPLIT_FN(p1qqbarv)
  DEFINE_SPLIT_FN(p1qqs)
  DEFINE_SPLIT_FN(p1qg)
  DEFINE_SPLIT_FN(p1gq)
  DEFINE_SPLIT_FN(p1gg)

  DEFINE_SPLIT_FN(p2gg)
  DEFINE_SPLIT_FN(p2qg2nf)
  DEFINE_SPLIT_FN(p2ps)
  DEFINE_SPLIT_FN(p2gq)
  DEFINE_SPLIT_FN(p2nsplus)
  DEFINE_SPLIT_FN(p2nsminus)
  DEFINE_SPLIT_FN(p2nss)

  DEFINE_SPLIT_FN(p3gg)
  DEFINE_SPLIT_FN(p3qg2nf)
  DEFINE_SPLIT_FN(p3ps)
  DEFINE_SPLIT_FN(p3gq)
  DEFINE_SPLIT_FN(p3nsplus)
  DEFINE_SPLIT_FN(p3nsminus)
  DEFINE_SPLIT_FN(p3nss)

end module hoppet_cpp_split_fns