! Dummy module to use when we are not compiling the exact n3lo
! splitting functions
module xpns3s_2604_exact
  use warnings_and_errors
  character(len=*), parameter :: name_xpns3s = "xpns3s_2604_exact"
contains

  FUNCTION P3NSSA_2604_exact (Y, NF)
    !
    IMPLICIT REAL*8 (A-Z)
    INTEGER nf

    call wae_error('P3NSSA_2604_exact: Exact n3lo splitting functions not compiled in.')
    
    RETURN
  END FUNCTION P3NSSA_2604_EXACT
  
end module xpns3s_2604_exact
