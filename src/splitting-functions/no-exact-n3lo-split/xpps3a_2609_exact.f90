! Dummy module to use when we are not compiling the exact n3lo
! splitting functions
module xpps3a_2609_exact
  use warnings_and_errors
  character(len=*), parameter :: name_xpps3a = "xpps3a_2609_exact"
contains
  FUNCTION P3PSA_2609_exact (Y, NF)
    !
    IMPLICIT REAL*8 (A-Z)
    INTEGER nf

    call wae_error('P3PSA_2609_exact: Exact n3lo splitting functions not compiled in.')

    
    RETURN
  END FUNCTION P3PSA_2609_exact
  !
  ! =================================================================av==
end module xpps3a_2609_exact
