! Dummy module to use when we are not compiling the exact n3lo
! splitting functions
module xpns3m_2604_exact
  use warnings_and_errors
  character(len=*), parameter :: name_xpns3m = "xpns3m_2604_exact"
contains

  FUNCTION P3NSMA_2604_exact (Y, NF)
    !
    IMPLICIT REAL*8 (A-Z)
    INTEGER nf

    call wae_error('P3NSMA_2604_exact: Exact n3lo splitting functions not compiled in.')
    
    RETURN
  END FUNCTION P3NSMA_2604_EXACT

  FUNCTION P3NSMB_2604_exact (Y, NF)
    !
    IMPLICIT REAL*8 (A-Z)
    INTEGER nf

    call wae_error('P3NSMB_2604_exact: Exact n3lo splitting functions not compiled in.')
    
    RETURN
  END FUNCTION P3NSMB_2604_EXACT

  FUNCTION P3NSMC_2604_exact (Y, NF)
    !
    IMPLICIT REAL*8 (A-Z)
    INTEGER nf

    call wae_error('P3NSMC_2604_exact: Exact n3lo splitting functions not compiled in.')
    
    RETURN
  END FUNCTION P3NSMC_2604_EXACT

end module xpns3m_2604_exact
