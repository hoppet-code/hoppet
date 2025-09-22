module daind_fake
  use types
  implicit none

contains

  !************************************************************************
  ! The following are the interface parameters from the original daind.f file
  ! provided to the hoppet authors by Johannes Blumlein, which itself
  ! quotes CF. Piessens, R., An Algorithm for Automatic Integration. 
  ! Angewandte Informatik, 9 (1973), 399--401.
  !
  !  INPUTPARAMETERS
  !  A,B      LIMITS OF THE INTEGRATION INTERVAL
  !  FUN      FUNCTION TO BE INTEGRATED (TO BE DECLARED EXTERNAL IN THE MAIN PR.)
  !  EPS      ABSOLUTE OR RELATIVE TOLERANCE,DEPENDING OF THE VALUE OF 'KEY'
  !  KEY      =1 THEN 'EPS' DENOTES AN ABSOLUTE, =2 THEN A RELATIVE TOLERANCE
  !  MAX      UPPER BOUND ON THE NUMBERS OF INTEGRAND EVALUATIONS (MAX.LE.10000)
  !
  !  OUTPUTPARAMETERS
  !  KOUNT    NUMBER OF INTEGRAND EVALUATIONS
  !  EST      ESTIMATION OF THE ABSOLUTE ERROR OF THE APPROXIMATION
  !
  ! Here, the interface just maps to hoppet's ig_LinWeight (CERNLIB dgauss variant)
  ! with the following adaptations relative to the original DAIND
  ! 
  !  - EPS is treated as in dgauss, i.e. as an absolute tolerance if |result| < 1
  !    relative otherwise
  !  - KEY and MAX are ignored
  !  - KOUNT is returned as 0
  !  - EST is returned as 0.0_dp
  function daind(a, b, fun, eps, key, max, kount, est)
    use integrator
    real(dp) :: daind
    real(dp), intent(in) :: a, b, eps
    external fun
    real(dp) :: fun
    !interface 
    !  function fun(x)
    !    use types; implicit none
    !    real(dp) :: x
    !    real(dp) :: fun
    !  end function fun
    !end interface
    integer, intent(in) :: key, max
    integer, intent(out) :: kount
    real(dp), intent(out) :: est

    kount = 0
    est = 0.0_dp

    daind = ig_LinWeight(fun, a, b, 1.0_dp, 1.0_dp, eps)
  end function daind
  
end module daind_fake