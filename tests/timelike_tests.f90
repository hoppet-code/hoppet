!! program to test coding of time-like splitting functions
!! against 3rd party (e.g. Apfel++)
!!
program timelike_tests
  use types; use consts_dp
  use convolution_communicator;
  use splitting_functions
  use qcd
  use sub_defs_io
  implicit none
  !---------------------
  integer  :: i
  integer  :: n=100
  real(dp) :: x, y, splitfn

  call qcd_SetNf(5)
  cc_piece = cc_REAL
  do i = 1, n-1
     x = (i-0.9_dp)*one/n
     y = -log(x)
     splitfn = 2*nf*sf_TP1gq(y)
     write(6,*) x, splitfn/x
  end do
  
end program timelike_tests


