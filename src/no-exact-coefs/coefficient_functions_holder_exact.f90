!=============================================================================
! dummy module that gets included if exact DIS coefficient functions are not
! compiled in
module coefficient_functions_holder_exact
  implicit none
end module coefficient_functions_holder_exact

!=============================================================================
! This just exits the program with an error message
subroutine InitCoefHolderExact(grid, ch)
  use types; use coefficient_functions_holder, dummy => InitCoefHolderExact
  use warnings_and_errors
  implicit none
  type(grid_def), intent(in) :: grid
  type(coef_holder), intent(inout), target :: ch

  call wae_error('InitCoefHolderExact: exact DIS coefficient functions not compiled in.')
end subroutine InitCoefHolderExact