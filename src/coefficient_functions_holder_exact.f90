
module coefficient_functions_holder_exact
  use types; use warnings_and_errors
  use consts_dp; use convolution_communicator
  use qcd
  use convolution
  use dglap_objects
  use coefficient_functions_holder 
  use coefficient_functions_holder_internal
  implicit none

contains

  subroutine InitCoefHolderExact_impl(grid, ch, nloop, nflcl)
    type(grid_def), intent(in) :: grid
    type(coef_holder), intent(inout), target :: ch
    integer, intent(in) :: nloop, nflcl
    !------------------------------

  end subroutine InitCoefHolderExact_impl

end module coefficient_functions_holder_exact


!=============================================================================
! This just calls the function in coefficient_functions_holder_exact
subroutine InitCoefHolderExact(grid, ch, nloop, nflcl)
  use types; use coefficient_functions_holder_exact, dummy => InitCoefHolderExact
  implicit none
  type(grid_def), intent(in) :: grid
  type(coef_holder), intent(inout), target :: ch
  integer, intent(in) :: nloop, nflcl
  call InitCoefHolderExact_impl(grid, ch, nloop, nflcl)
end subroutine InitCoefHolderExact