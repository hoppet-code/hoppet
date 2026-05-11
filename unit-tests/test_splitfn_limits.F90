module test_splitfn_limits  
  use unit_tests
  implicit none

contains
  !! a subroutine to make sure that the divergent splitting functions
  !! have the correct limits as y->0, even for extremely small y (where
  !! rounding errors could arise), which is important for some usage in
  !! PanScales
  subroutine test_splitfn_lims()
    use types
    use hoppet
    use splitting_functions
    use convolution_communicator
    real(dp) :: y, pii
    if (.not. do_test("test_splitfn_limits")) return

    y = 1e-40_dp

    cc_piece = cc_REAL
    pii = sf_Pqq(y)
    call check_approx_eq_0d("pqq REAL at y=1e-40", pii, +2*CF/y, 1.0e-12_dp, 1.0e-12_dp, tol_choice_or=.true.)

    cc_piece = cc_VIRT
    pii = sf_Pqq(y)
    call check_approx_eq_0d("pqq REAL at y=1e-40", pii, -2*CF/y, 1.0e-12_dp, 1.0e-12_dp, tol_choice_or=.true.)

    cc_piece = cc_REAL
    pii = sf_Pgg(y)
    call check_approx_eq_0d("pqq REAL at y=1e-40", pii, +2*CA/y, 1.0e-12_dp, 1.0e-12_dp, tol_choice_or=.true.)

    cc_piece = cc_VIRT
    pii = sf_Pgg(y)
    call check_approx_eq_0d("pqq REAL at y=1e-40", pii, -2*CA/y, 1.0e-12_dp, 1.0e-12_dp, tol_choice_or=.true.)

  end subroutine test_splitfn_lims
end module test_splitfn_limits  