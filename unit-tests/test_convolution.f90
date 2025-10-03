module test_convolution
  use integrator
  use convolution
  use streamlined_interface
  use unit_tests
  implicit none

  !procedure(ig_func), bind(C, name="pqq_c") :: pqq_c

  !interface
  !  function pqq_c(x) bind(C, name="pqq_c") result(res)
  !    use iso_c_binding, only: c_double
  !    real(c_double), intent(in) :: x
  !    real(c_double) :: res
  !  end function pqq_c
  !end interface
  procedure(ig_func_c), bind(C, name="pqq_c") :: pqq_c

contains

  subroutine test_InitGridConv()
    real(dp) :: res
    real(dp) :: aa= 0.1_dp, bb= 0.8_dp
    res = ig_LinWeight(pqq_c, aa,bb, 1.0_dp, 1.0_dp, 1e-10_dp)
    call check_approx_eq("test_InitGridConv-pqq_c-int", res, pqq_int_res(aa,bb), 1e-10_dp)

  contains
  
    function pqq_int_res(a, b) result(rr)
      real(dp), intent(in) :: a, b
      real(dp) :: rr
      rr = 0.5_dp * (a - b) * (2.0_dp + a + b) + 2.0_dp * log((-1.0_dp + a)/(-1.0_dp + b))
    end function pqq_int_res
    
end subroutine test_InitGridConv

end module test_convolution