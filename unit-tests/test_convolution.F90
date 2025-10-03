#include "inc/ftlMacros.inc"
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
  !procedure(ig_func), bind(C, name="pqq_c") :: pqq_c
  !procedure(ig_func), bind(C, name="xpqqy_c") :: xpqqy_c

!  interface; \
!    function CAT(NAME,_c)(x) result(func) bind(C, name=CAT(NAME,_c)) \
!      use iso_c_binding; real(c_double), value :: x; real(c_double) :: func; \
!    end function; \
!  end interface; 
  procedure(ig_func_c), bind(C, name="pqq_c") :: pqq_c
  procedure(ig_func_c), bind(C, name="xpqqy_c") :: xpqqy_c
  
#define MAKE_WRAPPER(NAME) \
  function CAT(NAME,_wrapper)(x) result(func); \
    use iso_c_binding, only: c_double; \
    real(dp), intent(in) :: x; real(dp) :: func; \
    func = real(CAT(NAME,_c)(real(x,c_double)),dp); \
  end function

  !MAKE_WRAPPER(xpqqy)

contains

  subroutine test_InitGridConv()
    real(dp) :: res

    ! first a test that the integrating a simple function works fine
    block
      real(dp) :: aa= 0.1_dp, bb= 0.8_dp
      res = ig_LinWeight(pqq_wrapper, aa,bb, 1.0_dp, 1.0_dp, 1e-10_dp)
      call check_approx_eq("test_InitGridConv-pqq_c-int", res, pqq_int_res(aa,bb), 1e-10_dp)
    end block

    block
      real(dp), pointer :: xq(:), xq_conv1(:), xq_conv2(:)
      type(grid_conv) :: pqq_gc
      call InitGridConv(grid, pqq_gc, xpqqy_wrapper)
      call AllocGridQuant(grid, xq)
      call AllocGridQuant(grid, xq_conv1)
      call AllocGridQuant(grid, xq_conv2)
      ! set up some simple x distribution on the grid
      xq = ff(yValues(grid))
      ! convolute with P_qq from the DGLAP holder and from the C++ version
      xq_conv1 = dh%P_LO%NS_plus * xq
      xq_conv2 = pqq_gc * xq
      ! check the results agree
      call check_approx_eq_1d("test_InitGridConv-conv", xq_conv1, xq_conv2, 1e-10_dp)
      deallocate(xq, xq_conv1, xq_conv2)
      call Delete(pqq_gc)
    end block

  contains
  
    function pqq_int_res(a, b) result(rr)
      real(dp), intent(in) :: a, b
      real(dp) :: rr
      rr = 0.5_dp * (a - b) * (2.0_dp + a + b) + 2.0_dp * log((-1.0_dp + a)/(-1.0_dp + b))
    end function pqq_int_res
    
    elemental function ff(y) result(res)
      real(dp), intent(in) :: y
      real(dp) :: res
      real(dp) :: x
      x = exp(-y)
      res = x**(-0.5_dp) * (1.0_dp - x)**3
    end function ff

  end subroutine test_InitGridConv

  MAKE_WRAPPER(pqq)
  MAKE_WRAPPER(xpqqy)


end module test_convolution