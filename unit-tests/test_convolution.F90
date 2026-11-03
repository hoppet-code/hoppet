#include "inc/ftlMacros.inc"
module test_convolution
  use integrator
  use convolution
  use streamlined_interface
  use unit_tests  
  implicit none

  ! common objects that will be allocated once and then
  ! reused in the various tests
  real(dp), pointer :: xq(:), xq_conv1(:), xq_conv2(:)

  ! interface to the C functions we want to test
  abstract interface
    function ig_func_c(x) result(func) bind(C)
      use iso_c_binding, only: c_double
      real(c_double), intent(in), value :: x
      real(c_double)             :: func
    end function ig_func_c
  end interface
  public :: ig_func_c

  procedure(ig_func_c), bind(C, name="pqq_c"  ) :: pqq_c
  procedure(ig_func_c), bind(C, name="xpqqy_c") :: xpqqy_c
  
  ! a wrapper that converts from dp to c_double and back;
  ! on most systems this is a redundant copy, but it is the
  ! recommended way to work when interfacing with C
  !
  ! (actual use is below in the "contains" section)
#define MAKE_WRAPPER(NAME) \
  function CAT(NAME,_wrapper)(x) result(func); \
    use iso_c_binding, only: c_double; \
    real(dp), intent(in) :: x; \
    real(dp) :: func; \
    func = real(CAT(NAME,_c)(real(x,c_double)),dp); \
  end function


contains

  MAKE_WRAPPER(pqq)
  MAKE_WRAPPER(xpqqy)


  subroutine test_InitGridConv()
    real(dp) :: res

    ! make sure we have a few PDF objects
    call AllocGridQuant(grid, xq)
    call AllocGridQuant(grid, xq_conv1)
    call AllocGridQuant(grid, xq_conv2)
    xq = ff(yValues(grid))
    
    call test_mvv_interfaces()
    call test_approxDeltaFn()

    if (.not. do_test("test_InitGridConv")) return

    ! first a test that the integrating a simple function works fine
    block
      real(dp) :: aa= 0.1_dp, bb= 0.8_dp
      res = ig_LinWeight(pqq_wrapper, aa,bb, 1.0_dp, 1.0_dp, 1e-10_dp)
      call check_approx_eq("test_InitGridConv-pqq_c-int", res, pqq_int_res(aa,bb), 1e-10_dp)
    end block

    ! now a test of InitGridConv with a simple C++ function for Pqq,
    ! comparing to the standard LO result that we have in the DGLAP holder
    block
      type(grid_conv) :: pqq_gc
      call InitGridConv(grid, pqq_gc, xpqqy_wrapper)
      ! set up some simple x distribution on the grid
      ! convolute with P_qq from the DGLAP holder and from the C++ version
      xq_conv1 = dh%P_LO%NS_plus * xq
      xq_conv2 = pqq_gc * xq
      ! check the results agree
      call check_approx_eq_1d("test_InitGridConv-conv", xq_conv1, xq_conv2, 1e-10_dp)
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


  !! some basic tests that the approxDeltaFn routine works as expected
  !! 
  !! The core test will be that the convolution of pqq with the approx delta function
  !! gives the same result as the original pqq function, for a few test values of
  !! y. The test values are chosen to be sensitive to the four separate grid pieces
  subroutine test_approxDeltaFn()
    use splitting_functions
    use convolution_communicator
    type(grid_def) :: unlocked_grid, gd(4), this_grid, test_grids(2)
    real(dp) :: dy=0.1_dp, ymax = 10.0_dp    
    integer  :: order = -5, i, igrid = 2
    type(grid_conv) :: pqq
    real(dp), allocatable :: delta(:), pqq_delta(:), ytest_vals(:)
    real(dp) ::   y
    call InitGridDef(gd(4),dy/27.0_dp,0.2_dp, order=order)
    call InitGridDef(gd(3),dy/9.0_dp,0.5_dp,  order=order)
    call InitGridDef(gd(2),dy/3.0_dp,2.0_dp,  order=order)
    call InitGridDef(gd(1),dy,       ymax  ,  order=order)
    call InitGridDef(unlocked_grid,gd(1:4),locked=.false.)

    ! try out a single grid and a nested grid, to make sure the
    ! approxDeltaFn works in both cases (different paths in the code)
    test_grids = [unlocked_grid, gd(1)]
    ytest_vals = [0.18_dp, 0.48_dp, 1.95_dp, 3.0_dp] ! values just inside the grid edges
    do igrid = 1, size(test_grids)
      this_grid = test_grids(igrid)
      call InitGridConv(this_grid, pqq, sf_Pqq)
      delta = approxDeltaFn(this_grid)
      pqq_delta = pqq * delta
      do i = 1, size(ytest_vals)
        if (igrid == 2 .and. i < 4) cycle ! only test the last value for the coarse grid
        y = ytest_vals(i)
        call check_approx_eq("test_approxDeltaFn", pqq_delta.aty.(y.with.this_grid), sf_Pqq(y), 1e-7_dp)
      end do
      call Delete(pqq)
    end do

    call Delete(unlocked_grid)
  end subroutine test_approxDeltaFn


  !! Routines to test the mvv_splitting_function type
  !! against their original implementation 
  subroutine test_mvv_interfaces()
    use hoppet_splitting_function_interfaces
    use splitting_functions_nnlo_p
    use convolution_communicator
    use xpij2p
    type(grid_conv) :: pgg_orig, pqg_orig
    type(grid_conv) :: pgg_OO,  pqg_OO
    real(dp), pointer :: xq1(:), xq2(:)
    real(dp) :: res_mvv, res_orig
    integer  :: piece
    type(mvv_splitting_function) :: pgg_mvv, pqg_mvv

    if (.not. do_test("test_mvv_interfaces")) return

    pgg_mvv = mvv_splitting_function(P2GGA, P2GGB , P2GGC , 0.5_dp**3)
    pqg_mvv = mvv_splitting_function(P2QGA, null(), null(), 0.5_dp**3)

    do piece = 1, 4
      cc_piece = piece
      res_mvv = pgg_mvv%f(0.5_dp, piece)
      res_orig = sf_P2gg(0.5_dp) 
      !write(6,*) "test_mvv_interfaces-pgg: piece=", piece, " res_mvv=", res_mvv, " res_orig=", res_orig
      call check_approx_eq("test_mvv_interfaces-pgg", res_mvv, res_orig, 1e-10_dp)

      res_mvv = pqg_mvv%f(0.5_dp, piece)
      res_orig = sf_P2qg2nf(0.5_dp) 
      !write(6,*) "test_mvv_interfaces-pqg: piece=", piece, " res_mvv=", res_mvv, " res_orig=", res_orig
      call check_approx_eq("test_mvv_interfaces-pqg", res_mvv, res_orig, 1e-10_dp)
    end do

    call InitGridConv(grid, pgg_orig, sf_P2gg)
    call InitGridConv(grid, pqg_orig, sf_P2qg2nf)

    call InitGridConv(grid, pgg_OO, mvv_splitting_function(P2GGA, P2GGB , P2GGC , 0.5_dp**3))
    call InitGridConv(grid, pqg_OO, mvv_splitting_function(P2QGA, null(), null(), 0.5_dp**3))

    xq_conv1 = pgg_orig * xq
    xq_conv2 = pgg_OO   * xq
    call check_approx_eq_1d("test_mvv_interfaces-pgg", xq_conv2, xq_conv1, 1e-10_dp)

    xq_conv1 = pqg_orig * xq
    xq_conv2 = pqg_OO   * xq
    call check_approx_eq_1d("test_mvv_interfaces-pqg", xq_conv2, xq_conv1, 1e-10_dp)

    call Delete(pgg_orig); call Delete(pqg_orig)
    call Delete(pgg_OO);   call Delete(pqg_OO)
  end subroutine test_mvv_interfaces

end module test_convolution