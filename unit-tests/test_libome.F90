#include "inc/ftlMacros.inc"

module test_libome
  use unit_tests
contains
  !! subroutine for testing the hoppet-libome interface
  !! 
  !! reference values can be found in https://arxiv.org/abs/2510.02175, table 2
  subroutine test_libome_interface
    use hoppet_libome_interfaces
    use hoppet_term ! for term colours
    use convolution_pieces
    use convolution
    use streamlined_interface
    use iso_c_binding
    implicit none
    real(c_double) :: as4pi, LM, NF, x, y
    real(c_double) :: res_reg, res_plus, res_delta
    integer(c_int) :: piece
    integer, parameter :: reg = cc_REALVIRT
    integer, parameter :: delta = cc_DELTA
    integer, parameter :: plus = cc_VIRT 

    if (.not. do_test("test_libome_interface")) return

    ! Example values for testing
    as4pi = 0.1_dp
    LM    = 10._dp
    NF    = 3.0_dp
    x     = 0.1_dp

    ! check all 3rd order entries from table2 of https://arxiv.org/pdf/2510.02175
    call sum_check("AqqQNSEven_reg" ,    6.4686126472_dp , AqqQNSEven_ptr , reg, 3)
    call sum_check("AqqQNSEven_plus" ,  -4.7428175642_dp , AqqQNSEven_ptr , plus, 3)
    call sum_check("AqqQNSEven_delta" , -2.6422012017_dp , AqqQNSEven_ptr , delta, 3)

    call sum_check("AqqQNSOdd_reg" ,    6.6927772815_dp  , AqqQNSOdd_ptr , reg, 3)
    call sum_check("AqqQNSOdd_plus" ,  -4.7428175642_dp  , AqqQNSOdd_ptr , plus, 3)
    call sum_check("AqqQNSOdd_delta" , -2.6422012017_dp  , AqqQNSOdd_ptr , delta, 3)

    call sum_check("AggQ_reg" ,   -262.65714922_dp   , AggQ_ptr , reg, 3)
    call sum_check("AggQ_plus" ,    -6.6285411655_dp , AggQ_ptr , plus, 3)
    call sum_check("AggQ_delta" ,   17.958634718_dp  , AggQ_ptr , delta, 3)

    call sum_check("AQqPS_reg" ,    91.950981088_dp  , AQqPS_ptr , reg, 3)
    call sum_check("AqqQPS_reg" ,  -38.624316410_dp  , AqqQPS_ptr , reg, 3)

    call sum_check("AQg_reg" ,    310.17900321_dp    , AQg_ptr , reg, 3)
    call sum_check("AqgQ_reg" ,   -73.710138327_dp   , AqgQ_ptr , reg, 3)
    call sum_check("AgqQ_reg" ,  -120.75198970_dp    , AgqQ_ptr , reg, 3)

    block
      type(grid_conv) :: A3ggQ
      call InitGridConv(grid, A3ggQ, conv_OME(AggQ_ptr,order=3))
    end block
    !! Call the libome functions through the interface
    !res_reg   = ome_AqqQNSEven_reg(as4pi, LM, NF, x)
    !res_plus  = ome_AqqQNSEven_plus(as4pi, LM, NF, x)
    !res_delta = ome_AqqQNSEven_delta(as4pi, LM, NF)

    !! Print results for verification
    !print *, "Results from libome interface:"
    !print *, "AqqQNSEven_reg: ", res_reg
    !print *, "AqqQNSEven_plus:", res_plus
    !print *, "AqqQNSEven_delta:", res_delta

  contains

    subroutine sum_check(name, refval, ptr, piece, max_order)
      use iso_c_binding
      character(len=*), intent(in) :: name
      real(dp)        , intent(in) :: refval
      type(c_ptr), intent(in) :: ptr
      integer(c_int), intent(in) :: piece
      integer(c_int), intent(in) :: max_order
      real(dp) :: res
      integer(c_int) :: order
      real(dp) :: as2pi

      ! integer(c_intptr_t) :: addr
      ! addr = transfer(ptr, addr)
      ! print '(a,I0)', "Pointer address: ", addr
      as2pi = as4pi * 2.0_dp
      res = 0.0_dp
      !print '(a)', "Computing sum for ", trim(name), "ptr = ", addr, " piece=", piece
      do order = 0, max_order
        res = res + (as2pi**order) * ome_piece_hoppet(ptr, log(1.0_dp/x), piece, order, LM, NF)
      end do
      if (piece /= cc_DELTA) res = res / x
      if (piece == cc_VIRT) res = -res     ! because we compare to the plus part, which has the opposite sign
      ! table gives 10 digits after the first digit
      ! actual accuracy (email from Arnd Behring) should be 2048 * epsilon(1.0_dp)
      call check_approx_eq(trim(name), res, refval, tol_abs = 1e-10_dp, tol_rel = 1e-10_dp, tol_choice_or = .true.)

    end subroutine sum_check

  end subroutine test_libome_interface    
end module test_libome