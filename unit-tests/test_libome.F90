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
    real(c_double) :: as4pi, LM, nf_light_cdouble, x, y
    real(c_double) :: res_reg, res_plus, res_delta
    integer(c_int) :: piece
    integer, parameter :: reg = cc_REALVIRT
    integer, parameter :: delta = cc_DELTA
    integer, parameter :: plus = cc_VIRT 

    if (.not. do_test("test_libome_interface")) return

    ! !y=0.46058326383423536798_dp, piece=1, order=3, LM=0, NF=5, NF=3
    ! res_reg = ome_piece_hoppet(AggQ_ptr, y=0.4775571658694307553_dp, piece=1, order=3, LM=0.0_dp, NF=5.0_dp)
    ! write(6,*) "res reg = ", res_reg
    ! return

    ! Example values for testing
    as4pi = 0.1_dp
    LM    = 10._dp
    nf_light_cdouble    = 3.0_dp
    x     = 0.1_dp

    ! check all 3rd order entries from table2 of https://arxiv.org/pdf/2510.02175
    call sum_check("AqqQNSEven_reg" ,    6.4686126472_dp , AqqQNSEven_ptr , reg,   3)
    call sum_check("AqqQNSEven_plus" ,  -4.7428175642_dp , AqqQNSEven_ptr , plus,  3)
    call sum_check("AqqQNSEven_delta" , -2.6422012017_dp , AqqQNSEven_ptr , delta, 3)

    call sum_check("AqqQNSOdd_reg" ,    6.6927772815_dp  , AqqQNSOdd_ptr , reg,   3)
    call sum_check("AqqQNSOdd_plus" ,  -4.7428175642_dp  , AqqQNSOdd_ptr , plus,  3)
    call sum_check("AqqQNSOdd_delta" , -2.6422012017_dp  , AqqQNSOdd_ptr , delta, 3)

    call sum_check("AggQ_reg" ,   -262.65714922_dp   , AggQ_ptr , reg, 3)
    call sum_check("AggQ_plus" ,    -6.6285411655_dp , AggQ_ptr , plus, 3)
    call sum_check("AggQ_delta" ,   17.958634718_dp  , AggQ_ptr , delta, 3)

    call sum_check("AQqPS_reg" ,    91.950981088_dp  , AQqPS_ptr , reg, 3)
    call sum_check("AqqQPS_reg" ,  -38.624316410_dp  , AqqQPS_ptr , reg, 3)

    call sum_check("AQg_reg" ,    310.17900321_dp    , AQg_ptr , reg, 3)
    call sum_check("AqgQ_reg" ,   -73.710138327_dp   , AqgQ_ptr , reg, 3)
    call sum_check("AgqQ_reg" ,  -120.75198970_dp    , AgqQ_ptr , reg, 3)

    call moment_check()
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

    !! sums the different orders from ome_piece_hoppet and compares to reference value
    !!
    !! \param name        = name of the OME (for printing)
    !! \param refval      = reference value to compare to
    !! \param ptr         = pointer to the OME object
    !! \param piece       = which piece (cc_REAL, cc_REALVIRT, cc_PLUS, cc_DELTA)
    !! \param max_order   = maximum order in alphas to sum up to
    !!
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

      as2pi = as4pi * 2.0_dp
      res = 0.0_dp
      do order = 0, max_order
        res = res + (as2pi**order) * ome_piece_hoppet(ptr, log(1.0_dp/x), piece, order, LM, nf_light_cdouble)
      end do
      if (piece /= cc_DELTA) res = res / x
      if (piece == cc_VIRT) res = -res     ! because we compare to the plus part, which has the opposite sign
      ! table gives 10 digits after the first digit
      ! actual accuracy (email from Arnd Behring) should be 2048 * epsilon(1.0_dp)
      call check_approx_eq(trim(name), res, refval, tol_abs = 1e-10_dp, tol_rel = 1e-10_dp, tol_choice_or = .true.)

    end subroutine sum_check


    subroutine moment_check()
      ! to be implemented
      use qcd, only: qcd_SetNf
      use convolution
      use streamlined_interface
      use mass_thresholds_n3lo
      type(grid_conv) :: A2ggQ, A3ggQ, A3piece
      real(dp), pointer :: xfx(:), res(:)
      real(dp) :: moment_N, moment_res, moment_res2
      type(mass_threshold_mat) :: mtm2, mtm3

      real(dp) :: as2pi
      call qcd_SetNf(int(nf_light_cdouble)+1)
      call SetNfDglapHolder(dh,int(nf_light_cdouble)+1)
      as2pi = as4pi * 2.0_dp
      call InitGridConv(grid, A3ggQ) ! initialises to zero

      call InitGridConv(grid, A3piece, conv_OME(AggQ_ptr,order=0, nf_light = nf_light_cdouble, LM=LM)) 
      call AddWithCoeff(A3ggQ, A3piece, as2pi**0)

      call InitGridConv(grid, A3piece, conv_OME(AggQ_ptr,order=1, nf_light = nf_light_cdouble, LM=LM)) 
      call AddWithCoeff(A3ggQ, A3piece, as2pi**1)

      call InitGridConv(grid, A3piece, conv_OME(AggQ_ptr,order=2, nf_light = nf_light_cdouble, LM=LM)) 
      call AddWithCoeff(A3ggQ, A3piece, as2pi**2)

      call InitGridConv(grid, A3piece, conv_OME(AggQ_ptr,order=3, nf_light = nf_light_cdouble, LM=LM)) 
      call AddWithCoeff(A3ggQ, A3piece, as2pi**3)

      moment_N = 2.0_dp
      moment_res = gc_moment(A3ggQ, moment_N)
      write(6,*) red//"Moment N=", moment_N, " of 3rd order AggQ OME = ", moment_res,reset

      call AllocGridQuant(grid,xfx)
      call AllocGridQuant(grid,res)

      
      xfx = exp(moment_N*yValues(grid))
      res = A3ggQ * xfx
      moment_res = res(grid%ny) / xfx(grid%ny)
      write(6,*) "Moment N=", moment_N, " of 3rd order AggQ OME = ", moment_res

      call InitGridConv(grid, A2ggQ, conv_OME(AggQ_ptr,order=2, LM=zero))
      moment_res = gc_moment(A2ggQ, moment_N)
      moment_res2 = gc_moment(dh%allMTM(3,4)%Sgg_H, moment_N)
      write(6,*) "Moment N=", moment_N, " of 2nd order"
      write(6,*) "  AggQ OME        = ", moment_res
      write(6,*) "  dh%allMTM%Sgg_H = ", moment_res2

      if (dh%nloop >= 4) then
        call InitGridConv(grid, A3ggQ, conv_OME(AggQ_ptr,order=3, nf_light = nf_light_cdouble, LM=zero))
        moment_res = gc_moment(A3ggQ, moment_N)
        moment_res2 = gc_moment(dh%allMTM(4,4)%Sgg_H, moment_N)
        write(6,*) "Moment N=", moment_N, " of 3rd order"
        write(6,*) "  AggQ OME        = ", moment_res
        write(6,*) "  dh%allMTM%Sgg_H = ", moment_res2
      end if


      call InitMTMLibOME(grid, mtm2, nloop=3)
      call InitMTMLibOME(grid, mtm3, nloop=4)
      write(6,*) "PShq     ", gc_moment(dh%allMTM(3,nf_int)%PShq, moment_N), gc_moment(mtm2%PShq, moment_N)
      write(6,*) "PShg     ", gc_moment(dh%allMTM(3,nf_int)%PShg, moment_N), gc_moment(mtm2%PShg, moment_N)
      write(6,*) "PSqq_H   ", gc_moment(dh%allMTM(3,nf_int)%PSqq_H, moment_N), gc_moment(mtm2%PSqq_H, moment_N)
      write(6,*) "NSqq_H   ", gc_moment(dh%allMTM(3,nf_int)%NSqq_H, moment_N), gc_moment(mtm2%NSqq_H, moment_N)
      write(6,*) "NSmqq_H  ", gc_moment(dh%allMTM(3,nf_int)%NSmqq_H, moment_N), gc_moment(mtm2%NSmqq_H, moment_N)      
      write(6,*) "Sgg_H    ", gc_moment(dh%allMTM(3,nf_int)%Sgg_H, moment_N), gc_moment(mtm2%Sgg_H, moment_N)
      write(6,*) "Sgq_H    ", gc_moment(dh%allMTM(3,nf_int)%Sgq_H, moment_N), gc_moment(mtm2%Sgq_H, moment_N)
      write(6,*) "Sqg_H    ", gc_moment(dh%allMTM(3,nf_int)%Sqg_H, moment_N), gc_moment(mtm2%Sqg_H, moment_N)

      if (dh%nloop >= 4) then
        write(6,*) "PShq     ", gc_moment(dh%allMTM(4,nf_int)%PShq, moment_N), gc_moment(mtm3%PShq, moment_N)
        write(6,*) "PShg     ", gc_moment(dh%allMTM(4,nf_int)%PShg, moment_N), gc_moment(mtm3%PShg, moment_N)
        write(6,*) "PSqq_H   ", gc_moment(dh%allMTM(4,nf_int)%PSqq_H, moment_N), gc_moment(mtm3%PSqq_H, moment_N)
        write(6,*) "NSqq_H   ", gc_moment(dh%allMTM(4,nf_int)%NSqq_H, moment_N), gc_moment(mtm3%NSqq_H, moment_N)
        write(6,*) "NSmqq_H  ", gc_moment(dh%allMTM(4,nf_int)%NSmqq_H, moment_N), gc_moment(mtm3%NSmqq_H, moment_N)      
        write(6,*) "Sgg_H    ", gc_moment(dh%allMTM(4,nf_int)%Sgg_H, moment_N), gc_moment(mtm3%Sgg_H, moment_N)
        write(6,*) "Sgq_H    ", gc_moment(dh%allMTM(4,nf_int)%Sgq_H, moment_N), gc_moment(mtm3%Sgq_H, moment_N)
        write(6,*) "Sqg_H    ", gc_moment(dh%allMTM(4,nf_int)%Sqg_H, moment_N), gc_moment(mtm3%Sqg_H, moment_N)
      end if

      call check_moment("nloop=3, PShq     ", moment_N,   mtm2%PShq   , dh%allMTM(3,nf_int)%PShq   )
      call check_moment("nloop=3, PShg     ", moment_N,   mtm2%PShg   , dh%allMTM(3,nf_int)%PShg   )
      call check_moment("nloop=3, PSqq_H   ", moment_N,   mtm2%PSqq_H , dh%allMTM(3,nf_int)%PSqq_H )
      call check_moment("nloop=3, NSqq_H   ", moment_N,   mtm2%NSqq_H , dh%allMTM(3,nf_int)%NSqq_H )
      call check_moment("nloop=3, NSmqq_H  ", moment_N,   mtm2%NSmqq_H, dh%allMTM(3,nf_int)%NSmqq_H)      
      call check_moment("nloop=3, Sgg_H    ", moment_N,   mtm2%Sgg_H  , dh%allMTM(3,nf_int)%Sgg_H  )
      call check_moment("nloop=3, Sgq_H    ", moment_N,   mtm2%Sgq_H  , dh%allMTM(3,nf_int)%Sgq_H  )
      call check_moment("nloop=3, Sqg_H    ", moment_N,   mtm2%Sqg_H  , dh%allMTM(3,nf_int)%Sqg_H  )
      if (dh%nloop >= 4) then
        call check_moment("nloop=4, PShq     ", moment_N,   mtm3%PShq   , dh%allMTM(4,nf_int)%PShq   )
        call check_moment("nloop=4, PShg     ", moment_N,   mtm3%PShg   , dh%allMTM(4,nf_int)%PShg   )
        call check_moment("nloop=4, PSqq_H   ", moment_N,   mtm3%PSqq_H , dh%allMTM(4,nf_int)%PSqq_H )
        call check_moment("nloop=4, NSqq_H   ", moment_N,   mtm3%NSqq_H , dh%allMTM(4,nf_int)%NSqq_H )
        call check_moment("nloop=4, NSmqq_H  ", moment_N,   mtm3%NSmqq_H, dh%allMTM(4,nf_int)%NSmqq_H)      
        call check_moment("nloop=4, Sgg_H    ", moment_N,   mtm3%Sgg_H  , dh%allMTM(4,nf_int)%Sgg_H  )
        call check_moment("nloop=4, Sgq_H    ", moment_N,   mtm3%Sgq_H  , dh%allMTM(4,nf_int)%Sgq_H  )
        call check_moment("nloop=4, Sqg_H    ", moment_N,   mtm3%Sqg_H  , dh%allMTM(4,nf_int)%Sqg_H  )
      end if

    end subroutine moment_check

    subroutine check_moment(name, momN, gc_test, gc_ref)
      use mass_thresholds_n3lo, only : gc_moment
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: momN
      type(grid_conv), intent(in) :: gc_test, gc_ref
      real(dp) :: res_test, res_ref

      res_test = gc_moment(gc_test, momN)
      res_ref  = gc_moment(gc_ref, momN)
      call check_approx_eq(name, answer=res_test, expected = res_ref,&
                                 tol_abs = 1e-7_dp, tol_rel = 1e-7_dp, tol_choice_or = .true.)
    end subroutine check_moment

  end subroutine test_libome_interface    
end module test_libome