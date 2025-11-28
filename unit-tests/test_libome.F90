#include "inc/ftlMacros.inc"
#include "inc/hoppet_config.inc"

module test_libome
  use unit_tests
  use qcd

contains
  !! subroutine for testing the hoppet-libome interface
  !! 
  !! reference values can be found in https://arxiv.org/abs/2510.02175, table 2
  subroutine test_libome_interface
    use hoppet_libome_fortran
    use hoppet_term ! for term colours
    use convolution_pieces
    use convolution
    use streamlined_interface
    use iso_c_binding
    implicit none
    real(c_double) :: as4pi, LM, nf_light_cdouble, x, y
    real(c_double) :: res_reg, res_plus, res_delta
    integer(c_int) :: piece
    integer, parameter :: reg   = cc_REALVIRT
    integer, parameter :: delta = cc_DELTA
    integer, parameter :: plus  = cc_VIRT
    integer(c_int), parameter  :: pieces(3) = [reg, plus, delta]
    integer :: ipiece

    if (.not. do_test("test_libome_interface")) return

    ! Example values for testing
    as4pi = 0.1_dp
    LM    = 10._dp
    nf_light_cdouble    = 3.0_dp
    x     = 0.1_dp

    ! check all 3rd order entries from table 2 of https://arxiv.org/pdf/2510.02175
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

    !! Run checks that the LM=0 and LM=eps results agree, since we take
    !! different paths within the libome interface in these two cases
    do ipiece = 1, size(pieces)
      piece = pieces(ipiece)
      call LM_zero_check("AqqQNSEven"//", piece="//to_string(piece), AqqQNSEven_ptr, piece, 3)
      call LM_zero_check("AqqQNSOdd "//", piece="//to_string(piece), AqqQNSOdd_ptr , piece, 3)
      call LM_zero_check("AggQ      "//", piece="//to_string(piece), AggQ_ptr      , piece, 3)
      call LM_zero_check("AQqPS     "//", piece="//to_string(piece), AQqPS_ptr    , piece, 3)
      call LM_zero_check("AqqQPS    "//", piece="//to_string(piece), AqqQPS_ptr   , piece, 3)
      call LM_zero_check("AQg       "//", piece="//to_string(piece), AQg_ptr      , piece, 3)
      call LM_zero_check("AqgQ      "//", piece="//to_string(piece), AqgQ_ptr     , piece, 3)
      call LM_zero_check("AgqQ      "//", piece="//to_string(piece), AgqQ_ptr     , piece, 3)
    end do

    call libome_moment_check_table4()

    call moment_check()

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


    subroutine LM_zero_check(name, ptr, piece, order)
      use iso_c_binding
      character(len=*), intent(in) :: name
      type(c_ptr),      intent(in) :: ptr
      integer(c_int),   intent(in) :: piece
      integer(c_int),   intent(in) :: order
      !--
      real(dp), parameter :: LM_eps = 1e-50_dp
      real(dp) :: res_LM0, res_LMeps
      ! to be implemented

      res_LM0   = ome_piece_hoppet(ptr, log(1.0_dp/x), piece, order, LM=0.0_dp , NF=nf_light_cdouble)
      res_LMeps = ome_piece_hoppet(ptr, log(1.0_dp/x), piece, order, LM=LM_eps, NF=nf_light_cdouble)
      call check_approx_eq(name//": LM=0 vs LM=eps", res_LM0, res_LMeps, tol_abs = 1e-10_dp, tol_rel = 1e-10_dp, tol_choice_or = .true.)
    end subroutine LM_zero_check


    subroutine moment_check()
      ! to be implemented
      use qcd, only: qcd_SetNf
      use convolution
      use streamlined_interface
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM      
      use mass_thresholds_n3lo
#endif
      real(dp) :: moment_N
      integer  :: imoment_N
      type(mass_threshold_mat) :: mtm2, mtm3, mtm3_fortran

      call InitMTMLibOME(grid, mtm2, nloop=3)
      call InitMTMLibOME(grid, mtm3, nloop=4)
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM      
      call InitMTMN3LOExactFortran(grid, mtm3_fortran) ! fortran exact except AQg
#endif


      do imoment_N = 2, 3
        moment_N = real(imoment_N, dp)

        call check_moment("nloop=3, PShq     ", moment_N,   mtm2%PShq   , dh%allMTM(3,nf_int)%PShq   )
        call check_moment("nloop=3, PShg     ", moment_N,   mtm2%PShg   , dh%allMTM(3,nf_int)%PShg   )
        call check_moment("nloop=3, PSqq_H   ", moment_N,   mtm2%PSqq_H , dh%allMTM(3,nf_int)%PSqq_H )
        call check_moment("nloop=3, NSqq_H   ", moment_N,   mtm2%NSqq_H , dh%allMTM(3,nf_int)%NSqq_H )
        call check_moment("nloop=3, NSmqq_H  ", moment_N,   mtm2%NSmqq_H, dh%allMTM(3,nf_int)%NSmqq_H)      
        call check_moment("nloop=3, Sgg_H    ", moment_N,   mtm2%Sgg_H  , dh%allMTM(3,nf_int)%Sgg_H  )
        call check_moment("nloop=3, Sgq_H    ", moment_N,   mtm2%Sgq_H  , dh%allMTM(3,nf_int)%Sgq_H  )
        call check_moment("nloop=3, Sqg_H    ", moment_N,   mtm2%Sqg_H  , dh%allMTM(3,nf_int)%Sqg_H  )
  
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM      
        ! AQg can differ at the 6e-6 relative level; in the fortran OME code, it is the only
        ! only one that is not exact. It appears to use expansions and, depending on the
        ! system and optimisation level, one can trigger differences, perhaps because of
        ! a combination of rounding errors and how they induce adaptive integration fluctuations
        call check_moment("nloop=4, PShg     ", moment_N,   mtm3%PShg   , mtm3_fortran%PShg, override_tol=1e-5_dp ) 
        ! all others should be precise.
        call check_moment("nloop=4, PShq     ", moment_N,   mtm3%PShq   , mtm3_fortran%PShq   )
        call check_moment("nloop=4, PSqq_H   ", moment_N,   mtm3%PSqq_H , mtm3_fortran%PSqq_H )
        call check_moment("nloop=4, NSqq_H   ", moment_N,   mtm3%NSqq_H , mtm3_fortran%NSqq_H )
        call check_moment("nloop=4, NSmqq_H  ", moment_N,   mtm3%NSmqq_H, mtm3_fortran%NSmqq_H)      
        call check_moment("nloop=4, Sgg_H    ", moment_N,   mtm3%Sgg_H  , mtm3_fortran%Sgg_H  )
        call check_moment("nloop=4, Sgq_H    ", moment_N,   mtm3%Sgq_H  , mtm3_fortran%Sgq_H  )
        call check_moment("nloop=4, Sqg_H    ", moment_N,   mtm3%Sqg_H  , mtm3_fortran%Sqg_H  )
#endif
      end do
    end subroutine moment_check

    subroutine check_moment(name, momN, gc_test, gc_ref, override_tol)
      use assertions, only : default_or_opt
      character(len=*),   intent(in) :: name
      real(dp),           intent(in) :: momN
      type(grid_conv),    intent(in) :: gc_test, gc_ref
      real(dp), optional, intent(in) :: override_tol
      !--
      real(dp) :: res_test, res_ref
      real(dp) :: tol

      tol = default_or_opt(1e-7_dp, override_tol)

      ! the hoppet convention is that the Nth moment is
      ! \int_0^1 dx x^(N-1) xpdf(x)
      ! which differs by one from the 
      res_test = gc_moment(gc_test, momN-one)
      res_ref  = gc_moment(gc_ref, momN-one)
      call check_approx_eq(name, answer=res_test, expected = res_ref, &
                                 tol_abs = tol, tol_rel = tol, tol_choice_or = .true.)
    end subroutine check_moment

  end subroutine test_libome_interface    


  !------------------------------------------------------------------------
  !! run checks against table 4 of https://arxiv.org/abs/2510.02175
  subroutine libome_moment_check_table4()
    use hoppet_libome_fortran
    use qcd, only: qcd_SetNf
    use convolution
    use iso_c_binding
    !--
    type(grid_def) :: grid_bigy
    real(dp) :: as2pi
    real(c_double) :: as4pi, LM, nf_light_cdouble


    as4pi = 0.10_dp
    LM    = 10._dp
    nf_light_cdouble    = 3.0_dp

    call qcd_SetNf(int(nf_light_cdouble)+1)
    as2pi = as4pi * 2.0_dp

    ! set up a grid with a big y range to get good accuracy for low moments
    ! and with a small dy to get good accuracy for high moments
    call InitGridDef(grid_bigy, dy=0.02_dp, ymax=24.0_dp, order=-6)

    call check_moment_sum("AqqQNSEven", AqqQNSEven_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[2.552050805_dp, 5.205091550_dp, 6.935755383_dp])

    call check_moment_sum("AQqPS     ", AQqPS_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[5.787728157_dp, 1.517439651_dp, 0.726319249_dp])

    call check_moment_sum("AqqQPS    ", AqqQPS_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[-1.927340435_dp, -0.287734919_dp, -0.111383030_dp])

    call check_moment_sum("AQg       ", AQg_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[-4.267201565_dp, -19.75022897_dp, -21.85267609_dp])

    call check_moment_sum("AqgQ      ", AqgQ_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[2.466090032_dp, 3.364276844_dp, 2.718842618_dp])

    call check_moment_sum("AgqQ      ", AgqQ_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[-5.412438528_dp, -0.735206982_dp, -0.219958419_dp])

    call check_moment_sum("AggQ      ", AggQ_ptr, &
                          moment_Ns=[2,4,6], &
                          refvals=[2.801111705_dp, 20.85194977_dp, 26.04572674_dp])

    call check_moment_sum("AqqQNSOdd ", AqqQNSOdd_ptr, &
                          moment_Ns=[3,5,7], &
                          refvals=[4.041765636_dp, 6.148438117_dp, 7.610661858_dp])

  contains

    subroutine check_moment_sum(name, ptr, moment_Ns, refvals)
      !use mass_thresholds_n3lo

      character(len=*), intent(in) :: name
      type(c_ptr),      intent(in) :: ptr
      integer,          intent(in) :: moment_Ns(:)
      real(dp),         intent(in) :: refvals(:)
      !----
      integer, parameter :: max_order = 3
      integer :: iorder, i
      integer(c_int) :: order
      type(grid_conv) :: A3sum, A3piece
      real(dp) :: moment_N, res


      call InitGridConv(grid_bigy, A3sum) !! sets it to zero
      do iorder = 0, max_order
        call InitGridConv(grid_bigy, A3piece, conv_OME(ptr, order=iorder, nf_light=nf_light_cdouble, LM=LM))        
        call AddWithCoeff(A3sum, A3piece, as2pi**iorder)
      end do

      do i = 1, size(moment_Ns)
        moment_N = real(moment_Ns(i), dp)
        res = gc_moment(A3sum, moment_N-one)
        call check_approx_eq(&
                    "libome "//name//", momN="//to_string(moment_Ns(i)), &
                     answer=res, expected=refvals(i), &
                     tol_abs=1e-6_dp, tol_rel=1e-6_dp, tol_choice_or=.true.)
      end do

      call Delete(A3sum)
      call Delete(A3piece)
    end subroutine check_moment_sum

  end subroutine libome_moment_check_table4

end module test_libome