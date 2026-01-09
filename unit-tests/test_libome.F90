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

    ! for additions in 2512.13508, there's no table, so we compare against 
    ! a reference value in libome/tests/AQqPSs.test.cpp
    as4pi = 0.25_dp
    LM    = -5.0_dp
    nf_light_cdouble = 3.0_dp
    x = 0.0625_dp
    call sum_check("AQqPSs_reg" ,  -1.42143670137779381366114343548_dp  , AQqPSs_ptr , reg, 3)

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
      call LM_zero_check("AQqPSs    "//", piece="//to_string(piece), AQqPSs_ptr   , piece, 3)
    end do

    call libome_moment_check_table4()

    call moment_check()

    call libome_QQbar_valence_check()

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
      type(grid_conv) :: mtm2_psqq, mtm3_psqq, dh_psqq, mtm3_fortran_psqq

      call InitMTMLibOME(grid, mtm2, nloop=3)
      call InitMTMLibOME(grid, mtm3, nloop=4)
      call InitGridConv(dh_psqq, dh%allMTM(3,nf_int)%P_light%qq)
      call AddWithCoeff(dh_psqq, dh%allMTM(3,nf_int)%P_light%NS_plus,-one)
      call InitGridConv(mtm2_psqq, mtm2%P_light%qq)
      call AddWithCoeff(mtm2_psqq, mtm2%P_light%NS_plus,-one)
      call InitGridConv(mtm3_psqq, mtm3%P_light%qq)
      call AddWithCoeff(mtm3_psqq, mtm3%P_light%NS_plus,-one)

#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM      
      call InitMTMN3LOExactFortran(grid, mtm3_fortran) ! fortran exact except AQg
      call InitGridConv(mtm3_fortran_psqq, mtm3_fortran%P_light%qq)
      call AddWithCoeff(mtm3_fortran_psqq, mtm3_fortran%P_light%NS_plus,-one)
#endif


      do imoment_N = 2, 3
        moment_N = real(imoment_N, dp)

        call check_moment("nloop=3, PShq     ", moment_N,   mtm2%PShq   , dh%allMTM(3,nf_int)%PShq   )
        call check_moment("nloop=3, PShg     ", moment_N,   mtm2%PShg   , dh%allMTM(3,nf_int)%PShg   )
        call check_moment("nloop=3, PSqq_H   ", moment_N,   mtm2_psqq , dh_psqq )
        call check_moment("nloop=3, NSqq_H   ", moment_N,   mtm2%P_light%NS_plus , dh%allMTM(3,nf_int)%P_light%NS_plus )
        call check_moment("nloop=3, NSmqq_H  ", moment_N,   mtm2%P_light%NS_minus, dh%allMTM(3,nf_int)%P_light%NS_plus)
        call check_moment("nloop=3, Sgg_H    ", moment_N,   mtm2%P_light%gg  , dh%allMTM(3,nf_int)%P_light%gg  )
        call check_moment("nloop=3, Sgq_H    ", moment_N,   mtm2%P_light%gq  , dh%allMTM(3,nf_int)%P_light%gq  )
        call check_moment("nloop=3, Sqg_H    ", moment_N,   mtm2%P_light%qg  , dh%allMTM(3,nf_int)%P_light%qg  )
  
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM      
        ! AQg can differ at the 6e-6 relative level; in the fortran OME code, it is the only
        ! only one that is not exact. It appears to use expansions and, depending on the
        ! system and optimisation level, one can trigger differences, perhaps because of
        ! a combination of rounding errors and how they induce adaptive integration fluctuations
        call check_moment("nloop=4, PShg     ", moment_N,   mtm3%PShg   , mtm3_fortran%PShg, override_tol=1e-5_dp ) 
        ! all others should be precise.
        call check_moment("nloop=4, PShq     ", moment_N,   mtm3%PShq   , mtm3_fortran%PShq   )
        call check_moment("nloop=4, PShq     ", moment_N,   mtm3%NShV   , mtm3_fortran%NShV   )
        call check_moment("nloop=4, PSqq_H   ", moment_N,   mtm3_psqq   , mtm3_fortran_psqq )
        call check_moment("nloop=4, NSqq_H   ", moment_N,   mtm3%P_light%NS_plus , mtm3_fortran%P_light%NS_plus )
        call check_moment("nloop=4, NSmqq_H  ", moment_N,   mtm3%P_light%NS_minus, mtm3_fortran%P_light%NS_minus)      
        call check_moment("nloop=4, Sgg_H    ", moment_N,   mtm3%P_light%gg  , mtm3_fortran%P_light%gg  )
        call check_moment("nloop=4, Sgq_H    ", moment_N,   mtm3%P_light%gq  , mtm3_fortran%P_light%gq  )
        call check_moment("nloop=4, Sqg_H    ", moment_N,   mtm3%P_light%qg  , mtm3_fortran%P_light%qg  )
#endif
      end do
      call Delete(dh_psqq)
      call Delete(mtm2_psqq)
      call Delete(mtm3_psqq)
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM
      call Delete(mtm3_fortran_psqq)
#endif
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


  !! This routine carries out checks of the $\alpha_s^3 L$ contribution
  !! to valence heavy-quark (and light-quark) evolution, checking
  !! consistency between the splitting-function based determination
  !! and the MTM-based determination in libOME.
  !!
  !! The essence of the test is look at
  !!
  !! - the difference between nf=4 and nf=3 evolution of a toy PDF with
  !!   only u-ubar valence content as the initial condition
  !! - the difference between the libOME MTM with LM=non-zero and LM=zero,
  !!   applied to the same toy PDF
  !!
  !! The two results should agree for the c-cbar valence component,
  !! and the s-sbar valence component should remain zero (by construction
  !! in the MTM case).
  !!
  subroutine libome_QQbar_valence_check()
    ! to be implemented
    use qcd, only: qcd_SetNf
    use pdf_representation
    use convolution
    use streamlined_interface

    real(dp) :: toy_pdf(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: dtoy_nf3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: dtoy_nf4(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: dtoy_mtm(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: xvals(0:grid%ny), x, xtest(3) = [0.01_dp, 0.1_dp, 0.5_dp]
    real(dp) :: ln_Q2_over_M2 = 4.3_dp
    type(mass_threshold_mat) :: mtm_lm0, mtm_lmQ2
    integer  :: ix

    xvals = xValues(grid)

    ! set up a toy PDF with u - ubar valence only
    toy_pdf = zero
    toy_pdf(:,iflv_u) = xvals**0.5_dp * (1.0_dp - xvals)**3
    toy_pdf(:,iflv_ubar) = -toy_pdf(:,iflv_u )

    ! indices are nloop,nf
    dtoy_nf3 = ln_Q2_over_M2 * (dh%allP(3,3) * toy_pdf)
    dtoy_nf4 = ln_Q2_over_M2 * (dh%allP(3,4) * toy_pdf)
    ! then get _just_ the part 
    dtoy_nf4 = dtoy_nf4 - dtoy_nf3

    call qcd_SetNf(4) ! set nf=4 means the number of flavours including the heavy one
    call InitMTMLibOME(grid, mtm_lm0 , nloop=4, LM= 0.0_dp)
    call InitMTMLibOME(grid, mtm_lmQ2, nloop=4, LM=-ln_Q2_over_M2)
    ! subtract off the two results
    call AddWithCoeff(mtm_lmQ2, mtm_lm0, -1.0_dp)
    dtoy_mtm = mtm_lmQ2 * toy_pdf

    x = 0.1_dp
    do ix = 1, size(xtest)
      x = xtest(ix)
      call check_approx_eq("libome v evolution LM!=0 Q-Qbar valence check at x="//to_string(x), &
                           answer=(dtoy_mtm(:, iflv_c) - dtoy_mtm(:, iflv_cbar)).atx.(x.with.grid),&
                           expected=(dtoy_nf4(:, iflv_c) - dtoy_nf4(:, iflv_cbar)).atx.(x.with.grid),&
                           tol_abs=1e-6_dp, tol_rel=1e-6_dp, tol_choice_or=.true.)

      ! this one isn't really a test of libOME, but rather that the
      ! absence of a specific term in libOME is consistent with a check
      ! that, starting from just u-ubar valencem the s-sbar valence
      ! evolves independently of nf at order $alpha_s^3 L$
      call check_approx_eq("NNLO evolution (nf4-nf3 from pure u-ubar) s-sbar valence=0 check at x="//to_string(x), &
                           answer=(dtoy_nf4(:, iflv_s) - dtoy_nf4(:, iflv_sbar)).atx.(x.with.grid),&
                           expected=zero,&
                           tol_abs=1e-8_dp, tol_choice_or=.true.)

      !write(6,*) "x =", x
      !write(6,*) "c-cbar (from split):", (dtoy_nf4(:, iflv_c) - dtoy_nf4(:, iflv_cbar)).atx.(x.with.grid)
      !write(6,*) "c-cbar (from MTM  ):", (dtoy_mtm(:, iflv_c) - dtoy_mtm(:, iflv_cbar)).atx.(x.with.grid)
      !write(6,*) "s-sbar (from split):", (dtoy_nf4(:, iflv_s) - dtoy_nf4(:, iflv_sbar)).atx.(x.with.grid)
      !write(6,*) "s-sbar (from MTM  ):", (dtoy_mtm(:, iflv_s) - dtoy_mtm(:, iflv_sbar)).atx.(x.with.grid)
    end do

  end subroutine libome_QQbar_valence_check
end module test_libome