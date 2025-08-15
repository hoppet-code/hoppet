!======================================================================
!!
!! Module for providing access to the deep inelastic scattering
!! structure functions up to N3LO
!!

module structure_functions
  use pdf_representation
  use streamlined_interface
  use splitting_functions
  use coefficient_functions_holder
  use qcd
  implicit none

  private
  !!
  !! public functions
  !! those for initialisation and I/O
  public :: InitStrFct, StartStrFct, write_f1, write_f2, write_f3
  !!
  !! those for evaluating the structure functions
  !!
  public :: F_LO, F_NLO, F_NNLO, F_N3LO, StrFct
  public :: F_LO_flav, F_NLO_flav, StrFct_flav
  !!
  !! those for querying the scale choices being used
  !!
  public :: sf_muR, sf_muF
  
  !!
  !! indices for the different structure functions (D=d,s,b, U=u,c)
  !!
  public :: iF1Wp, iF2Wp, iF3Wp, iF1Wm, iF2Wm, iF3Wm, iF1Z, iF2Z,&
       & iF3Z, iF1EM, iF2EM, iF1gZ, iF2gZ, iF3gZ
  integer, parameter :: iF1Wp = 1 !< F1 W+ : D + Ubar 
  integer, parameter :: iF2Wp = 2 !< F2 W+ : D + Ubar 
  integer, parameter :: iF3Wp = 3 !< F3 W+ : D + Ubar 
  integer, parameter :: iF1Wm =-1 !< F1 W- : Dbar + U 
  integer, parameter :: iF2Wm =-2 !< F2 W- : Dbar + U 
  integer, parameter :: iF3Wm =-3 !< F3 W- : Dbar + U 
  integer, parameter :: iF1Z  = 4 !< F1 Z  : (D + Dbar) * v_i^2a_i^2_down    + (U + Ubar) * v_i^2a_i^2_up
  integer, parameter :: iF2Z  = 5 !< F2 Z  : (D + Dbar) * v_i^2a_i^2_down    + (U + Ubar) * v_i^2a_i^2_up
  integer, parameter :: iF3Z  = 6 !< F3 Z  : (D + Dbar) * 2v_i_a_i_down      + (U + Ubar) * 2v_i_a_i_up
  integer, parameter :: iF1EM =-4 !< F1 γ  : (D + Dbar) * e2_down            + (U + Ubar) * e2_up
  integer, parameter :: iF2EM =-5 !< F2 γ  : (D + Dbar) * e2_down            + (U + Ubar) * e2_up
  integer, parameter :: iF1gZ = 0 !< F1 γZ : (D + Dbar) * e_down * 2v_i_down + (U + Ubar) * e_up * 2v_i_up
  integer, parameter :: iF2gZ =-6 !< F2 γZ : (D + Dbar) * e_down * 2v_i_down + (U + Ubar) * e_up * 2v_i_up
  integer, parameter :: iF3gZ = 7 !< F3 γZ : (D + Dbar) * e_down * 2a_i_down + (U + Ubar) * e_up * 2a_i_up

  !!
  !! coupling constants and fixed parameters
  !!
  ! public indices for scale choice options
  integer, parameter, public :: scale_choice_fixed     = 0 !< muR,muF scales predetermined in the StartStrFct call
  integer, parameter, public :: scale_choice_Q         = 1 !< muR,muF scales equal to xR,xQ * Q (xR,xQ to be set in StartStrFct)
  integer, parameter, public :: scale_choice_arbitrary = 2 !< muR,muF scales can be chosen freely in the F_LO (etc.) and StrFct calls

  !!
  !! holds the coefficient functions
  !!
  type(coef_holder), save :: ch

  !!
  !! Holds the structure functions. The tables are identical to the
  !! streamlined tables, but separate such that one can use the
  !! streamlined tables at the same time as the structure functions.
  !!
  !! These contain the SF summed over flavour, and dressed with
  !! quark couplings according to the indices above (iF1Wp etc)
  type(pdf_table), save :: sf_tables(0:max_table_index) 
  !! These contained the SF decomposed by flavour, WITHOUT
  !! couplings. For now functionality for these are somewhat
  !! limited. They are only available at LO/NLO.
  type(pdf_table), save :: sf_tables_flav(0:max_table_index)
  logical,         save :: sf_alloc_already_done = .false.

  !!
  !! constants and fixed parameters
  !!
  logical,  save      :: use_mass_thresholds
  real(dp), parameter :: viW = 1/sqrt(two), aiW = viW
  real(dp), save      :: vi2_ai2_avg_W, two_vi_ai_avg_W
  real(dp), save      :: vi2_ai2_Z_down, vi2_ai2_Z_up
  real(dp), save      :: two_vi_ai_Z_down, two_vi_ai_Z_up
  real(dp), save      :: two_vi_Z_down, two_vi_Z_up
  real(dp), save      :: two_ai_Z_down, two_ai_Z_up
  real(dp), parameter :: e_up = 2.0_dp/3.0_dp, e_down = - 1.0_dp/3.0_dp
  real(dp), parameter :: e2_up = 4.0_dp/9.0_dp, e2_down = 1.0_dp/9.0_dp

  !!
  !! these log terms are only used for scale_choice = 0,1, to speed up the code
  !!
  real(dp), save      :: log_muF2_over_Q2, log_muR2_over_Q2
  !real(dp), save      :: Qmin
  real(dp), public, save :: toy_alphas_Q0 = 0.35_dp ! Roughly αS(sqrt(2))
  type(running_coupling), public, save :: toy_coupling

  !!
  !! scale choices
  !!
  real(dp), save :: cst_mu, xmuR, xmuF
  integer, save  :: scale_choice_save
  
  logical, save  :: use_sep_orders, inc_flavour_decomposition
  logical, save  :: exact_coefs
  integer :: nf_lcl, order_setup

contains
  !> @brief Setup of constants and parameters needed for structure functions
  !!
  !! @param[in]      order_max      highest order in QCD to compute (1: LO, 2: NLO, etc)
  !! @param[opt]     nflav          integer number of flavours (if not present use variable flavour)
  !! @param[opt]     scale_choice   (0: fixed scale, 1: use Q, 2: use arbitrary scale)
  !! @param[opt]     constant_mu    if scale_choice = scale_choice_fixed (= 0) then this is the fixed scale
  !! @param[opt]     param_coefs    if .true. use parametrised coefficients functions
  !! @param[opt]     wmass          Mass of the W boson
  !! @param[opt]     zmass          Mass of the z boson
  !!
  subroutine StartStrFct(order_max, nflav, scale_choice, constant_mu,&
       & param_coefs, wmass, zmass)
    integer, intent(in)  :: order_max
    integer, optional    :: nflav, scale_choice
    real(dp), optional   :: constant_mu, wmass, zmass
    logical, optional    :: param_coefs
    !----------------------------------------------------------------------
    real(dp) :: sin_thw_sq, mw, mz

    if(.not.alloc_already_done) then
       call wae_error('StartStrFct', 'hoppetStart not called!')
    endif

    ! Read input
    ! mW and mZ defaults are taken from the 2022 PDG 
    mw                = default_or_opt(80.377_dp, wmass)
    mz                = default_or_opt(91.1876_dp, zmass)
    scale_choice_save = default_or_opt(scale_choice_Q, scale_choice)
    exact_coefs       = .not.default_or_opt(.true., param_coefs)
    cst_mu            = default_or_opt(zero, constant_mu)
        
    ! compute sin(\theta_w) from W/Z mass
    sin_thw_sq = 1.0_dp - (mw/mz)**2
    
    ! evaluate parameters needed for the structure functions
    ! cf. Eqs. (3.10+11) of 1109.3717
    vi2_ai2_Z_up     = one/four + (half - (four/three) * sin_thw_sq)**2
    vi2_ai2_Z_down   = one/four + (half - (two/three)  * sin_thw_sq)**2
    two_vi_ai_Z_up   = half - (four/three) * sin_thw_sq
    two_vi_ai_Z_down = half - (two/three)  * sin_thw_sq

    ! cf. Eq.3.20 + 3.17 & 3.18
    ! 1/nf \sum_j=1^nf (vj^2 + aj^2)
    vi2_ai2_avg_W = viW**2 + aiW**2
    ! 1/nf \sum_j=1^nf 2*vj*aj
    two_vi_ai_avg_W = two * viW * aiW

    ! these parameters are needed explicitly for the gamma-Z structure function
    two_vi_Z_up = one - (8.0_dp/three) * sin_thw_sq
    two_vi_Z_down   = -one + (four/three) * sin_thw_sq
    two_ai_Z_up = one
    two_ai_Z_down   = -one

    if (default_or_opt(-1, nflav) > 0) then
       ! if nflav is present, then use a fixed flavour number
       call qcd_SetNf(nflav)
       call InitCoefHolder(grid, ch, order_max, exact_coefs)
       use_mass_thresholds = .false.
       nf_lcl = nf_int
       write(6,'(a)'   ) "Starting the structure functions with fixed number of flavours"
       write(6,'(a,i2)') "nf = ", nf_lcl
    else
       ! otherwise, use variable flavour number scheme
       !    call InitCoefHolder(grid, ch, order_max, exact_coefs, nflo=3, nfhi=5)
       call InitCoefHolder(grid, ch, order_max, exact_coefs, nflo=3, nfhi=6)
       use_mass_thresholds = .true.
       
       write(6,'(a)'      ) "Starting the structure functions with mass thresholds at"
       write(6,'(a,f12.5)') "mc = ", masses(4) 
       write(6,'(a,f12.5)') "mb = ", masses(5)
       write(6,'(a,f12.5)') "mt = ", masses(6) 
       ! and start with a sensible local nf (which will be 5 here) 
       nf_lcl = nf_int
    endif
    
    ! if mass thresholds are used in the structure functions, only scale_choice 0 and 1 are allowed
    if ((use_mass_thresholds).and.(scale_choice_save.ge.scale_choice_arbitrary)) then
       call wae_error('StartStrFct', 'illegal value for scale_choice with mass thresholds turned on', intval = scale_choice_save)
    end if
 
  end subroutine StartStrFct

  !> @brief Initialize the structure functions up to specified order
  !!
  !! Initialize the structure functions up to specified order
  !! this requires the tabulated PDF to have been set up beforehand in the
  !! the streamlined interface. 
  !!
  !! By default, separate_orders = .false., which means that one
  !! can only access the sum over all orders. If separate_orders = .true.
  !! then one can access each order separately, but this is slower.
  !!
  !! @param[in]      order              order at which we initialise the structure functions
  !! @param[opt]     separate_orders    if .true. separate into individual orders rather than summing
  !! @param[opt]     xR             factor to multiply renormalisation scale
  !! @param[opt]     xF             factor to multiply factorisation scale
  !!
  subroutine InitStrFct(order, separate_orders, xR, xF, flavour_decomposition)
    integer, intent(in) :: order
    real(dp), optional   :: xR, xF
    logical, optional :: separate_orders, flavour_decomposition

    if(.not.setup_done(0)) then
       call wae_error('InitStrFct', 'hoppetEvolve not called!')
    endif

    xmuR              = default_or_opt(one, xR)
    xmuF              = default_or_opt(one, xF)

    inc_flavour_decomposition = .false.
    inc_flavour_decomposition = default_or_opt(inc_flavour_decomposition, flavour_decomposition)

    if(sf_alloc_already_done) then
       call Delete(sf_tables)
       if (inc_flavour_decomposition) call Delete(sf_tables_flav)
    endif

    ! Now we set up the tables. HoppetStart already called, so we can copy over the structure
    call AllocPdfTable(sf_tables, tables(0))
    if(inc_flavour_decomposition) call AllocPdfTable(sf_tables_flav, tables(0))
    ! Finally we need to set tab_iflv_max = 7. We don't change tables(0) as it contains the PDF
    !sf_tables(0) = tables(0)
    !sf_tables_flav(0) = tables(0)
    sf_tables(1:)%tab_iflv_max = 7
    if(inc_flavour_decomposition) sf_tables_flav(1:)%tab_iflv_max = 6
    sf_alloc_already_done = .true. ! Signals that the tables have been set up.

    ! To turn off b quarks completely (only for testing and comparison)
    ! uncomment following two lines:
    ! tables(0)%tab(:,-5,:) = zero
    ! tables(0)%tab(:,+5,:) = zero

    ! default is to sum each order before convolutions (it's faster)
    if (scale_choice_save == scale_choice_arbitrary) then
       use_sep_orders = .true.
    else
       use_sep_orders = .false.
    endif
    use_sep_orders = default_or_opt(use_sep_orders, separate_orders)

    if (use_sep_orders) then
       ! First we treat the case where we want to separate out each order in
       ! the final structure functions.
       ! This is slower, but is needed e.g. for VBFH
       if (scale_choice_save == scale_choice_Q .or. scale_choice_save == scale_choice_fixed) then
          ! if scale_choice = 0,1 use fast implementation
          ! tables is saved as an array in Q, and only tables(0),
          ! sf_tables(1), sf_tables(2), sf_tables(4), sf_tables(8) are non zero
          if (order.ge.1) call set_LO_structure_functions()
          if (order.ge.2) call set_NLO_structure_functions()
          if (order.ge.3) call set_NNLO_structure_functions()
          if (order.ge.4) call set_N3LO_structure_functions()
       else if(scale_choice_save == scale_choice_arbitrary) then
          ! if scale_choice >= 2 use slower implementation with full
          ! scale choices, such as sqrt(Q1*Q2)
          ! tables is saved as an array in muF now, and all components
          ! of sf_tables are non zero.
          if (order.ge.1) call set_LO_structure_functions_anyscale()
          if (order.ge.2) call set_NLO_structure_functions_anyscale()
          if (order.ge.3) call set_NNLO_structure_functions_anyscale()
          if (order.ge.4) call set_N3LO_structure_functions_anyscale()
       else
          call wae_error('InitStrFct', 'illegal value for scale_choice', intval = scale_choice_save)
       endif
    else
       ! Now set up the default case where we sum up everything right away.
       if (scale_choice_save >= scale_choice_arbitrary) & ! only allow for scale_choice = 0,1
            & call wae_error('InitStrFct', &
            & 'with scale_choice_arbitrary (or higher), separate_orders must be true, but was false; scale_choice=', &
            & intval = scale_choice_save)
       call set_structure_functions_upto(order)
    endif
    
    ! Now for flavour decomposed SFs
    if(inc_flavour_decomposition) then
       if (use_sep_orders) then
          ! First we treat the case where we want to separate out each order in
          ! the final structure functions.
          ! This is slower, but is needed e.g. for VBFH
          if (scale_choice_save == scale_choice_Q .or. scale_choice_save == scale_choice_fixed) then
             ! if scale_choice = 0,1 use fast implementation
             ! tables is saved as an array in Q, and only tables(0),
             ! sf_tables(1), sf_tables(2), sf_tables(4), sf_tables(8) are non zero
             if (order.ge.1) call set_LO_structure_functions_flav()
             if (order.ge.2) call set_NLO_structure_functions_flav()
          else if(scale_choice_save == scale_choice_arbitrary) then
             ! if scale_choice >= 2 use slower implementation with full
             ! scale choices, such as sqrt(Q1*Q2)
             ! tables is saved as an array in muF now, and all components
             ! of sf_tables are non zero.
             if (order.ge.1) call set_LO_structure_functions_anyscale_flav()
             if (order.ge.2) call set_NLO_structure_functions_anyscale_flav()
          else
             call wae_error('InitStrFct', 'illegal value for scale_choice', intval = scale_choice_save)
          endif
       else
          ! Now set up the default case where we sum up everything right away.
          if (scale_choice_save >= scale_choice_arbitrary) & ! only allow for scale_choice = 0,1
               & call wae_error('InitStrFct', &
               & 'with scale_choice_arbitrary (or higher), separate_orders must be true, but was false; scale_choice=', &
               & intval = scale_choice_save)
          call set_structure_functions_upto_flav(order)
       endif
    endif

    ! To rescale PDF by the N3LO F2 structure function (evaluated at 8 GeV)
    ! as a check of the size of missing N3LO PDFs, uncomment following line:
    ! call rescale_pdf_nnlo(8.0_dp,order)
    ! call rescale_pdf_n3lo(8.0_dp,order)

    setup_done(1:) = .true.
    order_setup = order
    
  end subroutine InitStrFct
  
  !> @brief write the F1 structure function to idev
  !!
  !! Writes the W and Z F1 structure functions to a file with equal spacing in log(x)
  !!
  !! @param[in]     idev     file/device to write to
  !! @param[in]     Qtest    Value of Q to use
  !! @param[in]     ymax     Largest value of y = - log(x) to print
  !! @param[in]     ny       number of bins in y to print
  !!
  subroutine write_f1 (idev, Qtest, ymax, ny, muF, muR)
    real(dp), intent(in) :: Qtest, ymax
    real(dp), intent(in), optional :: muF, muR
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO, res(-6:7)
    real(dp) :: F1EM_LO, F1EM_NLO, F1EM_NNLO, F1EM_N3LO, F1gZ_LO, F1gZ_NLO, F1gZ_NNLO, F1gZ_N3LO 
    integer  :: iy
    !F1 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    if (use_sep_orders) then
       write(idev,'(a,a)') '#     x        F1Wp(LO)    F1Wm(LO)    F1Wp(NLO)   F1Wm(NLO)  F1Wp(NNLO)  F1Wm(NNLO)', &
            & '  F1Wp(N3LO)  F1Wm(N3LO)    F1Z(LO)    F1Z(NLO)    F1Z(NNLO)&
            &   F1Z(N3LO)    F1γ(LO)    F1γ(NLO)    F1γ(NNLO)   F1γ(N3LO)   F1γZ(LO)&
            &   F1γZ(NLO)   F1γZ(NNLO)  F1γZ(N3LO)'
    else
       write(idev,'(a)') '# x F1Wp F1Wm F1Z F1γ F1γZ'
    endif
    mF = default_or_opt(sf_muF(Qtest),muF)
    mR = default_or_opt(sf_muR(Qtest),muR)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       if (use_sep_orders) then
          res = F_LO(xval, Qtest, mR, mF)
          write(idev,'(3es12.4)',advance='no') xval, res(iF1Wp),res(iF1Wm)
          F1Z_LO = res(iF1Z)
          F1EM_LO = res(iF1EM)
          F1gZ_LO = res(iF1gZ)
          res = F_NLO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF1Wp), res(iF1Wm)
          F1Z_NLO = res(iF1Z)
          F1EM_NLO = res(iF1EM)
          F1gZ_NLO = res(iF1gZ)
          res = F_NNLO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF1Wp), res(iF1Wm)
          F1Z_NNLO = res(iF1Z)
          F1EM_NNLO = res(iF1EM)
          F1gZ_NNLO = res(iF1gZ)
          res = F_N3LO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF1Wp), res(iF1Wm)
          F1Z_N3LO = res(iF1Z)
          F1EM_N3LO = res(iF1EM)
          F1gZ_N3LO = res(iF1gZ)
          write(idev,'(12es12.4)',advance='no') F1Z_LO, F1Z_NLO,&
               & F1Z_NNLO, F1Z_N3LO, F1EM_LO, F1EM_NLO, F1EM_NNLO,&
               & F1EM_N3LO, F1gZ_LO, F1gZ_NLO, F1gZ_NNLO, F1gZ_N3LO 
       else
          res = StrFct(xval, Qtest, mR, mF)
          write(idev,'(6es12.4)',advance='no') xval, res(iF1Wp),res(iF1Wm), res(iF1Z), res(iF1EM), res(iF1gZ)
       endif
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f1
  
  !> @brief write the F2 structure function to idev
  !!
  !! Writes the W and Z F2 structure functions to a file with equal spacing in log(x)
  !!
  !! @param[in]     idev     file/device to write to
  !! @param[in]     Qtest    Value of Q to use
  !! @param[in]     ymax     Largest value of y = - log(x) to print
  !! @param[in]     ny       number of bins in y to print
  !!
  subroutine write_f2 (idev, Qtest, ymax, ny, muF, muR)
    real(dp), intent(in) :: Qtest, ymax
    real(dp), intent(in), optional :: muF, muR
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO, res(-6:7)
    real(dp) :: F2EM_LO, F2EM_NLO, F2EM_NNLO, F2EM_N3LO, F2gZ_LO, F2gZ_NLO, F2gZ_NNLO, F2gZ_N3LO 
    integer  :: iy
    !F2 Wp Wm Z
    write(idev,'(a,f20.4,a,f20.4)') '# Q = ', Qtest
    if (use_sep_orders) then
       write(idev,'(a,a)') '#     x        F2Wp(LO)    F2Wm(LO)    F2Wp(NLO)   F2Wm(NLO)  F2Wp(NNLO)  F2Wm(NNLO)', &
            & '  F2Wp(N3LO)  F2Wm(N3LO)    F2Z(LO)    F2Z(NLO)    F2Z(NNLO)&
            &   F2Z(N3LO)    F2γ(LO)    F2γ(NLO)    F2γ(NNLO)   F2γ(N3LO)   F2γZ(LO)&
            &   F2γZ(NLO)   F2γZ(NNLO)  F2γZ(N3LO)'
    else
       write(idev,'(a)') '# x F2Wp F2Wm F2Z F2γ F2γZ'
    endif
    mF = default_or_opt(sf_muF(Qtest),muF)
    mR = default_or_opt(sf_muR(Qtest),muR)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       if (use_sep_orders) then
          res = F_LO(xval, Qtest, mR, mF)
          write(idev,'(3es12.4)',advance='no') xval, res(iF2Wp),res(iF2Wm)
          F2Z_LO = res(iF2Z)
          F2EM_LO = res(iF2EM)
          F2gZ_LO = res(iF2gZ)
          res = F_NLO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF2Wp), res(iF2Wm)
          F2Z_NLO = res(iF2Z)
          F2EM_NLO = res(iF2EM)
          F2gZ_NLO = res(iF2gZ)
          res = F_NNLO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF2Wp), res(iF2Wm)
          F2Z_NNLO = res(iF2Z)
          F2EM_NNLO = res(iF2EM)
          F2gZ_NNLO = res(iF2gZ)
          res = F_N3LO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF2Wp), res(iF2Wm)
          F2Z_N3LO = res(iF2Z)
          F2EM_N3LO = res(iF2EM)
          F2gZ_N3LO = res(iF2gZ)
          write(idev,'(12es12.4)',advance='no') F2Z_LO, F2Z_NLO,&
               & F2Z_NNLO, F2Z_N3LO, F2EM_LO, F2EM_NLO, F2EM_NNLO,&
               & F2EM_N3LO, F2gZ_LO, F2gZ_NLO, F2gZ_NNLO, F2gZ_N3LO 
       else
          res = StrFct(xval, Qtest, mR, mF)
          write(idev,'(6es12.4)',advance='no') xval, res(iF2Wp),res(iF2Wm), res(iF2Z), res(iF2EM), res(iF2gZ)
       endif
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f2
!  subroutine write_f2 (idev, Qtest, ymax, ny, muF, muR)
!    real(dp), intent(in) :: Qtest, ymax
!    real(dp), intent(in), optional :: muF, muR
!    integer, intent(in)  :: idev, ny
!    real(dp) :: ytest, xval, mR, mF, F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO, res(-6:7)
!    integer  :: iy
!    !F2 Wp Wm Z
!    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
!    if (use_sep_orders) then
!       write(idev,'(a,a)') '# x  F2Wp(LO) F2Wm(LO) F2Wp(NLO) F2Wm(NLO) F2Wp(NNLO) F2Wm(NNLO)', &
!            & ' F2Wp(N3LO) F2Wm(N3LO) F2Z(LO) F2Z(NLO) F2Z(NNLO) F2Z(N3LO)'
!    else
!       write(idev,'(a)') '# x F2Wp F2Wm F2Z'
!    endif
!    mF = default_or_opt(sf_muF(Qtest),muF)
!    mR = default_or_opt(sf_muR(Qtest),muR)
!    do iy = ny, 1, -1
!       ytest = iy * ymax / ny
!       xval = exp(-ytest)
!       if (use_sep_orders) then
!          res = F_LO(xval, Qtest, mR, mF)
!          write(idev,'(3es12.4)',advance='no') xval, res(iF2Wp),res(iF2Wm)
!          F2Z_LO = res(iF2Z)
!          res = F_NLO(xval, Qtest, mR, mF)
!          write(idev,'(2es12.4)',advance='no') res(iF2Wp), res(iF2Wm)
!          F2Z_NLO = res(iF2Z)
!          res = F_NNLO(xval, Qtest, mR, mF)
!          write(idev,'(2es12.4)',advance='no') res(iF2Wp), res(iF2Wm)
!          F2Z_NNLO = res(iF2Z)
!          res = F_N3LO(xval, Qtest, mR, mF)
!          write(idev,'(2es12.4)',advance='no') res(iF2Wp), res(iF2Wm)
!          F2Z_N3LO = res(iF2Z)
!          write(idev,'(4es12.4)',advance='no') F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO
!       else
!          res = StrFct(xval, Qtest, mR, mF)
!          write(idev,'(4es12.4)',advance='no') xval, res(iF2Wp),res(iF2Wm), res(iF2Z)
!       endif
!       write(idev,*)
!    end do
!    write(idev,*)
!    write(idev,*)
!  end subroutine write_f2

  !> @brief write the F3 structure function to idev
  !!
  !! Writes the W and Z F3 structure functions to a file with equal spacing in log(x)
  !!
  !! @param[in]     idev     file/device to write to
  !! @param[in]     Qtest    Value of Q to use
  !! @param[in]     ymax     Largest value of y = - log(x) to print
  !! @param[in]     ny       number of bins in y to print
  !!
    subroutine write_f3 (idev, Qtest, ymax, ny, muF, muR)
    real(dp), intent(in) :: Qtest, ymax
    real(dp), intent(in), optional :: muF, muR
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO, res(-6:7)
    real(dp) :: F3gZ_LO, F3gZ_NLO, F3gZ_NNLO, F3gZ_N3LO 
    integer  :: iy
    !F3 Wp Wm Z
    write(idev,'(a,f30.4,a,f30.4)') '# Q = ', Qtest
    if (use_sep_orders) then
       write(idev,'(a,a)') '#     x        F3Wp(LO)    F3Wm(LO)    F3Wp(NLO)   F3Wm(NLO)  F3Wp(NNLO)  F3Wm(NNLO)', &
            & '  F3Wp(N3LO)  F3Wm(N3LO)    F3Z(LO)    F3Z(NLO)    F3Z(NNLO)&
            &   F3Z(N3LO)    F3γZ(LO)&
            &   F3γZ(NLO)  F3γZ(NNLO)  F3γZ(N3LO)'
    else
       write(idev,'(a)') '# x F3Wp F3Wm F3Z F3γZ'
    endif
    mF = default_or_opt(sf_muF(Qtest),muF)
    mR = default_or_opt(sf_muR(Qtest),muR)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       if (use_sep_orders) then
          res = F_LO(xval, Qtest, mR, mF)
          write(idev,'(3es12.4)',advance='no') xval, res(iF3Wp),res(iF3Wm)
          F3Z_LO = res(iF3Z)
          F3gZ_LO = res(iF3gZ)
          res = F_NLO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF3Wp), res(iF3Wm)
          F3Z_NLO = res(iF3Z)
          F3gZ_NLO = res(iF3gZ)
          res = F_NNLO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF3Wp), res(iF3Wm)
          F3Z_NNLO = res(iF3Z)
          F3gZ_NNLO = res(iF3gZ)
          res = F_N3LO(xval, Qtest, mR, mF)
          write(idev,'(2es12.4)',advance='no') res(iF3Wp), res(iF3Wm)
          F3Z_N3LO = res(iF3Z)
          F3gZ_N3LO = res(iF3gZ)
          write(idev,'(8es12.4)',advance='no') F3Z_LO, F3Z_NLO,&
               & F3Z_NNLO, F3Z_N3LO, F3gZ_LO, F3gZ_NLO, F3gZ_NNLO,&
               & F3gZ_N3LO 
       else
          res = StrFct(xval, Qtest, mR, mF)
          write(idev,'(6es12.4)',advance='no') xval, res(iF3Wp),res(iF3Wm), res(iF3Z), res(iF3gZ)
       endif
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f3

!  subroutine write_f3 (idev, Qtest, ymax, ny, muF, muR)
!    real(dp), intent(in) :: Qtest, ymax
!    real(dp), intent(in), optional :: muF, muR
!    integer, intent(in)  :: idev, ny
!    real(dp) :: ytest, xval, mR, mF, F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO, res(-6:7)
!    integer  :: iy
!    !F3 Wp Wm Z
!    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
!    if (use_sep_orders) then
!       write(idev,'(a,a)') '# x  F3Wp(LO) F3Wm(LO) F3Wp(NLO) F3Wm(NLO) F3Wp(NNLO) F3Wm(NNLO)', &
!            & ' F3Wp(N3LO) F3Wm(N3LO) F3Z(LO) F3Z(NLO) F3Z(NNLO) F3Z(N3LO)'
!    else
!       write(idev,'(a)') '# x F3Wp F3Wm F3Z'
!    endif
!    mF = default_or_opt(sf_muF(Qtest),muF)
!    mR = default_or_opt(sf_muR(Qtest),muR)
!    do iy = ny, 1, -1
!       ytest = iy * ymax / ny
!       xval = exp(-ytest)
!       if (use_sep_orders) then
!          res = F_LO(xval, Qtest, mR, mF)
!          write(idev,'(3es12.4)',advance='no') xval, res(iF3Wp),res(iF3Wm)
!          F3Z_LO = res(iF3Z)
!          res = F_NLO(xval, Qtest, mR, mF)
!          write(idev,'(2es12.4)',advance='no') res(iF3Wp), res(iF3Wm)
!          F3Z_NLO = res(iF3Z)
!          res = F_NNLO(xval, Qtest, mR, mF)
!          write(idev,'(2es12.4)',advance='no') res(iF3Wp), res(iF3Wm)
!          F3Z_NNLO = res(iF3Z)
!          res = F_N3LO(xval, Qtest, mR, mF)
!          write(idev,'(2es12.4)',advance='no') res(iF3Wp), res(iF3Wm)
!          F3Z_N3LO = res(iF3Z)
!          write(idev,'(4es12.4)',advance='no') F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO
!       else
!          res = StrFct(xval, Qtest, mR, mF)
!          write(idev,'(4es12.4)',advance='no') xval, res(iF3Wp),res(iF3Wm), res(iF3Z)
!       endif
!       write(idev,*)
!    end do
!    write(idev,*)
!    write(idev,*)
!  end subroutine write_f3

  !----------------------------------------------------------------------
  ! set up structure functions for scale_choice = 0, 1, adding up each order
  subroutine set_structure_functions_upto (order)
    integer, intent(in) :: order
    integer :: iQ
    real(dp) :: Q, as2pi
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    
    do iQ = 0, tables(0)%nQ
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       call set_scale_logs(Q)
       
       as2pi = alphasLocal(sf_muR(Q)) / (twopi)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
   
       ! start with LO
       if (order.ge.1) then
          sf_tables(1)%tab(:,:,iQ) = structure_function_general(ch%C2LO*f, ch%CLLO*f, ch%C3LO*f)
       endif
       
       ! now add NLO terms
       if (order.ge.2) then
          if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
             ! For central scale with scale_choice=1, we don't need to do all
             ! the convolutions of the splitting functions
             f2 = ch%C2NLO * f
             fL = ch%CLNLO * f
             f3 = ch%C3NLO * f
          else
             ! do the convolution with the coefficient functions and also the
             ! corresponding splitting-function contributions when scales
             ! are not equal to Q
             PLO_f = dh%P_LO * f
             f2 = CxNLO_with_logs(ch%C2LO, ch%C2NLO, f, PLO_f)
             fL = CxNLO_with_logs(ch%CLLO, ch%CLNLO, f, PLO_f)
             f3 = CxNLO_with_logs(ch%C3LO, ch%C3NLO, f, PLO_f)
          endif

          sf_tables(1)%tab(:,:,iQ) = sf_tables(1)%tab(:,:,iQ) + &
               & as2pi * structure_function_general(f2, fL, f3)
       endif

       ! now add NNLO terms
       if (order.ge.3) then
          if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
             ! For central scale with scale_choice=1, we don't need to do all
             ! the convolutions of the splitting functions
             f2 = ch%C2NNLO * f
             fL = ch%CLNNLO * f
             f3 = ch%C3NNLO * f
          else
             ! do the convolution with the coefficient functions and also the
             ! corresponding splitting-function contributions when scales
             ! are not equal to Q
             PLO_f   = dh%P_LO  * f
             PNLO_f  = dh%P_NLO * f
             PLO2_f  = dh%P_LO  * PLO_f
             f2 = CxNNLO_with_logs(ch%C2LO, ch%C2NLO, ch%C2NNLO, f, PLO_f, PNLO_f, PLO2_f)
             fL = CxNNLO_with_logs(ch%CLLO, ch%CLNLO, ch%CLNNLO, f, PLO_f, PNLO_f, PLO2_f)
             f3 = CxNNLO_with_logs(ch%C3LO, ch%C3NLO, ch%C3NNLO, f, PLO_f, PNLO_f, PLO2_f)
          endif

          sf_tables(1)%tab(:,:,iQ) = sf_tables(1)%tab(:,:,iQ) + &
               (as2pi**2) * structure_function_general(f2, fL, f3)
       endif

       ! and finally, add the N3LO terms
       if (order.ge.4) then
          if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
             ! For central scale with scale_choice=1 (i.e. mur=muf=Q), we don't
             ! need to do all the convolutions of the splitting functions
             f2 = ch%C2N3LO * f
             fL = ch%CLN3LO * f
             f3 = ch%C3N3LO * f
          else
             ! do the convolution with the coefficient functions and also the
             ! corresponding splitting-function contributions when scales
             ! are not equal to Q
             PLO_f    = dh%P_LO   * f
             PNLO_f   = dh%P_NLO  * f
             PLO2_f   = dh%P_LO   * PLO_f
             PNNLO_f  = dh%P_NNLO * f
             PLONLO_f = dh%P_LO   * PNLO_f
             PNLOLO_f = dh%P_NLO  * PLO_f
             PLO3_f   = dh%P_LO   * PLO2_f

             f2 = CxN3LO_with_logs(ch%C2LO, ch%C2NLO, ch%C2NNLO, ch%C2N3LO, f, PLO_f, PNLO_f, PLO2_f, &
                  &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
             fL = CxN3LO_with_logs(ch%CLLO, ch%CLNLO, ch%CLNNLO, ch%CLN3LO, f, PLO_f, PNLO_f, PLO2_f, &
                  &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
             f3 = CxN3LO_with_logs(ch%C3LO, ch%C3NLO, ch%C3NNLO, ch%C3N3LO, f, PLO_f, PNLO_f, PLO2_f, &
                  &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
          endif

          ! now the fl_11 piece
          f2_fl11 = ch%C2N3LO_fl11 * f
          fL_fl11 = ch%CLN3LO_fl11 * f

          sf_tables(1)%tab(:,:,iQ) = sf_tables(1)%tab(:,:,iQ) + &
               & (as2pi**3) * structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)
       endif
    end do
  end subroutine set_structure_functions_upto


  !----------------------------------------------------------------------
  ! set up structure functions for scale_choice = 0, 1, adding up each order
  subroutine set_structure_functions_upto_flav (order)
    integer, intent(in) :: order
    integer :: iQ
    real(dp) :: Q, as2pi
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    
    do iQ = 0, tables(0)%nQ
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       call set_scale_logs(Q)
       
       as2pi = alphasLocal(sf_muR(Q)) / (twopi)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
   
       ! start with LO
       if (order.ge.1) then
          sf_tables_flav(1)%tab(:,:,iQ) = ch%CLLO*f
          sf_tables_flav(2)%tab(:,:,iQ) = ch%C2LO*f
          sf_tables_flav(3)%tab(:,:,iQ) = ch%C3LO*f
       endif
       
       ! now add NLO terms
       if (order.ge.2) then
          if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
             ! For central scale with scale_choice=1, we don't need to do all
             ! the convolutions of the splitting functions
             f2 = ch%C2NLO * f
             fL = ch%CLNLO * f
             f3 = ch%C3NLO * f
          else
             ! do the convolution with the coefficient functions and also the
             ! corresponding splitting-function contributions when scales
             ! are not equal to Q
             PLO_f = dh%P_LO * f
             f2 = CxNLO_with_logs(ch%C2LO, ch%C2NLO, f, PLO_f)
             fL = CxNLO_with_logs(ch%CLLO, ch%CLNLO, f, PLO_f)
             f3 = CxNLO_with_logs(ch%C3LO, ch%C3NLO, f, PLO_f)
          endif

          sf_tables_flav(1)%tab(:,:,iQ) = sf_tables_flav(1)%tab(:,:,iQ) + &
               & as2pi * fl
          sf_tables_flav(2)%tab(:,:,iQ) = sf_tables_flav(2)%tab(:,:,iQ) + &
               & as2pi * f2
          sf_tables_flav(3)%tab(:,:,iQ) = sf_tables_flav(3)%tab(:,:,iQ) + &
               & as2pi * f3
       endif

       if(order.ge.3) call wae_error("Flavour decomposed SFs only implemented up to NLO.")
    end do
  end subroutine set_structure_functions_upto_flav

  
  !----------------------------------------------------------------------
  ! set up LO structure functions for scale_choice = 0, 1
  subroutine set_LO_structure_functions()
    integer :: iQ
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       Q = tables(0)%Q_vals(iQ)
       ! explicitly evaluate the PDF at scale sf_muF(Q)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       sf_tables(1)%tab(:,:,iQ) = structure_function_general(ch%C2LO*f, ch%CLLO*f, ch%C3LO*f)
    end do
    
  end subroutine set_LO_structure_functions
  
  !----------------------------------------------------------------------
  ! set up LO structure functions for scale_choice = 0, 1 decomposed in flavour
  subroutine set_LO_structure_functions_flav()
    integer :: iQ
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       Q = tables(0)%Q_vals(iQ)
       ! explicitly evaluate the PDF at scale sf_muF(Q)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       sf_tables_flav(1)%tab(:,:,iQ) = ch%CLLO*f
       sf_tables_flav(2)%tab(:,:,iQ) = ch%C2LO*f
       sf_tables_flav(3)%tab(:,:,iQ) = ch%C3LO*f
    end do
    
  end subroutine set_LO_structure_functions_flav
  
  !----------------------------------------------------------------------
  ! Set up the LO structure functions for any scale choice
  subroutine set_LO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       f = tables(0)%tab(:,:,iQ)
       
       if (use_mass_thresholds) then
          call use_vfns(f, tables(0)%Q_vals(iQ))
       endif
       
       sf_tables(1)%tab(:,:,iQ) = structure_function_general(&
            & f * ch%C2LO, &
            & f * ch%CLLO, &
            & f * ch%C3LO)       
    end do
    
  end subroutine set_LO_structure_functions_anyscale
  !----------------------------------------------------------------------
  ! Set up the LO structure functions for any scale choice flavour decomposed
  subroutine set_LO_structure_functions_anyscale_flav()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       f = tables(0)%tab(:,:,iQ)
       
       if (use_mass_thresholds) then
          call use_vfns(f, tables(0)%Q_vals(iQ))
       endif
       
       sf_tables_flav(1)%tab(:,:,iQ) = ch%CLLO*f
       sf_tables_flav(2)%tab(:,:,iQ) = ch%C2LO*f
       sf_tables_flav(3)%tab(:,:,iQ) = ch%C3LO*f
    end do
    
  end subroutine set_LO_structure_functions_anyscale_flav
  
  !----------------------------------------------------------------------
  ! set up NLO structure functions for scale_choice = 0, 1
  subroutine set_NLO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif

       if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
       ! For central scale with scale_choice=1, we don't need to do all
       ! the convolutions of the splitting functions
          f2 = ch%C2NLO * f
          fL = ch%CLNLO * f
          f3 = ch%C3NLO * f
       else
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
          PLO_f = dh%P_LO * f
          f2 = CxNLO_with_logs(ch%C2LO, ch%C2NLO, f, PLO_f)
          fL = CxNLO_with_logs(ch%CLLO, ch%CLNLO, f, PLO_f)
          f3 = CxNLO_with_logs(ch%C3LO, ch%C3NLO, f, PLO_f)
       endif
          
       sf_tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine set_NLO_structure_functions

  !----------------------------------------------------------------------
  ! set up NLO structure functions for scale_choice = 0, 1 with flavour decomposition
  subroutine set_NLO_structure_functions_flav()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif

       if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
       ! For central scale with scale_choice=1, we don't need to do all
       ! the convolutions of the splitting functions
          f2 = ch%C2NLO * f
          fL = ch%CLNLO * f
          f3 = ch%C3NLO * f
       else
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
          PLO_f = dh%P_LO * f
          f2 = CxNLO_with_logs(ch%C2LO, ch%C2NLO, f, PLO_f)
          fL = CxNLO_with_logs(ch%CLLO, ch%CLNLO, f, PLO_f)
          f3 = CxNLO_with_logs(ch%C3LO, ch%C3NLO, f, PLO_f)
       endif
          
       sf_tables_flav(4)%tab(:,:,iQ) = fL
       sf_tables_flav(5)%tab(:,:,iQ) = f2
       sf_tables_flav(6)%tab(:,:,iQ) = f3

    end do

  end subroutine set_NLO_structure_functions_flav

  !----------------------------------------------------------------------
  ! Set up the NLO structure functions for any scale choice
  subroutine set_NLO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       ! Internal Q value effectively corresponds to sf_muF(Q1,Q2)
       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)

       ! Save the NLO pieces in sf_tables(2) and sf_tables(3)
      
       ! Get the NLO coefficient function, (C_NLO x f) 
       f2 = (ch%C2NLO * f)
       fL = (ch%CLNLO * f)
       f3 = (ch%C3NLO * f)
       sf_tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now compute the (C_LO x P_LO x f) term
       PLO_f = dh%P_LO * f
       f2 = (ch%C2LO * PLO_f)
       fL = (ch%CLLO * PLO_f)
       f3 = (ch%C3LO * PLO_f)
       sf_tables(3)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
    end do

  end subroutine set_NLO_structure_functions_anyscale
  
  !----------------------------------------------------------------------
  ! Set up the NLO structure functions for any scale choice
  subroutine set_NLO_structure_functions_anyscale_flav()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       ! Internal Q value effectively corresponds to sf_muF(Q1,Q2)
       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)

       ! Save the NLO pieces in sf_tables(2) and sf_tables(3)
      
       ! Get the NLO coefficient function, (C_NLO x f) 
       f2 = (ch%C2NLO * f)
       fL = (ch%CLNLO * f)
       f3 = (ch%C3NLO * f)
       sf_tables_flav(4)%tab(:,:,iQ) = fL
       sf_tables_flav(5)%tab(:,:,iQ) = f2
       sf_tables_flav(6)%tab(:,:,iQ) = f3

       ! Now compute the (C_LO x P_LO x f) term
       PLO_f = dh%P_LO * f
       f2 = (ch%C2LO * PLO_f)
       fL = (ch%CLLO * PLO_f)
       f3 = (ch%C3LO * PLO_f)
       sf_tables_flav(7)%tab(:,:,iQ) = fL
       sf_tables_flav(8)%tab(:,:,iQ) = f2
       sf_tables_flav(9)%tab(:,:,iQ) = f3
    end do

  end subroutine set_NLO_structure_functions_anyscale_flav
  
  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for NLO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxNLO_with_logs(CxLO, CxNLO, f, PLO_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------

    res = CxNLO * f - (CxLO * log_muF2_over_Q2) * PLO_f
  end function CxNLO_with_logs
  
  !----------------------------------------------------------------------
  ! set up the NNLO structure functions for scale_choice = 0, 1
  subroutine set_NNLO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ

       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       
       if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
       ! For central scale with scale_choice=1, we don't need to do all
       ! the convolutions of the splitting functions
          f2 = ch%C2NNLO * f
          fL = ch%CLNNLO * f
          f3 = ch%C3NNLO * f
       else
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
          PLO_f   = dh%P_LO  * f
          PNLO_f  = dh%P_NLO * f
          PLO2_f  = dh%P_LO  * PLO_f
          f2 = CxNNLO_with_logs(ch%C2LO, ch%C2NLO, ch%C2NNLO, f, PLO_f, PNLO_f, PLO2_f)
          fL = CxNNLO_with_logs(ch%CLLO, ch%CLNLO, ch%CLNNLO, f, PLO_f, PNLO_f, PLO2_f)
          f3 = CxNNLO_with_logs(ch%C3LO, ch%C3NLO, ch%C3NNLO, f, PLO_f, PNLO_f, PLO2_f)
       endif

       sf_tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
     
    end do

  end subroutine set_NNLO_structure_functions
  
  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for NNLO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxNNLO_with_logs(CxLO, CxNLO, CxNNLO, f, PLO_f, PNLO_f, PLO2_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    type(split_mat), intent(in) :: CxNNLO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------
    real(dp) :: f_NLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f_NNLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: LR, LF

    LR = log_muR2_over_Q2
    LF = log_muF2_over_Q2
    
    f_NLO  = - LF * PLO_f

    f_NNLO = + half * LF**2 * PLO2_f &
         &   - LF * PNLO_f &
         &   - (twopi*beta0*(LR*LF - half * LF**2)) * PLO_f
    
    res = CxNNLO * f  +  CxNLO * f_NLO  +  CxLO * f_NNLO  +  (twopi*beta0*LR) * (CxNLO * f)
  end function CxNNLO_with_logs


  !----------------------------------------------------------------------
  ! set up the NNLO structure functions
  subroutine set_NNLO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    do iQ = 0, tables(0)%nQ

       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)
       
       ! save the NNLO pieces in sf_tables(4:7)
       
       PLO2_f  = dh%P_LO  * (dh%P_LO * f)
       PLO_f   = dh%P_LO  * f
       PNLO_f  = dh%P_NLO * f

       ! first calculate the pure NNLO term, (C_NNLO x f) 
       f2 = (ch%C2NNLO * f)
       fL = (ch%CLNNLO * f)
       f3 = (ch%C3NNLO * f)
       sf_tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO^2 x f) term
       f2 =  (ch%C2LO * PLO2_f)
       fL =  (ch%CLLO * PLO2_f)
       f3 =  (ch%C3LO * PLO2_f)
       sf_tables(5)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO) term
       f2 = (ch%C2LO * PNLO_f)
       fL = (ch%CLLO * PNLO_f)
       f3 = (ch%C3LO * PNLO_f)
       sf_tables(6)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO) term
       f2 = (ch%C2NLO * PLO_f)
       fL = (ch%CLNLO * PLO_f)
       f3 = (ch%C3NLO * PLO_f)
       sf_tables(7)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine set_NNLO_structure_functions_anyscale

  !----------------------------------------------------------------------
  ! set up the N3LO structure functions for scale_choice = 0, 1
  ! Warning : for now factorisation and renormalisation scale are set to Q
  subroutine set_N3LO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ

       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),sf_muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       
       if ((scale_choice_save.eq.scale_choice_Q).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
       ! For central scale with scale_choice=1 (i.e. mur=muf=Q), we don't
       ! need to do all the convolutions of the splitting functions
          f2 = ch%C2N3LO * f
          fL = ch%CLN3LO * f
          f3 = ch%C3N3LO * f
       else
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
          PLO_f    = dh%P_LO   * f
          PNLO_f   = dh%P_NLO  * f
          PLO2_f   = dh%P_LO   * PLO_f
          PNNLO_f  = dh%P_NNLO * f
          PLONLO_f = dh%P_LO   * PNLO_f
          PNLOLO_f = dh%P_NLO  * PLO_f
          PLO3_f   = dh%P_LO   * PLO2_f
          
          f2 = CxN3LO_with_logs(ch%C2LO, ch%C2NLO, ch%C2NNLO, ch%C2N3LO, f, PLO_f, PNLO_f, PLO2_f, &
               &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
          fL = CxN3LO_with_logs(ch%CLLO, ch%CLNLO, ch%CLNNLO, ch%CLN3LO, f, PLO_f, PNLO_f, PLO2_f, &
               &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
          f3 = CxN3LO_with_logs(ch%C3LO, ch%C3NLO, ch%C3NNLO, ch%C3N3LO, f, PLO_f, PNLO_f, PLO2_f, &
               &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       endif
       
       ! now the fl_11 piece
       f2_fl11 = ch%C2N3LO_fl11 * f
       fL_fl11 = ch%CLN3LO_fl11 * f

       !sf_tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       sf_tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)
     
    end do

  end subroutine set_N3LO_structure_functions
  

  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for N3LO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxN3LO_with_logs(CxLO, CxNLO, CxNNLO, CxN3LO, f, PLO_f, PNLO_f, PLO2_f, &
       &                    PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    type(split_mat), intent(in) :: CxNNLO
    type(split_mat), intent(in) :: CxN3LO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------
    real(dp) :: f_NLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f_NNLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f_N3LO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: LR, LF

    LR = log_muR2_over_Q2
    LF = log_muF2_over_Q2
    
    f_NLO  = - LF * PLO_f

    f_NNLO = + half * LF**2 * PLO2_f &
         &   - LF * PNLO_f &
         &   - (twopi*beta0*(LR*LF - half * LF**2)) * PLO_f

    f_N3LO = -(1.0_dp/6.0_dp) * LF * ( &
         & - three * twopi**2 * beta1 * (LF - two * LR) * PLO_f &
         & + two * (twopi*beta0)**2 * (LF**2 - three * LF * LR + three * LR**2) * PLO_f &
         & + LF**2 * PLO3_f &
         & + three * twopi * beta0 * (LF - two * LR) * (LF * PLO2_f - two * PNLO_f) &
         & - 3.0_dp * LF * (PLONLO_f + PNLOLO_f) + 6.0_dp * PNNLO_f)
    
    res = CxN3LO * f  &
         & + (1.0_dp/6.0_dp) * LR * (6.0_dp * twopi**2 * beta1 * (CxNLO * f) &
         & + twopi*beta0 * (12.0_dp * (CxNNLO * f) + 6.0_dp * twopi * beta0 * LR * (CxNLO * f))) &
         & + CxNNLO * f_NLO + twopi * beta0 * LR * (CxNLO * f_NLO) &
         & + CxNLO * f_NNLO + CxLO * f_N3LO
  end function CxN3LO_with_logs

  !----------------------------------------------------------------------
  ! set up the N3LO structure functions
  subroutine set_N3LO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    do iQ = 0, tables(0)%nQ

       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)
       
       ! save the N3LO pieces in sf_tables(8:15)
       
       PLO_f    = dh%P_LO  * f
       PLO2_f   = dh%P_LO  * PLO_f
       PNLO_f   = dh%P_NLO * f
       PNNLO_f  = dh%P_NNLO * f
       PLONLO_f = dh%P_LO   * PNLO_f
       PNLOLO_f = dh%P_NLO  * PLO_f
       PLO3_f   = dh%P_LO   * PLO2_f

       ! first calculate the pure N3LO term, (C_N3LO x f) 
       f2 = (ch%C2N3LO * f)
       fL = (ch%CLN3LO * f)
       f3 = (ch%C3N3LO * f)
       f2_fl11 = (ch%C2N3LO_fl11 * f)
       fL_fl11 = (ch%CLN3LO_fl11 * f)
       ! sf_tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       sf_tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)

       ! Now calculate the (C_LO x P_LO^3 x f) term
       f2 =  (ch%C2LO * PLO3_f)
       fL =  (ch%CLLO * PLO3_f)
       f3 =  (ch%C3LO * PLO3_f)
       sf_tables(9)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO x P_NLO x f) term
       f2 =  (ch%C2LO * PLONLO_f)
       fL =  (ch%CLLO * PLONLO_f)
       f3 =  (ch%C3LO * PLONLO_f)
       sf_tables(10)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO x P_LO x f) term
       f2 =  (ch%C2LO * PNLOLO_f)
       fL =  (ch%CLLO * PNLOLO_f)
       f3 =  (ch%C3LO * PNLOLO_f)
       sf_tables(11)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO^2 x f) term
       f2 =  (ch%C2NLO * PLO2_f)
       fL =  (ch%CLNLO * PLO2_f)
       f3 =  (ch%C3NLO * PLO2_f)
       sf_tables(12)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_NLO) term
       f2 = (ch%C2NLO * PNLO_f)
       fL = (ch%CLNLO * PNLO_f)
       f3 = (ch%C3NLO * PNLO_f)
       sf_tables(13)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NNLO x P_LO) term
       f2 = (ch%C2NNLO * PLO_f)
       fL = (ch%CLNNLO * PLO_f)
       f3 = (ch%C3NNLO * PLO_f)
       sf_tables(14)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NNLO) term
       f2 = (ch%C2LO * PNNLO_f)
       fL = (ch%CLLO * PNNLO_f)
       f3 = (ch%C3LO * PNNLO_f)
       sf_tables(15)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine set_N3LO_structure_functions_anyscale


 !----------------------------------------------------------------------
 ! Structure function valid up to NNLO
  function structure_function_general(C2_f, CL_f, C3_f) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp) :: C2_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: res(0:ubound(C2_f,dim=1), -6:7)
    C2_f_fl11 = zero
    CL_f_fl11 = zero
    res = structure_function_general_full(C2_f, CL_f, C3_f, C2_f_fl11, CL_f_fl11)
  end function structure_function_general


 !----------------------------------------------------------------------
 ! Structure function including N3LO fl_11 piece for neutral current
  function structure_function_general_full(C2_f, CL_f, C3_f, C2_f_fl11, CL_f_fl11) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp), intent(in) :: C2_f_fl11(0:,ncompmin:), CL_f_fl11(0:,ncompmin:)
    real(dp) :: C2_f_NC(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_NC(0:grid%ny,ncompmin:ncompmax)
    real(dp)             :: res(0:ubound(C2_f,dim=1), -6:7)
    !----------------------
    ! not just up and down, but the sum 
    real(dp) :: u(0:grid%ny), d(0:grid%ny), ubar(0:grid%ny), dbar(0:grid%ny) 
    real(dp) :: two_xvals(0:grid%ny)
    integer  :: nf_save
    
    two_xvals = two*xValues(grid) ! move this outside at some point
    C2_f_NC = C2_f + C2_f_fl11
    CL_f_NC = CL_f + CL_f_fl11

    res = zero
    !--- deal with Z case -----------------------------------------
    res(:,iF2Z) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(C2_f_NC) + ubarlike(C2_f_NC))*vi2_ai2_Z_up

    ! temporarily put FL into F1;
    res(:,iF1Z) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(CL_f_NC) + ubarlike(CL_f_NC))*vi2_ai2_Z_up
    ! then convert to F1
    res(:,iF1Z) = (res(:,iF2Z) - res(:,iF1Z)) / two_xvals

    res(:,iF3Z) = (dlike(C3_f) - dbarlike(C3_f))*two_vi_ai_Z_down + &
         &       (ulike(C3_f) - ubarlike(C3_f))*two_vi_ai_Z_up
    res(:,iF3Z) = res(:,iF3Z)/xValues(grid)

    !--- deal with EM case -----------------------------------------
    res(:,iF2EM) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e2_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e2_up
    
    ! temporarily put FL into F1;
    res(:,iF1EM) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e2_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e2_up
    ! then convert to F1
    res(:,iF1EM) = (res(:,iF2EM) - res(:,iF1EM)) / two_xvals
    
    !--- deal with gamma-Z case -----------------------------------------
    ! equations taken from section 19 of PDG (19.18)
    res(:,iF2gZ) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e_down*two_vi_Z_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e_up * two_vi_Z_up
    
    ! temporarily put FL into F1;
    res(:,iF1gZ) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e_down*two_vi_Z_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e_up * two_vi_Z_up
    ! then convert to F1
    res(:,iF1gZ) = (res(:,iF2gZ) - res(:,iF1gZ)) / two_xvals
    
    res(:,iF3gZ) = (dlike(C3_f) - dbarlike(C3_f))*e_down*two_ai_Z_down + &
         &        (ulike(C3_f) - ubarlike(C3_f))*e_up * two_ai_Z_up
    res(:,iF3gZ) = res(:,iF3gZ)/xValues(grid)

    ! for the W cases, it only makes sense to sum over an even number
    ! of light flavours; so save the actual number of flavours, switch
    ! the module-local nf_lcl variable to the nearest (lower) event number
    ! for our W calculations; switch back later
    nf_save = nf_lcl
    nf_lcl = (nf_lcl/2) * 2

    ! AK It looks like we have swapped the meaning of Wp and
    !Wm. Looking at eq. 18.19 in the PDF and eq. 8 in 1206.7007 W-
    !should couple to u + dbar

    ! --- deal with Wm case
    !-----------------------------------------
    res(:,iF2Wm) = (ulike(C2_f) + dbarlike(C2_f))*vi2_ai2_avg_W
    ! temporarily put FL into F1;
    res(:,iF1Wm) = (ulike(CL_f) + dbarlike(CL_f))*vi2_ai2_avg_W
    ! then convert to F1
    res(:,iF1Wm) = (res(:,iF2Wm) - res(:,iF1Wm)) / two_xvals
    res(:,iF3Wm) = (ulike(C3_f) - dbarlike(C3_f))*two_vi_ai_avg_W
    res(:,iF3Wm) = res(:,iF3Wm)/xValues(grid)

    !--- deal with Wp case -----------------------------------------
    res(:,iF2Wp) = (dlike(C2_f) + ubarlike(C2_f))*vi2_ai2_avg_W
    ! temporarily put FL into F1;
    res(:,iF1Wp) = (dlike(CL_f) + ubarlike(CL_f))*vi2_ai2_avg_W
    ! then convert to F1
    res(:,iF1Wp) = (res(:,iF2Wp) - res(:,iF1Wp)) / two_xvals
    res(:,iF3Wp) = (dlike(C3_f) - ubarlike(C3_f))*two_vi_ai_avg_W
    res(:,iF3Wp) = res(:,iF3Wp)/xValues(grid)

!    ! --- deal with Wp case
!    !-----------------------------------------
!    res(:,iF2Wp) = (ulike(C2_f) + dbarlike(C2_f))*vi2_ai2_avg_W
!    ! temporarily put FL into F1;
!    res(:,iF1Wp) = (ulike(CL_f) + dbarlike(CL_f))*vi2_ai2_avg_W
!    ! then convert to F1
!    res(:,iF1Wp) = (res(:,iF2Wp) - res(:,iF1Wp)) / two_xvals
!    res(:,iF3Wp) = (ulike(C3_f) - dbarlike(C3_f))*two_vi_ai_avg_W
!    res(:,iF3Wp) = res(:,iF3Wp)/xValues(grid)
!
!    !--- deal with Wm case -----------------------------------------
!    res(:,iF2Wm) = (dlike(C2_f) + ubarlike(C2_f))*vi2_ai2_avg_W
!    ! temporarily put FL into F1;
!    res(:,iF1Wm) = (dlike(CL_f) + ubarlike(CL_f))*vi2_ai2_avg_W
!    ! then convert to F1
!    res(:,iF1Wm) = (res(:,iF2Wm) - res(:,iF1Wm)) / two_xvals
!    res(:,iF3Wm) = (dlike(C3_f) - ubarlike(C3_f))*two_vi_ai_avg_W
!    res(:,iF3Wm) = res(:,iF3Wm)/xValues(grid)

    ! reset nf_lcl to the full (possibly odd) saved value
    nf_lcl = nf_save

    ! AK 16-08-2023: The factor of two is included in SF in the VBF
    ! paper by Maltoni et al. It was kept here until today since the
    ! structure functions were only used with VBF. Now it has been
    ! removed. WARNING: If running with a version of proVBFH pre this
    ! fix then a factor of 2 will be missing...

    ! overall factor of two that we still haven't fully looked into as
    ! of 2015-02-24 [GPS TMP] 
    ! res = res * two
    !! GPS+AK TMP: we have just included a factor of 2 but this plainly
    !! should not be present for the electromagnetic part; so here
    !! we eliminate it again...
    !res(:,iF2EM) = half * res(:,iF2EM)
    !res(:,iF1EM) = half * res(:,iF1EM)

  end function structure_function_general_full
  

  !----------------------------------------------------------------------
  function ulike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:, 2: nf_lcl: 2), dim=2)
  end function ulike
  !----------------------------------------------------------------------
  function dlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:, 1: nf_lcl: 2), dim=2)
  end function dlike
  !----------------------------------------------------------------------
  function ubarlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:,-2:-nf_lcl:-2), dim=2)
  end function ubarlike
  !----------------------------------------------------------------------
  function dbarlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:,-1:-nf_lcl:-2), dim=2)
  end function dbarlike
  
  
  !----------------------------------------------------------------------
  real(dp) function alphasLocal(muR)
    real(dp), intent(in) :: muR
    !real(dp) :: muR_lcl
    !muR_lcl = max(muR,Qmin)
    ! we use alphas from the LHAPDF PDF
    ! alphasLocal = alphasPDF(muR_lcl)
    ! we use alphas from HOPPET
    alphasLocal = Value(coupling, muR)
  end function alphasLocal
  

  subroutine GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)
   real(dp), intent(in)  :: Q
   real(dp), intent(in),  optional :: muR, muF
   real(dp), intent(out) :: muR_lcl, muF_lcl, mu_table

   ! run a bunch of checks to make sure that the requested scales are consistent
   ! with what is needed here
   if (scale_choice_save == scale_choice_Q) then
      muR_lcl = xmuR * Q
      muF_lcl = xmuF * Q
      mu_table = Q
      if (present(muR)) then
         if (abs(muR/muR_lcl-1).gt.1e-10) call wae_error("GetStrFctScales (scale_choice_Q): ",&
               "requested muR inconsistent with initialisation, muR/(xmuR*Q) = ", &
               dbleval = muR/muR_lcl)
      endif 
      if (present(muF)) then
         if (abs(muF/muF_lcl-1).gt.1e-10) call wae_error("GetStrFctScales (scale_choice_Q): ",&
               "requested muR inconsistent with initialisation, muF/(xmuF*Q) = ", &
               dbleval = muF/muF_lcl)
      endif 
    else if (scale_choice_save == scale_choice_fixed) then
      muR_lcl = xmuR * cst_mu
      muF_lcl = xmuF * cst_mu
      mu_table = Q
      if (present(muR)) then
         if (abs(muR/muR_lcl-1).gt.1e-10) call wae_error("GetStrFctScales (scale_choice_fixed): ",&
               "requested muR inconsistent with initialisation, muR/(xmuR*cst_mu) = ", &
               dbleval = muR/muR_lcl)
      endif 
      if (present(muF)) then
         if (abs(muF/muF_lcl-1).gt.1e-10) call wae_error("GetStrFctScales (scale_choice_fixed): ",&
               "requested muR inconsistent with initialisation, muF/(xmuF*cst_mu) = ", &
               dbleval = muF/muF_lcl)
      endif 
    else if (scale_choice_save == scale_choice_arbitrary) then
      if (.not. present(muR)) call wae_error("GetStrFctScales: muR not present but required for scale_choice_arbitrary")
      if (.not. present(muF)) call wae_error("GetStrFctScales: muF not present but required for scale_choice_arbitrary")
      muR_lcl = muR
      muF_lcl = muF
      mu_table = muF
    else
      call wae_error("GetStrFctScales: scale_choice_save not recognised", intval = scale_choice_save)
    endif

    ! consistency check for use_sep_orders and scale choice
    if (.not. use_sep_orders .and. mu_table .ne. Q) then 
      call wae_error("GetStrFctScales: mu_table /= Q but use_sep_orders = .false.")
    end if
  end subroutine GetStrFctScales

  !---------------------------------------------------------------------
  !> @brief calculate the structure function at x, Q, muR, muF summed over all orders
  !!
  !! Calculate the structure function at x, Q, muR, muF summed over
  !! all orders. muR and muF are only needed if we are using the
  !! scale_choice_arbitrary, as otherwise they are already included in
  !! the sf_tables.
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @return          an array of all structure functions summed over orders
  !!
  function StrFct (x, Q, muR, muF) result(res)
    real(dp), intent(in)  :: x, Q
    real(dp), intent(in), optional  :: muR, muF
    real(dp) :: res(-6:7)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if (use_sep_orders) then
      ! if we kept each order separate, then add up all the fixed order terms one by one
      if (order_setup.ge.1) res =       F_LO  (x, Q, muR_lcl, muF_lcl)
      if (order_setup.ge.2) res = res + F_NLO (x, Q, muR_lcl, muF_lcl)
      if (order_setup.ge.3) res = res + F_NNLO(x, Q, muR_lcl, muF_lcl)
      if (order_setup.ge.4) res = res + F_N3LO(x, Q, muR_lcl, muF_lcl)
    else
      ! if we haven't kept each order separate, then everything is in sf_tables(1)
      call EvalPdfTable_xQ(sf_tables(1), x, mu_table, res)
    endif
    
  end function StrFct

  !---------------------------------------------------------------------
  !> @brief calculate the structure function at x, Q, muR, muF summed over all orders
  !!
  !! Calculate the structure function at x, Q, muR, muF summed over
  !! all orders. muR and muF are only needed if we are using the
  !! scale_choice_arbitrary, as otherwise they are already included in
  !! the sf_tables.
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in,opt]   muR        renormalisation scale 
  !! @param[in,opt]   muF        factorisation scale
  !! @param[in]       iflav      parton-flavour
  !! @return          an array of FL, F2, F3 for the given flavour without couplings to the parton
  !!
  function StrFct_flav (x, Q, muR, muF, iflav) result(res)
    real(dp), intent(in) :: x, Q
    real(dp), intent(in), optional :: muR, muF
    integer, intent(in) :: iflav
    real(dp) :: res(1:3)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    if(.not.inc_flavour_decomposition) call wae_error('StrFct_flav', &
               'You did not initialise the Structure Functions with flavour decomposition. Exiting.')

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if (use_sep_orders) then
      ! if we kept each order separate, then add up all the fixed order terms one by one
      if (order_setup.ge.1) res =       F_LO_flav  (x, Q, muR_lcl, muF_lcl, iflav)
      if (order_setup.ge.2) res = res + F_NLO_flav (x, Q, muR_lcl, muF_lcl, iflav)
      if (order_setup.ge.3) call wae_error('StrFct_flav: Flavour decomposed structure functions only implemneted up to NLO.')
    else
      ! if we haven't kept each order separate, then everything is in sf_tables(1:3)
      res(1) = EvalPdfTable_xQf(sf_tables_flav(1), x, mu_table, iflav)
      res(2) = EvalPdfTable_xQf(sf_tables_flav(2), x, mu_table, iflav)
      res(3) = EvalPdfTable_xQf(sf_tables_flav(3), x, mu_table, iflav)
    endif
    
  end function StrFct_Flav

  !> @brief calculate the leading order structure function at x, Q, muR, muF 
  !!
  !! Calculate the leading order structure function at x, Q, muR,
  !! muF. muR and muF are only needed if we are using the
  !! scale_choice_arbitrary, as otherwise they are already included in
  !! the sf_tables.
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @return          an array of all the leading order structure functions
  !!
  function F_LO (x, Q, muR, muF) result(res)
    real(dp), intent(in)  :: x, Q, muR, muF
    real(dp) :: res(-6:7)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if(.not.use_sep_orders) then 
      call wae_error('F_LO: you did not initialise the Structure Functions with separate orders. Exiting.')
    endif
    
    call EvalPdfTable_xQ(sf_tables(1), x, mu_table, res)
    
  end function F_LO

  !> @brief calculate the leading order structure function at x, Q, muR, muF 
  !!
  !! Calculate the leading order structure function at x, Q, muR,
  !! muF. muR and muF are only needed if we are using the
  !! scale_choice_arbitrary, as otherwise they are already included in
  !! the sf_tables.
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @param[in]       iflav      parton-flavour
  !! @return          an array of FL, F2, F3 for the given flavour without couplings to the parton
  !!
  function F_LO_flav (x, Q, muR, muF, iflav) result(res)
    real(dp), intent(in)  :: x, Q, muR, muF
    integer, intent(in) :: iflav
    real(dp) :: res(1:3)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if(.not.use_sep_orders) call wae_error('F_LO_flav: you did not   &
         &       initialise the Structure Functions with separate&
         & orders.                   Exiting.')
    
    if(.not.inc_flavour_decomposition) call wae_error('F_LO_flav',&
         & 'You          did not initialise the Structure Functions&
         & with flavour          decomposition. Exiting.')

    res(1) = EvalPdfTable_xQf(sf_tables_flav(1), x, mu_table, iflav)
    res(2) = EvalPdfTable_xQf(sf_tables_flav(2), x, mu_table, iflav)
    res(3) = EvalPdfTable_xQf(sf_tables_flav(3), x, mu_table, iflav)
    
  end function F_LO_FLAV

  !> @brief calculate the NLO structure function at x, Q, muR, muF 
  !!
  !! Calculate the pure NLO contribution to the structure function at x, Q, muR, muF. muR and
  !! muF are only needed if we are using the scale_choice_arbitrary,
  !! as otherwise they are already included in the sf_tables.
  !!
  !! The result includes a factor of (as/2pi)
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @return          an array of all the NLO structure functions
  !!
  function F_NLO (x, Q, muR, muF) result(res)
    real(dp), intent(in)  :: x, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if(.not.use_sep_orders) call wae_error('F_NLO', &
               'You did not initialise the Structure Functions with separate orders. Exiting.')
    
    if (order_setup < 2) then
      res = zero
      return
    end if

    as2pi = alphasLocal(muR) / (twopi)
    
    ! C_NLO x f (x) in C1f(:) 
    call EvalPdfTable_xQ(sf_tables(2), x, mu_table, C1f)
    res = C1f
    
    ! if scale_choice = 0,1 then this term is already taken care of
    if (scale_choice_save.ge.scale_choice_arbitrary) then
       LFQ2 = two*log(muF/Q)
       ! C_LO x P_LO x f (x) in C0P0f(:)
       call EvalPdfTable_xQ(sf_tables(3), x, mu_table, C0P0f)
       res = res - C0P0f * LFQ2
    endif
    
    res = res * as2pi
    
  end function F_NLO

  !> @brief calculate the NLO structure function at x, Q, muR, muF 
  !!
  !! Calculate the pure NLO contribution to the structure function at x, Q, muR, muF. muR and
  !! muF are only needed if we are using the scale_choice_arbitrary,
  !! as otherwise they are already included in the sf_tables.
  !!
  !! The result includes a factor of (as/2pi)
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @param[in]       iflav      parton-flavour
  !! @return          an array of FL, F2, F3 for the given flavour without couplings to the parton
  !!
  function F_NLO_flav (x, Q, muR, muF, iflav) result(res)
    real(dp), intent(in)  :: x, Q, muR, muF
    integer, intent(in) :: iflav
    real(dp) :: res(1:3), as2pi, LFQ2
    real(dp) :: C1f(1:3), C0P0f(1:3)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if(.not.use_sep_orders) call wae_error('F_NLO_flav', 'You did not&
         & initialise the Structure Functions with separate orders.&
         & Exiting.')
    
    if(.not.inc_flavour_decomposition) call wae_error('F_NLO_flav',&
         & 'You did not initialise the Structure Functions with&
         & flavour decomposition. Exiting.')
    
    if (order_setup < 2) then
      res = zero
      return
    end if

    as2pi = alphasLocal(muR) / (twopi)
    
    ! C_NLO x f (x) in C1f(:)
    res(1) = EvalPdfTable_xQf(sf_tables_flav(4), x, mu_table, iflav)
    res(2) = EvalPdfTable_xQf(sf_tables_flav(5), x, mu_table, iflav)
    res(3) = EvalPdfTable_xQf(sf_tables_flav(6), x, mu_table, iflav)
    
    ! if scale_choice = 0,1 then this term is already taken care of
    if (scale_choice_save.ge.scale_choice_arbitrary) then
       LFQ2 = two*log(muF/Q)
       ! C_LO x P_LO x f (x) in C0P0f(:)
       C0P0f(1) = EvalPdfTable_xQf(sf_tables_flav(7), x, mu_table, iflav)
       C0P0f(2) = EvalPdfTable_xQf(sf_tables_flav(8), x, mu_table, iflav)
       C0P0f(3) = EvalPdfTable_xQf(sf_tables_flav(9), x, mu_table, iflav)
       res = res - C0P0f * LFQ2
    endif
    
    res = res * as2pi
    
  end function F_NLO_flav


  !> @brief calculate the NNLO structure function at x, Q, muR, muF 
  !!
  !! Calculate the pure NNLO contribution to the structure function at x, Q, muR, muF. muR and
  !! muF are only needed if we are using the scale_choice_arbitrary,
  !! as otherwise they are already included in the sf_tables.
  !!
  !! The result includes a factor of (as/2pi)**2
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @return          an array of all the NNLO structure functions
  !!
  function F_NNLO (x, Q, muR, muF) result(res)
    real(dp), intent(in)  :: x, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7), C2f(-6:7), C0P0sqf(-6:7)
    real(dp) :: C0P1f(-6:7), C1P0f(-6:7)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)

    if(.not.use_sep_orders) call wae_error('F_NNLO',&
            'You did not initialise the Structure Functions with separate orders. Exiting.')
        
    if (order_setup < 3) then
      res = zero
      return
    end if

    as2pi = alphasLocal(muR) / (twopi)
    
    ! C_NNLO x f (x) in C2f(:,3) 
    call EvalPdfTable_xQ(sf_tables(4), x, mu_table, C2f)
    res = C2f

    ! if scale_choice = 0,1 then these terms are already taken care of
    if (scale_choice_save.ge.scale_choice_arbitrary) then
       LRQ2 = two*log(muR/Q)
       LFQ2 = two*log(muF/Q)
       ! C_NLO x f (x) in C1f
       call EvalPdfTable_xQ(sf_tables(2), x, mu_table, C1f)
       ! C_LO x P_LO x f (x) in C0P0f
       call EvalPdfTable_xQ(sf_tables(3), x, mu_table, C0P0f)
       ! C_LO x P_LO^2 x f (x) in C0P0sqf
       call EvalPdfTable_xQ(sf_tables(5), x, mu_table, C0P0sqf)
       ! C_LO x P_NLO x f (x) in C0P1f
       call EvalPdfTable_xQ(sf_tables(6), x, mu_table, C0P1f)
       ! C_NLO x P_LO x f (x) in C1P1f
       call EvalPdfTable_xQ(sf_tables(7), x, mu_table, C1P0f)
       ! add up all the different pieces
       res = res - C1P0f * LFQ2 + C0P0sqf * half * LFQ2**2 &
            & - twopi * beta0 * C0P0f * (LRQ2*LFQ2 - half*LFQ2**2) &
            & - C0P1f * LFQ2 + twopi * beta0 * C1f * LRQ2
    endif
    
    res = res * (as2pi)**2
    
  end function F_NNLO

  
  !> @brief calculate the pure N3LO contribution to the structure function at x, Q, muR, muF 
  !!
  !! Calculate the pure N3LO contribution to the structure function at x, Q, muR, muF. 
  !! muR and muF are only needed if we are using the scale_choice_arbitrary,
  !! as otherwise they are already included in the sf_tables.
  !!
  !! The result includes a factor of (as/2pi)**3
  !!
  !! @param[in]       x          x value
  !! @param[in]       Q          Q value
  !! @param[in]       muR        renormalisation scale 
  !! @param[in]       muF        factorisation scale
  !! @return          an array of all the N3LO structure functions
  !!
  function F_N3LO (x, Q, muR, muF) result(res)
    real(dp), intent(in)  :: x, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7), C3f(-6:7), C0P0sqf(-6:7), C2f(-6:7)
    real(dp) :: C0P1f(-6:7), C1P0f(-6:7), C0P0cbf(-6:7), C0P01f(-6:7), C0P10f(-6:7)
    real(dp) :: C1P0sqf(-6:7), C1P1f(-6:7), C2P0f(-6:7), C0P2f(-6:7)
    real(dp) :: muR_lcl, muF_lcl, mu_table

    call GetStrFctScales(Q, muR, muF, muR_lcl, muF_lcl, mu_table)
    
    if(.not.use_sep_orders) call wae_error('F_N3LO',&
          'You did not initialise the Structure Functions with separate orders. Exiting.')
    
    if (order_setup < 4) then
      res = zero
      return
    end if

    as2pi = alphasLocal(muR) / (twopi)

    ! C_N3LO x f (x) in C2f(:,8) 
    call EvalPdfTable_xQ(sf_tables(8), x, mu_table, C3f)
    res = C3f

    ! if scale_choice = 0,1 then these terms are already taken care of
    if (scale_choice_save.ge.scale_choice_arbitrary) then
       LRQ2 = two*log(muR/Q)
       LFQ2 = two*log(muF/Q)
       ! C_NLO x f (x) in C1f
       call EvalPdfTable_xQ(sf_tables(2), x, mu_table, C1f)
       ! C_LO x P_LO x f (x) in C0P0f
       call EvalPdfTable_xQ(sf_tables(3), x, mu_table, C0P0f)
       ! C_NNLO x f (x) in C2f(:,3) 
       call EvalPdfTable_xQ(sf_tables(4), x, mu_table, C2f)
       ! C_LO x P_LO^2 x f (x) in C0P0sqf
       call EvalPdfTable_xQ(sf_tables(5), x, mu_table, C0P0sqf)
       ! C_LO x P_NLO x f (x) in C0P1f
       call EvalPdfTable_xQ(sf_tables(6), x, mu_table, C0P1f)
       ! C_NLO x P_LO x f (x) in C1P1f
       call EvalPdfTable_xQ(sf_tables(7), x, mu_table, C1P0f)
       
       ! C_LO x P_LO^3 x f (x) in C0P0cbf
       call EvalPdfTable_xQ(sf_tables(9), x, mu_table, C0P0cbf)
       ! C_LO x P_LO x P_NLO x f (x) in C0P01f
       call EvalPdfTable_xQ(sf_tables(10), x, mu_table, C0P01f)
       ! C_LO x P_LO x P_NLO x f (x) in C0P10f
       call EvalPdfTable_xQ(sf_tables(11), x, mu_table, C0P10f)
       ! C_NLO x P_LO^2 x f (x) in C1P0sqf
       call EvalPdfTable_xQ(sf_tables(12), x, mu_table, C1P0sqf)
       ! C_NLO x P_NLO x f (x) in C1P1f
       call EvalPdfTable_xQ(sf_tables(13), x, mu_table, C1P1f)
       ! C_NNLO x P_LO x f (x) in C2P0f
       call EvalPdfTable_xQ(sf_tables(14), x, mu_table, C2P0f)
       ! C_LO x P_NNLO x f (x) in C0P2f
       call EvalPdfTable_xQ(sf_tables(15), x, mu_table, C0P2f)
    
       ! add up all the different pieces
       ! The commented lines are copy/pasted from mathematica
       ! CpN3LO fmuF + b1 CpNLO fmuF LRQ + 2 b0 CpNNLO fmuF LRQ 
       res = res + twopi**2 * beta1 * C1f * LRQ2 + two * twopi * beta0 * C2f * LRQ2 &
       ! + b0^2 CpNLO fmuF LRQ^2 - CpNNLO fmuF LFQ PLO 
           & + (twopi*beta0)**2 * C1f * LRQ2**2 - C2P0f * LFQ2 &
       ! + 1/2 b1 CpLO fmuF LFQ^2 PLO + 1/2 b0 CpNLO fmuF LFQ^2 PLO  
           & + half * twopi**2 * beta1 * C0P0f * LFQ2**2 + half * twopi * beta0 * C1P0f * LFQ2**2 &
       ! - 1/3 b0^2 CpLO fmuF LFQ^3 PLO - b1 CpLO fmuF LFQ LRQ PLO 
           & - (one/three)*(twopi*beta0)**2 * C0P0f * LFQ2**3 - twopi**2 * beta1 * C0P0f * LFQ2 * LRQ2 &
       ! - 2 b0 CpNLO fmuF LFQ LRQ PLO + b0^2 CpLO fmuF LFQ^2 LRQ PLO 
           & - two * twopi * beta0 * C1P0f * LFQ2 * LRQ2 + (twopi*beta0)**2 * C0P0f * LFQ2**2 * LRQ2 &
       ! - b0^2 CpLO fmuF LFQ LRQ^2 PLO + 1/2 CpNLO fmuF LFQ^2 PLO^2 
           & - (twopi*beta0)**2 * C0P0f * LFQ2 * LRQ2**2 + half * C1P0sqf * LFQ2**2 &
       ! - 1/2 b0 CpLO fmuF LFQ^3 PLO^2 + b0 CpLO fmuF LFQ^2 LRQ PLO^2  
           & - half * twopi*beta0 * C0P0sqf * LFQ2**3 + twopi*beta0 * C0P0sqf * LRQ2 * LFQ2**2 &
       ! - 1/6 CpLO fmuF LFQ^3 PLO^3 - CpNLO fmuF LFQ PNLO 
           & - (one/6.0_dp) * C0P0cbf * LFQ2**3 - C1P1f * LFQ2 &
       ! + b0 CpLO fmuF LFQ^2 PNLO - 2 b0 CpLO fmuF LFQ LRQ PNLO
           & + twopi*beta0 * C0P1f * LFQ2**2 - two*twopi*beta0 * C0P1f * LFQ2 * LRQ2 &
       ! + CpLO fmuF LFQ^2 PLO PNLO - CpLO fmuF LFQ PNNLO
           & + half*(C0P01f+C0P10f) * LFQ2**2 - C0P2f * LFQ2

    endif

    res = res * (as2pi)**3
    
  end function F_N3LO

  !----------------------------------------------------------------------
  ! set_scale_logs is only used for scale_choice = 0,1
  subroutine set_scale_logs(Q)
    real(dp), intent(in) :: Q

    log_muF2_over_Q2 = two * log(sf_muF(Q) / Q)
    log_muR2_over_Q2 = two * log(sf_muR(Q) / Q)
  end subroutine set_scale_logs
  
  !> @brief returns the renormalisation scale used in the structure function sf_tables at Q
  !!
  !! Returns the renormalisation scale used in the structure function
  !! sf_tables at Q. In the case of scale_choice_arbitrary this function
  !! returns Q, as that is the scale at which the structure functions
  !! have been set up.
  !!
  !! @param[in]       Q          Q value
  !! @return          the renormalisation scale used in the structure function sf_tables at Q
  !!
  real(dp) function sf_muR(Q)
    real(dp), intent(in) :: Q
    sf_muR = zero
    if (scale_choice_save.eq.scale_choice_fixed) then
       ! scale_choice = 0 : muF = xmuF * cst_mu
       sf_muR = xmuR * cst_mu
    elseif (scale_choice_save.eq.scale_choice_Q) then
       ! scale_choice = 1 : mu_R = xmuR * Q
       sf_muR = xmuR * Q
    elseif (scale_choice_save.eq.scale_choice_arbitrary) then
       ! scale_choice = 2 : arbitary, return Q
       sf_muR = Q
    else
       call wae_error('sf_muR(Q)', 'illegal value for scale_choice', intval = scale_choice_save)
    endif
  end function sf_muR

  !> @brief returns the factorisation scale used in the structure function sf_tables at Q
  !!
  !! Returns the factorisation scale used in the structure function
  !! sf_tables at Q. In the case of scale_choice_arbitrary this function
  !! returns Q, as that is the scale at which the structure functions
  !! have been set up.
  !!
  !! @param[in]       Q          Q value
  !! @return          the factorisation scale used in the structure function sf_tables at Q
  !!
  real(dp) function sf_muF(Q)
    real(dp), intent(in) :: Q
    sf_muF = zero
    if (scale_choice_save.eq.scale_choice_fixed) then
       ! scale_choice = 0 : muF = xmuF * cst_mu
       sf_muF = xmuF * cst_mu
    elseif (scale_choice_save.eq.scale_choice_Q) then
       ! scale_choice = 1 : muF = xmuF * Q
       sf_muF = xmuF * Q
    elseif (scale_choice_save.eq.scale_choice_arbitrary) then
       ! scale_choice = 2 : arbitary, return Q
       sf_muF = Q
    else
       call wae_error('sf_muF(Q)', 'illegal value for scale_choice', intval = scale_choice_save)
    endif
  end function sf_muF
  
  subroutine use_vfns (f, Q)
    implicit none
    real(dp), intent(in) :: Q
    real(dp), intent(inout) :: f(0:grid%ny,ncompmin:ncompmax)
    integer :: inf, current_nf, dummy!, nfhi

    if(ch%nflo.lt.3) call wae_error('use_vfns',&
      'Do you really want to run with fewer than 3 light flavours?? ch%nflo=', intval=ch%nflo)

    ! First set current nf equal to the smallest value it can take
    current_nf = ch%nflo
    if(Q.lt.masses(4)) then   
      current_nf = 3
    elseif(Q.lt.masses(5)) then 
      current_nf = 4
    elseif(Q.lt.masses(6)) then  
      current_nf = 5
    elseif(Q.ge.masses(6)) then  
      current_nf = 6
    endif      

    current_nf = min(current_nf,ch%nfhi) ! In case we run we force fewer than 6 flavours

    do inf = current_nf + 1, 6 
      f(:,-inf) = zero
      f(:, inf) = zero
    enddo
    
    nf_lcl = current_nf

    ! set nf in coeff fct to current number of flavours above threshold
    call SetNfCoefHolder(ch, current_nf)
    
  end subroutine use_vfns

  
  !----------------------------------------------------------------------
  ! Rescale the PDF by F2 NNLO / F2 N3LO
  subroutine rescale_pdf_n3lo (Qval, order)
    real(dp), intent(in) :: Qval
    integer, intent(in)  :: order
    real(dp) :: mF, mR, factor
    real(dp) :: str_fct(-6:7), f2nnlo, f2n3lo, y(0:grid%ny)
    integer :: iy, ii
    mF = sf_muF(Qval)
    mR = sf_muR(Qval)
    y = yValues(grid)
    do iy = 0, grid%ny
       str_fct(:) = F_LO(y(iy), Qval, mR, mF) + F_NLO(y(iy), Qval, mR, mF) + F_NNLO(y(iy), Qval, mR, mF)
       f2nnlo = str_fct(iF2Z)
       str_fct(:) = str_fct(:) + F_N3LO(y(iy), Qval, mR, mF)
       f2n3lo = str_fct(iF2Z)
       factor = 1.0_dp
       if (f2n3lo.gt.0.0_dp) factor = f2nnlo/f2n3lo
       do ii = ncompmin, ncompmax
          tables(0)%tab(iy,ii,:) = tables(0)%tab(iy,ii,:) * factor
       enddo
       
    enddo

    ! Re-set the structure functions using updated PDF
    if (scale_choice_save.eq.scale_choice_Q.or.scale_choice_save.eq.scale_choice_fixed) then
       if (order.ge.1) call set_LO_structure_functions()
       if (order.ge.2) call set_NLO_structure_functions()
       if (order.ge.3) call set_NNLO_structure_functions()
       if (order.ge.4) call set_N3LO_structure_functions()
    else if(scale_choice_save.eq.scale_choice_arbitrary) then
       if (order.ge.1) call set_LO_structure_functions_anyscale()
       if (order.ge.2) call set_NLO_structure_functions_anyscale()
       if (order.ge.3) call set_NNLO_structure_functions_anyscale()
       if (order.ge.4) call set_N3LO_structure_functions_anyscale()
    else
       call wae_error('rescale_pdf_n3lo', 'illegal value for scale_choice', intval = scale_choice_save)
    endif

  end subroutine rescale_pdf_n3lo
  
  !---------------------------------------------------------------------- 
  ! Rescale the PDF by F2 NLO / F2 NNLO
  subroutine rescale_pdf_nnlo (Qval, order)
    real(dp), intent(in) :: Qval
    integer, intent(in)  :: order
    real(dp) :: mF, mR, factor
    real(dp) :: str_fct(-6:7), f2nlo, f2nnlo, y(0:grid%ny) 
    integer :: iy, ii
    mF = sf_muF(Qval) 
    mR = sf_muR(Qval)
    y = yValues(grid)
    do iy = 0, grid%ny              
       str_fct(:) = F_LO(y(iy), Qval, mR, mF) + F_NLO(y(iy), Qval, mR, mF)
       f2nlo = str_fct(iF2Z) 
       str_fct(:) = str_fct(:) + F_NNLO(y(iy), Qval, mR, mF)
       f2nnlo = str_fct(iF2Z)
       factor = 1.0_dp
       if (f2nnlo.gt.0.0_dp) factor = f2nlo/f2nnlo
       do ii = ncompmin, ncompmax
          tables(0)%tab(iy,ii,:) = tables(0)%tab(iy,ii,:) * factor
       enddo
    enddo
    
    ! Re-set the structure functions using updated PDF
    if (scale_choice_save.eq.scale_choice_Q.or.scale_choice_save.eq.scale_choice_fixed) then
       if (order.ge.1) call set_LO_structure_functions()
       if (order.ge.2) call set_NLO_structure_functions()
       if (order.ge.3) call set_NNLO_structure_functions()
       if (order.ge.4) call set_N3LO_structure_functions()
    else if(scale_choice_save.eq.scale_choice_arbitrary) then
       if (order.ge.1) call set_LO_structure_functions_anyscale()
       if (order.ge.2) call set_NLO_structure_functions_anyscale()
       if (order.ge.3) call set_NNLO_structure_functions_anyscale()
       if (order.ge.4) call set_N3LO_structure_functions_anyscale()
    else
       call wae_error('rescale_pdf_n3lo', 'illegal value for scale_choice', intval = scale_choice_save)
    endif
  end subroutine rescale_pdf_nnlo


end module structure_functions


!> @brief  Minimal setup of structure functions
!!
!! @param[in]      order_max      highest order in QCD to compute (1: LO, 2: NLO, etc)
!!
subroutine hoppetStartStrFct(order_max)
  use streamlined_interface; use structure_functions
  implicit none
  integer, intent(in)  :: order_max
  !----------------------------------------------------------------------
  
  call StartStrFct(order_max)
  
end subroutine hoppetStartStrFct

!> @brief  Setup of constants and parameters needed for structure functions
!!
!! @param[in]      order_max      highest order in QCD to compute (1: LO, 2: NLO, etc)
!! @param[opt]     nflav          integer number of flavours (if not present use variable flavour)
!! @param[opt]     scale_choice   (0: fixed scale, 1: use Q, 2: use arbitrary scale)
!! @param[opt]     constant_mu    if scale_choice = scale_choice_fixed (= 0) then this is the fixed scale
!! @param[opt]     param_coefs    if .true. use parametrised coefficients functions
!! @param[opt]     wmass          Mass of the W boson
!! @param[opt]     zmass          Mass of the z boson
!!
subroutine hoppetStartStrFctExtended(order_max, nflav, scale_choice,&
     & constant_mu, param_coefs, wmass, zmass)
  use streamlined_interface; use structure_functions
  implicit none
  real(dp), intent(in) :: constant_mu, wmass, zmass
  integer, intent(in)  :: order_max, nflav, scale_choice
  logical , intent(in) :: param_coefs
  !----------------------------------------------------------------------

call StartStrFct(order_max, nflav, scale_choice, constant_mu,&
     & param_coefs, wmass, zmass)
  
end subroutine hoppetStartStrFctExtended

!> @brief Initialize the structure functions up to specified order
!!
!! Initialize the structure functions up to specified order
!! this requires the tabulated PDF to have been set up beforehand in the
!! the streamlined interface. 
!!
!! By default, separate_orders = .false., which means that one
!! can only access the sum over all orders. If separate_orders = .true.
!! then one can access each order separately, but this is slower.
!!
!! @param[in]      order              order at which we initialise the structure functions
!! @param[opt]     separate_orders    if .true. separate into individual orders rather than summing
!! @param[opt]     xR             factor to multiply renormalisation scale
!! @param[opt]     xF             factor to multiply factorisation scale
!!
subroutine hoppetInitStrFct(order, separate_orders, xR, xF)
  use streamlined_interface; use structure_functions
  implicit none
  integer, intent(in) :: order
  logical, intent(in) :: separate_orders
  real(dp), intent(in) :: xR, xF
  
  call InitStrFct(order, separate_orders, xR, xF, .false.)

end subroutine hoppetInitStrFct

!> @brief Initialize the structure functions up to specified order
!!
!! Initialize the structure functions up to specified order
!! this requires the tabulated PDF to have been set up beforehand in the
!! the streamlined interface. 
!!
!! By default, separate_orders = .false., which means that one
!! can only access the sum over all orders. If separate_orders = .true.
!! then one can access each order separately, but this is slower.
!!
!! @param[in]      order              order at which we initialise the structure functions
!! @param[opt]     separate_orders    if .true. separate into individual orders rather than summing
!! @param[opt]     xR             factor to multiply renormalisation scale
!! @param[opt]     xF             factor to multiply factorisation scale
!!
subroutine hoppetInitStrFctFlav(order, separate_orders, xR, xF, flavour_decomposition)
  use streamlined_interface; use structure_functions
  implicit none
  integer, intent(in) :: order
  logical, intent(in) :: separate_orders, flavour_decomposition
  real(dp), intent(in) :: xR, xF
  
  call InitStrFct(order, separate_orders, xR, xF, flavour_decomposition)

end subroutine hoppetInitStrFctFlav

!> @brief calculate the structure function at x, Q, muR, muF summed over all orders
!!
!! Calculate the structure function at x, Q, muR, muF summed over
!! all orders. If using a scale choice other than scale_choice_arbitrary,
!! muR and muF must be consistent with the scale choice made in the 
!! hoppetStartStrFct call. See also hoppetStrFctNoMu for a variant
!! that does not take muR and muF arguments.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @return          an array of all structure functions summed over orders
!!
subroutine hoppetStrFct(x, Q, muR_in, muF_in, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  real(dp) :: res(-6:7)

  res = StrFct(x, Q, muR_in, muF_in)
end subroutine hoppetStrFct

!> @brief calculate the structure function at x, Q, muR, muF summed over all orders
!!
!! Calculate the structure function at x, Q, muR, muF summed over
!! all orders. If using a scale choice other than scale_choice_arbitrary,
!! muR and muF must be consistent with the scale choice made in the 
!! hoppetStartStrFct call. See also hoppetStrFctNoMu for a variant
!! that does not take muR and muF arguments.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @param[in]       iflav      parton flavour
!! @return          an array of all structure functions summed over orders decomposed in flavour
!!
subroutine hoppetStrFctFlav(x, Q, muR_in, muF_in, iflav, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  integer  :: iflav
  real(dp) :: res(1:3)

  res = StrFct_flav(x, Q, muR_in, muF_in, iflav)
end subroutine hoppetStrFctFlav

!> @brief calculate the structure function at x, Q, summed over all orders
!!
!! Calculate the structure function at x, Q, for the scale choice as indicated
!! in hoppetStartStrFct. This can only be used with scale_choice_Q and scale_choice_fixed. 
!! See also hoppetStrFct for a variant with muR and muF choice. 
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @return          an array of all structure functions summed over orders
!!
subroutine hoppetStrFctNoMu(x, Q, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: res(-6:7)
  
  res = StrFct(x, Q)
end subroutine hoppetStrFctNoMu

!> @brief calculate the structure function at x, Q, summed over all orders
!!
!! Calculate the structure function at x, Q, for the scale choice as indicated
!! in hoppetStartStrFct. This can only be used with scale_choice_Q and scale_choice_fixed. 
!! See also hoppetStrFct for a variant with muR and muF choice. 
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       iflav      parton flavour
!! @return          an array of all structure functions summed over orders decomposed in flavour
!!
subroutine hoppetStrFctNoMuFlav(x, Q, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: res(1:3)
  
  res = StrFct_flav(x, Q, iflav = iflav)
end subroutine hoppetStrFctNoMuFlav


!> @brief calculate the leading order structure function at x, Q, muR, muF 
!!
!! Calculate the leading order structure function at x, Q, muR,
!! muF. muR and muF are only needed if we are using the
!! scale_choice_arbitrary, as otherwise they are already included in
!! the sf_tables.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @return          an array of all the leading order structure functions
!!
subroutine hoppetStrFctLO(x, Q, muR_in, muF_in, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  real(dp) :: res(-6:7)

  res = F_LO(x, Q, muR_in, muF_in)
end subroutine hoppetStrFctLO

!> @brief calculate the leading order structure function at x, Q, muR, muF 
!!
!! Calculate the leading order structure function at x, Q, muR,
!! muF. muR and muF are only needed if we are using the
!! scale_choice_arbitrary, as otherwise they are already included in
!! the sf_tables.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @param[in]       iflav      partonflavour
!! @return          an array of all the leading order structure functions summed over orders
!!
subroutine hoppetStrFctLOFlav(x, Q, muR_in, muF_in, iflav, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  integer  :: iflav
  real(dp) :: res(1:3)

  res = F_LO_flav(x, Q, muR_in, muF_in, iflav)
end subroutine hoppetStrFctLOFlav

!> @brief calculate the NLO structure function at x, Q, muR, muF 
!!
!! Calculate the NLO structure function at x, Q, muR, muF. muR and
!! muF are only needed if we are using the scale_choice_arbitrary,
!! as otherwise they are already included in the sf_tables.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @return          an array of all the NLO structure functions
!!
subroutine hoppetStrFctNLO(x, Q, muR_in, muF_in, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  real(dp) :: res(-6:7)

  res = F_NLO(x, Q, muR_in, muF_in)
end subroutine hoppetStrFctNLO

!> @brief calculate the NLO structure function at x, Q, muR, muF 
!!
!! Calculate the NLO structure function at x, Q, muR, muF. muR and
!! muF are only needed if we are using the scale_choice_arbitrary,
!! as otherwise they are already included in the sf_tables.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @return          an array of all the NLO structure functions
!!
subroutine hoppetStrFctNLOFlav(x, Q, muR_in, muF_in, iflav, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  integer  :: iflav
  real(dp) :: res(1:3)

  res = F_NLO_flav(x, Q, muR_in, muF_in, iflav)
end subroutine hoppetStrFctNLOFlav

!> @brief calculate the NNLO structure function at x, Q, muR, muF 
!!
!! Calculate the NNLO structure function at x, Q, muR, muF. muR and
!! muF are only needed if we are using the scale_choice_arbitrary,
!! as otherwise they are already included in the sf_tables.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @return          an array of all the NNLO structure functions
!!
subroutine hoppetStrFctNNLO(x, Q, muR_in, muF_in, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  real(dp) :: res(-6:7)

  res = F_NNLO(x, Q, muR_in, muF_in)
end subroutine hoppetStrFctNNLO

!> @brief calculate the N3LO structure function at x, Q, muR, muF 
!!
!! Calculate the N3LO structure function at x, Q, muR, muF. muR and
!! muF are only needed if we are using the scale_choice_arbitrary,
!! as otherwise they are already included in the sf_tables.
!!
!! @param[in]       x          x value
!! @param[in]       Q          Q value
!! @param[in]       muR        renormalisation scale 
!! @param[in]       muF        factorisation scale
!! @return          an array of all the N3LO structure functions
!!
subroutine hoppetStrFctN3LO(x, Q, muR_in, muF_in, res) 
  use streamlined_interface; use structure_functions
  real(dp) :: x, Q
  real(dp) :: muR_in, muF_in
  real(dp) :: res(-6:7)

  res = F_N3LO(x, Q, muR_in, muF_in)
end subroutine hoppetStrFctN3LO
