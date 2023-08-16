!======================================================================
!
! Summary of our understanding of the coefficient functions
! ---------------------------------------------------------
!
! - F1 = (F2 - FL)/2x
!
! - hep-ph/0504042 (MMV): gives F2 and FL results (electromagnetic)
!
!   (1/x) F_a = C_{a,ns} \otimes q_{ns}
!                   +  <e^2> (C_{a,q} \otimes q_s  + C_{a,g} \otimes g)
!
!   with z = 2,L, where:
!
!   * <e^2> is the average electric charge
!   * q_{ns} = ???? (we're supposed to deduce it from Eq.(4.2))
!   * C_{a,q} = C_{a,ns} + C_{a,ps};
!   * C_{2,ps}^{(0)} = C_{2,ps}^{(1)} = 0 (and presumably for FL too)
!
! - http://www.sciencedirect.com/science/article/pii/055032139290087R#
!   (Zijlstra & van Neerven) has the second-order coefficient
!   functions. 
!
! - from 1109.3717 (Zaro & Co):
!
!   * q_{ns,i}^+ = (q_i + qbar_i) - q_s
!   * q_{ns,i}^- = (q_i - qbar_i) - q_{ns}^v
!   * q_{ns}^v   = \sum_{i=1}^{n_f} (q_i - qbar_i)
!   * q_s        = \sum_{i=1}^{n_f} (q_i + qbar_i)
!
!   That, together with [C_{a,q} = C_{a,ns} + C_{a,ps}] means that the combination
! 
!       q_{ns,j}^+ * C_{i,ns}^+ + q_s * C_{i,q}
!
!   reduces to 
!
!       (q_j+qbar_j) * C_{i,ns}^+ + q_s * C_{i,ps}
!
module structure_functions
  use pdf_representation
  use streamlined_interface
  use splitting_functions
  use coefficient_functions_holder
  use qcd
  implicit none

  private
  public :: InitStrFct, StartStrFct, write_f1, write_f2, write_f3
  public :: F1Wp, F2Wp, F3Wp, F1Wm, F2Wm, F3Wm, F1Z, F2Z, F3Z, F1EM, F2EM, F1gZ, F2gZ, F3gZ
  public :: F_LO, F_NLO, F_NNLO, F_N3LO
  public :: muR, muF, quark_masses_sf
  
  ! holds the coefficient functions
  type(coef_holder), save :: ch
  
  ! indices for the different structure functions
  integer, parameter :: F1Wp= 1, F2Wp= 2, F3Wp= 3
  integer, parameter :: F1Wm=-1, F2Wm=-2, F3Wm=-3
  integer, parameter :: F1Z = 4, F2Z = 5, F3Z = 6
  integer, parameter :: F1EM = -4, F2EM = -5
  integer, parameter :: F1gZ = 0, F2gZ = -6, F3gZ = 7
  
  ! constants and fixed parameters
  real(dp), save      :: quark_masses_sf(4:6) = quark_masses_def(4:6)
  logical,  save      :: use_mass_thresholds
  real(dp), parameter :: viW = 1/sqrt(two), aiW = viW
  real(dp), save      :: vi2_ai2_avg_W, two_vi_ai_avg_W
  real(dp), save      :: vi2_ai2_Z_down, vi2_ai2_Z_up
  real(dp), save      :: two_vi_ai_Z_down, two_vi_ai_Z_up
  real(dp), save      :: two_vi_Z_down, two_vi_Z_up
  real(dp), save      :: two_ai_Z_down, two_ai_Z_up
  real(dp), parameter :: e_up = 2.0_dp/3.0_dp, e_down = - 1.0_dp/3.0_dp
  real(dp), parameter :: e2_up = 4.0_dp/9.0_dp, e2_down = 1.0_dp/9.0_dp
  ! these log terms are only used for scale_choice = 0,1, to speed up the code
  real(dp), save      :: log_muF2_over_Q2, log_muR2_over_Q2
  real(dp), save      :: Qmin
  real(dp), public, save :: toy_alphas_Q0 = 0.35_dp ! Roughly αS(sqrt(2))
  type(running_coupling), public, save :: toy_coupling

  ! scale choices
  real(dp), save :: cst_mu, xmuR, xmuF
  integer, save  :: scale_choice

  logical, save  :: use_sep_orders
  logical, save  :: exact_coefs
  integer :: nf_lcl, order_setup


contains

  !----------------------------------------------------------------------
  ! Setup of constants and parameters needed for structure functions
  subroutine StartStrFct(rts, order_max, nflav, xR, xF, sc_choice, cmu, param_coefs, &
       & Qmin_PDF, wmass, zmass)
    real(dp), intent(in) :: rts
    integer, intent(in)  :: order_max
    integer, optional    :: nflav, sc_choice
    real(dp), optional   :: xR, xF, cmu, Qmin_PDF, wmass, zmass
    logical, optional    :: param_coefs
    !----------------------------------------------------------------------
    real(dp) :: sin_thw_sq, ymax, dy, minQval, maxQval, dlnlnQ, mw, mz
    integer  :: nloop, order

    ! take sensible default value for mw and mz
    mw = 80.398_dp
    mz = 91.187_dp
    ! otherwise use user input
    if(present(wmass)) mw = wmass
    if(present(zmass)) mz = zmass
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
    
    ! default settings
    xmuR = one
    xmuF = one
    scale_choice = 1
    exact_coefs  = .false.
    cst_mu       = zero
    Qmin         = one
    
    ! change to user input if specified
    if(present(xR)) xmuR=xR
    if(present(xF)) xmuF=xF
    if(present(sc_choice)) scale_choice=sc_choice
    if(present(cmu)) cst_mu=cmu
    if(present(param_coefs)) exact_coefs = .not.param_coefs
    if(present(Qmin_PDF)) Qmin=Qmin_PDF

    ! Streamlined initialization
    ! including  parameters for x-grid
    order = -6 
    ymax  = 16.0_dp
    dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0_dp
    nloop = 3
    minQval = min(xmuF*Qmin, Qmin)
    maxQval = max(xmuF*rts, rts)

    ! initialise the grid and dglap holder
    call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)

    if (present(nflav)) then
       ! if nflav is present, then use a fixed flavour number
       call qcd_SetNf(nflav)
       call InitCoefHolder(grid, ch, order_max, exact_coefs)
       use_mass_thresholds = .false.
       nf_lcl = nf_int
       print*, "Starting the structure functions with fixed number of flavours"
       write(6,*) "nf = ", nf_lcl
    else
       ! otherwise, use variable flavour number scheme
   !    call InitCoefHolder(grid, ch, order_max, exact_coefs, nflo=3, nfhi=5)
       call InitCoefHolder(grid, ch, order_max, exact_coefs, nflo=3, nfhi=6)
       use_mass_thresholds = .true.
       
       !call InitCoefHolder(grid, ch, order_max, exact_coefs)
       !use_mass_thresholds = .false.
       quark_masses_sf(4:6) = masses(4:6) ! Takes the thresholds from the streamlined 
                                          ! interface which should be filled first!
       !quark_masses_sf(5) = 4.18_dp
       print*, "Starting the structure functions with mass thresholds at"
       print*, "mc = ", quark_masses_sf(4) 
       print*, "mb = ", quark_masses_sf(5) 
       print*, "mt = ", quark_masses_sf(6) 
       ! and start with a sensible local nf (which will be 5 here) 
       nf_lcl = nf_int
    endif
    
    ! if mass thresholds are used in the structure functions, only scale_choice 0 and 1 are allowed
    if ((use_mass_thresholds).and.(scale_choice.gt.1)) then
       call wae_error('StartStrFct', 'illegal value for scale_choice with mass thresholds turned on', intval = scale_choice)
    end if

    ! AK: Finally we need to set tab_iflv_max = 7. We don't change tables(0) as it contains the PDF
    tables(1:)%tab_iflv_max = 7

  end subroutine StartStrFct

  
  !----------------------------------------------------------------------
  ! Initialize the structure functions up to specified order
  ! this requires the PDF to have been set up beforehand, and filled in tables(0)
  subroutine InitStrFct(order, separate_orders)
    integer, intent(in) :: order
    logical, optional :: separate_orders

    ! To turn off b quarks completely (only for testing and comparison)
    ! uncomment following two lines:
    ! tables(0)%tab(:,-5,:) = zero
    ! tables(0)%tab(:,+5,:) = zero

    ! default is to sum each order before convolutions (it's faster)
    if (present(separate_orders)) then
       use_sep_orders = separate_orders
    else
       use_sep_orders = .false.
    endif

    if (use_sep_orders) then
       ! First we treat the case where we want to separate out each order in
       ! the final structure functions.
       ! This is slower, but is needed e.g. for VBFH
       if (scale_choice.le.1) then
          ! if scale_choice = 0,1 use fast implementation
          ! tables is saved as an array in Q, and only tables(0),
          ! tables(1), tables(2), tables(4), tables(8) are non zero
          if (order.ge.1) call set_LO_structure_functions()
          if (order.ge.2) call set_NLO_structure_functions()
          if (order.ge.3) call set_NNLO_structure_functions()
          if (order.ge.4) call set_N3LO_structure_functions()
       else
          ! if scale_choice >= 2 use slower implementation with full
          ! scale choices, such as sqrt(Q1*Q2)
          ! tables is saved as an array in muF now, and all components
          ! of tables are non zero.
          if (order.ge.1) call set_LO_structure_functions_anyscale()
          if (order.ge.2) call set_NLO_structure_functions_anyscale()
          if (order.ge.3) call set_NNLO_structure_functions_anyscale()
          if (order.ge.4) call set_N3LO_structure_functions_anyscale()
       endif
    else
       ! Now set up the default case where we sum up everything right away.
       if (scale_choice.gt.1) & ! only allow for scale_choice = 0,1
            & call wae_error('InitStrFct', 'illegal value for scale_choice', intval = scale_choice)
       call set_structure_functions_upto(order)
    endif

    ! To rescale PDF by the N3LO F2 structure function (evaluated at 8 GeV)
    ! as a check of the size of missing N3LO PDFs, uncomment following line:
    ! call rescale_pdf_nnlo(8.0_dp,order)
    ! call rescale_pdf_n3lo(8.0_dp,order)

    setup_done(1:) = .true.
    order_setup = order
    
  end subroutine InitStrFct

  !----------------------------------------------------------------------
  ! Rescale the PDF by F2 NNLO / F2 N3LO
  subroutine rescale_pdf_n3lo (Qval, order)
    real(dp), intent(in) :: Qval
    integer, intent(in)  :: order
    real(dp) :: mF, mR, factor
    real(dp) :: str_fct(-6:7), f2nnlo, f2n3lo, y(0:grid%ny)
    integer :: iy, ii
    mF = muF(Qval)
    mR = muR(Qval)
    y = yValues(grid)
    do iy = 0, grid%ny
       str_fct(:) = F_LO(y(iy), Qval, mR, mF) + F_NLO(y(iy), Qval, mR, mF) + F_NNLO(y(iy), Qval, mR, mF)
       f2nnlo = str_fct(F2Z)
       str_fct(:) = str_fct(:) + F_N3LO(y(iy), Qval, mR, mF)
       f2n3lo = str_fct(F2Z)
       factor = 1.0_dp
       if (f2n3lo.gt.0.0_dp) factor = f2nnlo/f2n3lo
       do ii = ncompmin, ncompmax
          tables(0)%tab(iy,ii,:) = tables(0)%tab(iy,ii,:) * factor
       enddo
       
    enddo

    ! Re-set the structure functions using updated PDF
    if (scale_choice.le.1) then
       if (order.ge.1) call set_LO_structure_functions()
       if (order.ge.2) call set_NLO_structure_functions()
       if (order.ge.3) call set_NNLO_structure_functions()
       if (order.ge.4) call set_N3LO_structure_functions()
    else
       if (order.ge.1) call set_LO_structure_functions_anyscale()
       if (order.ge.2) call set_NLO_structure_functions_anyscale()
       if (order.ge.3) call set_NNLO_structure_functions_anyscale()
       if (order.ge.4) call set_N3LO_structure_functions_anyscale()
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
    mF = muF(Qval) 
    mR = muR(Qval)
    y = yValues(grid)
    do iy = 0, grid%ny              
       str_fct(:) = F_LO(y(iy), Qval, mR, mF) + F_NLO(y(iy), Qval, mR, mF)
       f2nlo = str_fct(F2Z) 
       str_fct(:) = str_fct(:) + F_NNLO(y(iy), Qval, mR, mF)
       f2nnlo = str_fct(F2Z)
       factor = 1.0_dp
       if (f2nnlo.gt.0.0_dp) factor = f2nlo/f2nnlo
       do ii = ncompmin, ncompmax
          tables(0)%tab(iy,ii,:) = tables(0)%tab(iy,ii,:) * factor
       enddo
    enddo
    
    ! Re-set the structure functions using updated PDF
    if (scale_choice.le.1) then
       if (order.ge.1) call set_LO_structure_functions()
       if (order.ge.2) call set_NLO_structure_functions()
       if (order.ge.3) call set_NNLO_structure_functions()
       if (order.ge.4) call set_N3LO_structure_functions()
    else
       if (order.ge.1) call set_LO_structure_functions_anyscale()
       if (order.ge.2) call set_NLO_structure_functions_anyscale()
       if (order.ge.3) call set_NNLO_structure_functions_anyscale()
       if (order.ge.4) call set_N3LO_structure_functions_anyscale()
    endif
  end subroutine rescale_pdf_nnlo
  
  !----------------------------------------------------------------------
  ! write the F1 structure function to idev
  subroutine write_f1 (idev, Qtest, ymax, ny)
    real(dp), intent(in) :: Qtest, ymax
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO, res(-6:7)
    integer  :: iy
    !F1 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    if (use_sep_orders) then
       write(idev,'(a,a)') '# x  F1Wp(LO) F1Wm(LO) F1Wp(NLO) F1Wm(NLO) F1Wp(NNLO) F1Wm(NNLO)', &
            & ' F1Wp(N3LO) F1Wm(N3LO) F1Z(LO) F1Z(NLO) F1Z(NNLO) F1Z(N3LO)'
    else
       write(idev,'(a)') '# x F1Wp F1Wm F1Z'
    endif
    mF = muF(Qtest)
    mR = muR(Qtest)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       if (use_sep_orders) then
          res = F_LO(ytest, Qtest, mR, mF)
          write(idev,'(3es22.12)',advance='no') xval, res(F1Wp),res(F1Wm)
          F1Z_LO = res(F1Z)
          res = F_NLO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F1Wp), res(F1Wm)
          F1Z_NLO = res(F1Z)
          res = F_NNLO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F1Wp), res(F1Wm)
          F1Z_NNLO = res(F1Z)
          res = F_N3LO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F1Wp), res(F1Wm)
          F1Z_N3LO = res(F1Z)
          write(idev,'(4es22.12)',advance='no') F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO
       else
          res = StrFct(ytest, Qtest, mR, mF)
          write(idev,'(4es22.12)',advance='no') xval, res(F1Wp),res(F1Wm), res(F1Z)
       endif
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f1
  
  !----------------------------------------------------------------------
  ! write the F2 structure function to idev
  subroutine write_f2 (idev, Qtest, ymax, ny)
    real(dp), intent(in) :: Qtest, ymax
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO, res(-6:7)
    integer  :: iy
    !F2 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    if (use_sep_orders) then
       write(idev,'(a,a)') '# x  F2Wp(LO) F2Wm(LO) F2Wp(NLO) F2Wm(NLO) F2Wp(NNLO) F2Wm(NNLO)', &
            & ' F2Wp(N3LO) F2Wm(N3LO) F2Z(LO) F2Z(NLO) F2Z(NNLO) F2Z(N3LO)'
    else
       write(idev,'(a)') '# x F2Wp F2Wm F2Z'
    endif
    mF = muF(Qtest)
    mR = muR(Qtest)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       if (use_sep_orders) then
          res = F_LO(ytest, Qtest, mR, mF)
          write(idev,'(3es22.12)',advance='no') xval, res(F2Wp),res(F2Wm)
          F2Z_LO = res(F2Z)
          res = F_NLO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F2Wp), res(F2Wm)
          F2Z_NLO = res(F2Z)
          res = F_NNLO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F2Wp), res(F2Wm)
          F2Z_NNLO = res(F2Z)
          res = F_N3LO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F2Wp), res(F2Wm)
          F2Z_N3LO = res(F2Z)
          write(idev,'(4es22.12)',advance='no') F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO
       else
          res = StrFct(ytest, Qtest, mR, mF)
          write(idev,'(4es22.12)',advance='no') xval, res(F2Wp),res(F2Wm), res(F2Z)
       endif
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f2

  !----------------------------------------------------------------------
  ! write the F3 structure function to idev
  subroutine write_f3 (idev, Qtest, ymax, ny)
    real(dp), intent(in) :: Qtest, ymax
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO, res(-6:7)
    integer  :: iy
    !F3 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    if (use_sep_orders) then
       write(idev,'(a,a)') '# x  F3Wp(LO) F3Wm(LO) F3Wp(NLO) F3Wm(NLO) F3Wp(NNLO) F3Wm(NNLO)', &
            & ' F3Wp(N3LO) F3Wm(N3LO) F3Z(LO) F3Z(NLO) F3Z(NNLO) F3Z(N3LO)'
    else
       write(idev,'(a)') '# x F3Wp F3Wm F3Z'
    endif
    mF = muF(Qtest)
    mR = muR(Qtest)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       if (use_sep_orders) then
          res = F_LO(ytest, Qtest, mR, mF)
          write(idev,'(3es22.12)',advance='no') xval, res(F3Wp),res(F3Wm)
          F3Z_LO = res(F3Z)
          res = F_NLO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F3Wp), res(F3Wm)
          F3Z_NLO = res(F3Z)
          res = F_NNLO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F3Wp), res(F3Wm)
          F3Z_NNLO = res(F3Z)
          res = F_N3LO(ytest, Qtest, mR, mF)
          write(idev,'(2es22.12)',advance='no') res(F3Wp), res(F3Wm)
          F3Z_N3LO = res(F3Z)
          write(idev,'(4es22.12)',advance='no') F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO
       else
          res = StrFct(ytest, Qtest, mR, mF)
          write(idev,'(4es22.12)',advance='no') xval, res(F3Wp),res(F3Wm), res(F3Z)
       endif
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f3

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
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       
       as2pi = alphasLocal(muR(Q)) / (twopi)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
   
       ! start with LO
       if (order.ge.1) then
          tables(1)%tab(:,:,iQ) = structure_function_general(ch%C2LO*f, ch%CLLO*f, ch%C3LO*f)
       endif
       
       ! now add NLO terms
       if (order.ge.2) then
          if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
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

          tables(1)%tab(:,:,iQ) = tables(1)%tab(:,:,iQ) + &
               & as2pi * structure_function_general(f2, fL, f3)
       endif

       ! now add NNLO terms
       if (order.ge.3) then
          if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
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

          tables(1)%tab(:,:,iQ) = tables(1)%tab(:,:,iQ) + &
               (as2pi**2) * structure_function_general(f2, fL, f3)
       endif

       ! and finally, add the N3LO terms
       if (order.ge.4) then
          if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
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

          tables(1)%tab(:,:,iQ) = tables(1)%tab(:,:,iQ) + &
               & (as2pi**3) * structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)
       endif
    end do
  end subroutine set_structure_functions_upto

  
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
       ! explicitly evaluate the PDF at scale muF(Q)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       tables(1)%tab(:,:,iQ) = structure_function_general(ch%C2LO*f, ch%CLLO*f, ch%C3LO*f)
    end do
    
  end subroutine set_LO_structure_functions
  
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
       
       tables(1)%tab(:,:,iQ) = structure_function_general(&
            & f * ch%C2LO, &
            & f * ch%CLLO, &
            & f * ch%C3LO)       
    end do
    
  end subroutine set_LO_structure_functions_anyscale
  
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
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif

       if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
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
          
       tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine set_NLO_structure_functions

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
       ! Internal Q value effectively corresponds to muF(Q1,Q2)
       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)

       ! Save the NLO pieces in tables(2) and tables(3)
      
       ! Get the NLO coefficient function, (C_NLO x f) 
       f2 = (ch%C2NLO * f)
       fL = (ch%CLNLO * f)
       f3 = (ch%C3NLO * f)
       tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now compute the (C_LO x P_LO x f) term
       PLO_f = dh%P_LO * f
       f2 = (ch%C2LO * PLO_f)
       fL = (ch%CLLO * PLO_f)
       f3 = (ch%C3LO * PLO_f)
       tables(3)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
    end do

  end subroutine set_NLO_structure_functions_anyscale
  
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
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       
       if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
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

       tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
     
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
       
       ! save the NNLO pieces in tables(4:7)
       
       PLO2_f  = dh%P_LO  * (dh%P_LO * f)
       PLO_f   = dh%P_LO  * f
       PNLO_f  = dh%P_NLO * f

       ! first calculate the pure NNLO term, (C_NNLO x f) 
       f2 = (ch%C2NNLO * f)
       fL = (ch%CLNNLO * f)
       f3 = (ch%C3NNLO * f)
       tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO^2 x f) term
       f2 =  (ch%C2LO * PLO2_f)
       fL =  (ch%CLLO * PLO2_f)
       f3 =  (ch%C3LO * PLO2_f)
       tables(5)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO) term
       f2 = (ch%C2LO * PNLO_f)
       fL = (ch%CLLO * PNLO_f)
       f3 = (ch%C3LO * PNLO_f)
       tables(6)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO) term
       f2 = (ch%C2NLO * PLO_f)
       fL = (ch%CLNLO * PLO_f)
       f3 = (ch%C3NLO * PLO_f)
       tables(7)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

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
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       
       if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
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

       !tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)
     
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
       
       ! save the N3LO pieces in tables(8:15)
       
       PLO2_f   = dh%P_LO  * (dh%P_LO * f)
       PLO_f    = dh%P_LO  * f
       PNLO_f   = dh%P_NLO * f
       PNNLO_f  = dh%P_NNLO * f
       PLONLO_f = dh%P_LO   * PNLO_f
       PNLOLO_f = dh%P_NLO  * PLO_f
       PLO3_f   = dh%P_LO   * (dh%P_LO * PLO_f)

       ! first calculate the pure N3LO term, (C_N3LO x f) 
       f2 = (ch%C2N3LO * f)
       fL = (ch%CLN3LO * f)
       f3 = (ch%C3N3LO * f)
       f2_fl11 = (ch%C2N3LO_fl11 * f)
       fL_fl11 = (ch%CLN3LO_fl11 * f)
       ! tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)

       ! Now calculate the (C_LO x P_LO^3 x f) term
       f2 =  (ch%C2LO * PLO3_f)
       fL =  (ch%CLLO * PLO3_f)
       f3 =  (ch%C3LO * PLO3_f)
       tables(9)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO x P_NLO x f) term
       f2 =  (ch%C2LO * PLONLO_f)
       fL =  (ch%CLLO * PLONLO_f)
       f3 =  (ch%C3LO * PLONLO_f)
       tables(10)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO x P_LO x f) term
       f2 =  (ch%C2LO * PNLOLO_f)
       fL =  (ch%CLLO * PNLOLO_f)
       f3 =  (ch%C3LO * PNLOLO_f)
       tables(11)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO^2 x f) term
       f2 =  (ch%C2NLO * PLO2_f)
       fL =  (ch%CLNLO * PLO2_f)
       f3 =  (ch%C3NLO * PLO2_f)
       tables(12)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_NLO) term
       f2 = (ch%C2NLO * PNLO_f)
       fL = (ch%CLNLO * PNLO_f)
       f3 = (ch%C3NLO * PNLO_f)
       tables(13)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NNLO x P_LO) term
       f2 = (ch%C2NNLO * PLO_f)
       fL = (ch%CLNNLO * PLO_f)
       f3 = (ch%C3NNLO * PLO_f)
       tables(14)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NNLO) term
       f2 = (ch%C2LO * PNNLO_f)
       fL = (ch%CLLO * PNNLO_f)
       f3 = (ch%C3LO * PNNLO_f)
       tables(15)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

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
    res(:,F2Z) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(C2_f_NC) + ubarlike(C2_f_NC))*vi2_ai2_Z_up

    ! temporarily put FL into F1;
    res(:,F1Z) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(CL_f_NC) + ubarlike(CL_f_NC))*vi2_ai2_Z_up
    ! then convert to F1
    res(:,F1Z) = (res(:,F2Z) - res(:,F1Z)) / two_xvals

    res(:,F3Z) = (dlike(C3_f) - dbarlike(C3_f))*two_vi_ai_Z_down + &
         &       (ulike(C3_f) - ubarlike(C3_f))*two_vi_ai_Z_up
    res(:,F3Z) = res(:,F3Z)/xValues(grid)

    !--- deal with EM case -----------------------------------------
    res(:,F2EM) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e2_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e2_up
    
    ! temporarily put FL into F1;
    res(:,F1EM) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e2_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e2_up
    ! then convert to F1
    res(:,F1EM) = (res(:,F2EM) - res(:,F1EM)) / two_xvals
    
    !--- deal with gamma-Z case -----------------------------------------
    ! equations taken from section 19 of PDG (19.18)
    res(:,F2gZ) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e_down*two_vi_Z_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e_up * two_vi_Z_up
    
    ! temporarily put FL into F1;
    res(:,F1gZ) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e_down*two_vi_Z_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e_up * two_vi_Z_up
    ! then convert to F1
    res(:,F1gZ) = (res(:,F2gZ) - res(:,F1gZ)) / two_xvals
    
    res(:,F3gZ) = (dlike(C3_f) - dbarlike(C3_f))*e_down*two_ai_Z_down + &
         &        (ulike(C3_f) - ubarlike(C3_f))*e_up * two_ai_Z_up
    res(:,F3gZ) = res(:,F3gZ)/xValues(grid)

    ! for the W cases, it only makes sense to sum over an even number
    ! of light flavours; so save the actual number of flavours, switch
    ! the module-local nf_lcl variable to the nearest (lower) event number
    ! for our W calculations; switch back later
    nf_save = nf_lcl
    nf_lcl = (nf_lcl/2) * 2
    
    !--- deal with Wp case -----------------------------------------
    res(:,F2Wp) = (ulike(C2_f) + dbarlike(C2_f))*vi2_ai2_avg_W
    ! temporarily put FL into F1;
    res(:,F1Wp) = (ulike(CL_f) + dbarlike(CL_f))*vi2_ai2_avg_W
    ! then convert to F1
    res(:,F1Wp) = (res(:,F2Wp) - res(:,F1Wp)) / two_xvals
    res(:,F3Wp) = (ulike(C3_f) - dbarlike(C3_f))*two_vi_ai_avg_W
    res(:,F3Wp) = res(:,F3Wp)/xValues(grid)

    !--- deal with Wm case -----------------------------------------
    res(:,F2Wm) = (dlike(C2_f) + ubarlike(C2_f))*vi2_ai2_avg_W
    ! temporarily put FL into F1;
    res(:,F1Wm) = (dlike(CL_f) + ubarlike(CL_f))*vi2_ai2_avg_W
    ! then convert to F1
    res(:,F1Wm) = (res(:,F2Wm) - res(:,F1Wm)) / two_xvals
    res(:,F3Wm) = (dlike(C3_f) - ubarlike(C3_f))*two_vi_ai_avg_W
    res(:,F3Wm) = res(:,F3Wm)/xValues(grid)

    ! reset nf_lcl to the full (possibly odd) saved value
    nf_lcl = nf_save

    ! AK 16-08-2023: The factor of two is included in SF in the VBF
    ! paper by Maltoni et al. It was kept here until today since the
    ! structure functions were only used with VBF. Now it has been
    ! removed. WARNING: If running with a version of proVBFH pre this
    ! fix then a factor of 2 will be missing...

    ! overall factor of two that we still haven't fully looked into as
    ! of 2015-02-24 [GPS TMP] res = res * two
    !! GPS+AK TMP: we have just included a factor of 2 but this plainly
    !! should not be present for the electromagnetic part; so here
    !! we eliminate it again...
    !res(:,F2EM) = half * res(:,F2EM)
    !res(:,F1EM) = half * res(:,F1EM)

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
    real(dp) :: muR_lcl
    muR_lcl = max(muR,Qmin)
    ! we use alphas from the LHAPDF PDF
    ! alphasLocal = alphasPDF(muR_lcl)
    ! we use alphas from HOPPET
    alphasLocal = Value(coupling, muR_lcl)
  end function alphasLocal
  

  !----------------------------------------------------------------------
  ! F
  ! calculate the structure function at x, muF
  ! this is the sum over all orders
  function StrFct (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7)
    
    if (use_sep_orders) then
       ! if we kept each order separate, then add up all the fixed order terms one by one
       if (order_setup.ge.1) res = F_LO(y, Q, muR, muF)
       if (order_setup.ge.2) res = res + F_NLO(y, Q, muR, muF)
       if (order_setup.ge.3) res = res + F_NNLO(y, Q, muR, muF)
       if (order_setup.ge.4) res = res + F_N3LO(y, Q, muR, muF)
    else
       ! if we haven't kept each order separate, then everything is in tables(1)
       call EvalPdfTable_yQ(tables(1), y, Q, res)
    endif
    
  end function StrFct


  !----------------------------------------------------------------------
  ! F_LO
  ! calculate the leading order structure function at x, muF
  !
  function F_LO (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7)
    real(dp) :: Q_or_muF
    
    ! if scale_choice = 0,1, then evaluate tables at Q
    ! (since it is saved as an array in Q)
    Q_or_muF = Q
    ! if scale_choice >= 2, then evaluate tables at muF
    ! (since it is saved as an array in muF)
    if (scale_choice.ge.2) Q_or_muF = muF

    call EvalPdfTable_yQ(tables(1), y, Q_or_muF, res)
    
  end function F_LO

  !----------------------------------------------------------------------
  ! F_NLO
  ! calculate the next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  !
  function F_NLO (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7)
    real(dp) :: Q_or_muF

    as2pi = alphasLocal(muR) / (twopi)

    ! if scale_choice = 0,1, then evaluate tables at Q
    ! (since it is saved as an array in Q)
    Q_or_muF = Q
    ! if scale_choice >= 2, then evaluate tables at muF
    ! (since it is saved as an array in muF)
    if (scale_choice.ge.2) Q_or_muF = muF
    
    ! C_NLO x f (x) in C1f(:) 
    call EvalPdfTable_yQ(tables(2), y, Q_or_muF, C1f)
    res = C1f
    
    ! if scale_choice = 0,1 then this term is already taken care of
    if (scale_choice.ge.2) then
       LFQ2 = two*log(muF/Q)
       ! C_LO x P_LO x f (x) in C0P0f(:)
       call EvalPdfTable_yQ(tables(3), y, Q_or_muF, C0P0f)
       res = res - C0P0f * LFQ2
    endif
    
    res = res * as2pi
    
  end function F_NLO


  !----------------------------------------------------------------------
  ! F_NNLO
  ! calculate the next-to-next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  function F_NNLO (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7), C2f(-6:7), C0P0sqf(-6:7)
    real(dp) :: C0P1f(-6:7), C1P0f(-6:7)
    real(dp) :: Q_or_muF
    
    as2pi = alphasLocal(muR) / (twopi)

    ! if scale_choice = 0,1, then evaluate tables at Q
    ! (since it is saved as an array in Q)
    Q_or_muF = Q
    ! if scale_choice >= 2, then evaluate tables at muF
    ! (since it is saved as an array in muF)
    if (scale_choice.ge.2) Q_or_muF = muF
    
    ! C_NNLO x f (x) in C2f(:,3) 
    call EvalPdfTable_yQ(tables(4), y, Q_or_muF, C2f)
    res = C2f

    ! if scale_choice = 0,1 then these terms are already taken care of
    if (scale_choice.ge.2) then
       LRQ2 = two*log(muR/Q)
       LFQ2 = two*log(muF/Q)
       ! C_NLO x f (x) in C1f
       call EvalPdfTable_yQ(tables(2), y, Q_or_muF, C1f)
       ! C_LO x P_LO x f (x) in C0P0f
       call EvalPdfTable_yQ(tables(3), y, Q_or_muF, C0P0f)
       ! C_LO x P_LO^2 x f (x) in C0P0sqf
       call EvalPdfTable_yQ(tables(5), y, Q_or_muF, C0P0sqf)
       ! C_LO x P_NLO x f (x) in C0P1f
       call EvalPdfTable_yQ(tables(6), y, Q_or_muF, C0P1f)
       ! C_NLO x P_LO x f (x) in C1P1f
       call EvalPdfTable_yQ(tables(7), y, Q_or_muF, C1P0f)
       ! add up all the different pieces
       res = res - C1P0f * LFQ2 + C0P0sqf * half * LFQ2**2 &
            & - twopi * beta0 * C0P0f * (LRQ2*LFQ2 - half*LFQ2**2) &
            & - C0P1f * LFQ2 + twopi * beta0 * C1f * LRQ2
    endif
    
    res = res * (as2pi)**2
    
  end function F_NNLO

  
  !----------------------------------------------------------------------
  ! F_N3LO
  ! calculate the next-to-next-to-next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  function F_N3LO (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7), C3f(-6:7), C0P0sqf(-6:7), C2f(-6:7)
    real(dp) :: C0P1f(-6:7), C1P0f(-6:7), C0P0cbf(-6:7), C0P01f(-6:7), C0P10f(-6:7)
    real(dp) :: C1P0sqf(-6:7), C1P1f(-6:7), C2P0f(-6:7), C0P2f(-6:7)
    real(dp) :: Q_or_muF
    
    as2pi = alphasLocal(muR) / (twopi)

    ! if scale_choice = 0,1, then evaluate tables at Q
    ! (since it is saved as an array in Q)
    Q_or_muF = Q
    ! if scale_choice >= 2, then evaluate tables at muF
    ! (since it is saved as an array in muF)
    if (scale_choice.ge.2) Q_or_muF = muF
    
    ! C_N3LO x f (x) in C2f(:,8) 
    call EvalPdfTable_yQ(tables(8), y, Q_or_muF, C3f)
    res = C3f

    ! if scale_choice = 0,1 then these terms are already taken care of
    if (scale_choice.ge.2) then
       LRQ2 = two*log(muR/Q)
       LFQ2 = two*log(muF/Q)
       ! C_NLO x f (x) in C1f
       call EvalPdfTable_yQ(tables(2), y, Q_or_muF, C1f)
       ! C_LO x P_LO x f (x) in C0P0f
       call EvalPdfTable_yQ(tables(3), y, Q_or_muF, C0P0f)
       ! C_NNLO x f (x) in C2f(:,3) 
       call EvalPdfTable_yQ(tables(4), y, Q_or_muF, C2f)
       ! C_LO x P_LO^2 x f (x) in C0P0sqf
       call EvalPdfTable_yQ(tables(5), y, Q_or_muF, C0P0sqf)
       ! C_LO x P_NLO x f (x) in C0P1f
       call EvalPdfTable_yQ(tables(6), y, Q_or_muF, C0P1f)
       ! C_NLO x P_LO x f (x) in C1P1f
       call EvalPdfTable_yQ(tables(7), y, Q_or_muF, C1P0f)
       
       ! C_LO x P_LO^3 x f (x) in C0P0cbf
       call EvalPdfTable_yQ(tables(9), y, Q_or_muF, C0P0cbf)
       ! C_LO x P_LO x P_NLO x f (x) in C0P01f
       call EvalPdfTable_yQ(tables(10), y, Q_or_muF, C0P01f)
       ! C_LO x P_LO x P_NLO x f (x) in C0P10f
       call EvalPdfTable_yQ(tables(11), y, Q_or_muF, C0P10f)
       ! C_NLO x P_LO^2 x f (x) in C1P0sqf
       call EvalPdfTable_yQ(tables(12), y, Q_or_muF, C1P0sqf)
       ! C_NLO x P_NLO x f (x) in C1P1f
       call EvalPdfTable_yQ(tables(13), y, Q_or_muF, C1P1f)
       ! C_NNLO x P_LO x f (x) in C2P0f
       call EvalPdfTable_yQ(tables(14), y, Q_or_muF, C2P0f)
       ! C_LO x P_NNLO x f (x) in C0P2f
       call EvalPdfTable_yQ(tables(15), y, Q_or_muF, C0P2f)
    
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

    log_muF2_over_Q2 = two * log(muF(Q) / Q)
    log_muR2_over_Q2 = two * log(muR(Q) / Q)
  end subroutine set_scale_logs
  
  !----------------------------------------------------------------------
  ! mu_R 
  real(dp) function muR(Q)
    real(dp), intent(in) :: Q
    muR = zero
    if (scale_choice.eq.0) then
       ! scale_choice = 0 : muF = xmuF * cst_mu
       muR = xmuR * cst_mu
    elseif (scale_choice.eq.1) then
       ! scale_choice = 1 : mu_R = xmuR * Q
       muR = xmuR * Q
    else
       call wae_error('muR(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muR

  !----------------------------------------------------------------------
  ! mu_F 
  real(dp) function muF(Q)
    real(dp), intent(in) :: Q
    muF = zero
    if (scale_choice.eq.0) then
       ! scale_choice = 0 : muF = xmuF * cst_mu
       muF = xmuF * cst_mu
    elseif (scale_choice.eq.1) then
       ! scale_choice = 1 : muF = xmuF * Q
       muF = xmuF * Q
    else
       call wae_error('muF(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muF
  
  subroutine use_vfns (f, Q)
    implicit none
    real(dp), intent(in) :: Q
    real(dp), intent(inout) :: f(0:grid%ny,ncompmin:ncompmax)
    integer :: inf, current_nf, dummy!, nfhi

    if(ch%nflo.lt.3) stop 'Do you really want to run with fewer than 3 light flavours??'

    ! First set current nf equal to the smallest value it can take
    current_nf = ch%nflo
    if(Q.lt.quark_masses_sf(4)) then   
      current_nf = 3
    elseif(Q.lt.quark_masses_sf(5)) then 
      current_nf = 4
    elseif(Q.lt.quark_masses_sf(6)) then  
      current_nf = 5
    elseif(Q.ge.quark_masses_sf(6)) then  
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

end module structure_functions
