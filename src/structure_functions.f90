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
  use qcd_coupling
  use dglap_holders; use dglap_choices
  use pdf_representation
  use pdf_tabulate
  use dummy_pdfs
  use types; use consts_dp; use splitting_functions
  use convolution
  use qcd; use warnings_and_errors
  use coefficient_functions_holder
  implicit none

  private
  public :: InitStrFct
  public :: read_PDF
  public :: set_structure_functions_up_to
  public :: write_f1
  public :: write_f2
  public :: write_f3
  
  real(dp), external :: alphasPDF
  type(running_coupling), save :: coupling
  !! holds information about the grid
  type(grid_def),     save :: grid, gdarray(4)
  type(grid_def),     save :: grid_n3lo, gdarray_n3lo(4)

  !! holds the splitting functions
  type(dglap_holder), save :: dh

  !! 0 is main pdf table, while i=1:8 contain convolutions with the
  !! splitting function
  type(pdf_table), save :: tables(0:15)
  ! indices for the different structure functions
  integer, parameter :: F1Wp= 1, F2Wp= 2, F3Wp= 3
  integer, parameter :: F1Wm=-1, F2Wm=-2, F3Wm=-3
  integer, parameter :: F1Z = 4, F2Z = 5, F3Z = 6
  integer, parameter :: F1EM = -4, F2EM = -5
  
  ! constants and fixed parameters
  real(dp), parameter :: viW = 1/sqrt(two), aiW = viW
  real(dp), save      :: vi2_ai2_avg_W, two_vi_ai_avg_W
  real(dp), save      :: vi2_ai2_Z_down, vi2_ai2_Z_up
  real(dp), save      :: two_vi_ai_Z_down, two_vi_ai_Z_up
  real(dp), parameter :: e2_up = 4.0_dp/9.0_dp, e2_down = 1.0_dp/9.0_dp
  ! these log terms are only used for scale_choice = 0,1, to speed up the code
  real(dp), save      :: log_muF2_over_Q2, log_muR2_over_Q2
  real(dp), save      :: Qmin, toy_Q0, dglap_Q0
  integer, save       :: nflav
  real(dp), public, save :: toy_alphas_Q0 = 0.1185_dp
  type(running_coupling), public, save :: toy_coupling

  ! scale choices
  real(dp), save :: cst_muR, cst_muF, xmuR, xmuF, xmuR_PDF
  integer, save  :: scale_choice

  real(dp), save        :: C2LO,   CLLO,   C3LO
  type(split_mat), save :: C2NLO,  CLNLO,  C3NLO
  type(split_mat), save :: C2NNLO, CLNNLO, C3NNLO
  type(split_mat), save :: C2N3LO, CLN3LO, C3N3LO
  type(split_mat), save :: C2N3LO_fl11, CLN3LO_fl11
  integer :: nf_lcl


contains

  !----------------------------------------------------------------------
  ! This routine assumes that LHAPDF has already been initialised with a PDF
  subroutine InitStrFct(rts, order_max, nf, xR, xF, sc_choice, muR, muF, toyQ0, dglapQ0, xR_PDF)
    real(dp), intent(in) :: rts, xR, xF
    integer, intent(in)  :: order_max, nf
    integer, optional    :: sc_choice
    real(dp), optional   :: muR, muF, toyQ0, dglapQ0, xR_PDF
    !----------------------------------------------------------------------
    real(dp) :: ymax, dy, minQval, maxQval, dlnlnQ
    real(dp) :: sin_thw
    integer  :: nloop, order
    if(present(sc_choice))then
       scale_choice=sc_choice
    else
       scale_choice=1
    endif
    if(present(muR)) cst_muR=muR
    if(present(muF)) cst_muF=muF
    if(present(toyQ0)) toy_Q0=toyQ0
    if(present(dglapQ0)) dglap_Q0=dglapQ0
    if(present(xR_PDF)) xmuR_PDF=xR_PDF
    xmuR = xR
    xmuF = xF
    
    ! where should Qmin and sin_thw go ?
    ! computed from W,Z mass
    sin_thw = 0.22289722252391826839_dp
    Qmin = 1.0_dp ! TEMPORARY 
    ! call getQ2min(0,Qmin)
    ! Qmin = sqrt(Qmin)
    
    ! evaluate parameters needed for the structure functions
    ! cf. Eqs. (3.10+11) of 1109.3717
    vi2_ai2_Z_up     = one/four + (half - (four/three) * sin_thw)**2
    vi2_ai2_Z_down   = one/four + (half - (two/three)  * sin_thw)**2
    two_vi_ai_Z_up   = half - (four/three) * sin_thw
    two_vi_ai_Z_down = half - (two/three)  * sin_thw

    ! cf. Eq.3.20 + 3.17 & 3.18
    ! 1/nf \sum_j=1^nf (vj^2 + aj^2)
    vi2_ai2_avg_W = viW**2 + aiW**2
    ! 1/nf \sum_j=1^nf 2*vj*aj
    two_vi_ai_avg_W = two * viW * aiW

    ! Streamlined initialization
    ! including  parameters for x-grid
    order = -6 
    ymax  = 16.0_dp
    dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0_dp
    nloop = 3
    minQval = min(xmuF*Qmin, Qmin)
    maxQval = max(xmuF*rts, rts)
    
    call hoppetStartExtendedLocal(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)

    call qcd_SetNf(nflav)
    nf_lcl = nf_int
    write(6,*) "nf_lcl = ", nf_lcl

    ! first initialise the LO coefficient "functions" (just a number, since delta-fn)
    C2LO = one
    CLLO = zero
    C3LO = one
    
    ! now initialise some NLO coefficient functions
    call InitC2NLO(grid, C2NLO)
    call InitCLNLO(grid, CLNLO)
    call InitC3NLO(grid, C3NLO) 

    ! and the NNLO ones
    call InitC2NNLO(grid, C2NNLO)
    call InitCLNNLO(grid, CLNNLO)
    call InitC3NNLO(grid, C3NNLO) 

    ! and the N3LO ones
    call InitC2N3LO(grid_n3lo, C2N3LO)
    call InitCLN3LO(grid_n3lo, CLN3LO)
    call InitC3N3LO(grid_n3lo, C3N3LO)
    ! including the fl11 terms for Z/photon exchanges
    call InitC2N3LO_fl11(grid_n3lo, C2N3LO_fl11)
    call InitCLN3LO_fl11(grid_n3lo, CLN3LO_fl11)

    ! read the PDF in
    call read_PDF()
    ! set up all the strucutre functions
    call set_structure_functions_up_to(order_max)
    
  end subroutine InitStrFct

  !----------------------------------------------------------------------
  subroutine hoppetStartExtendedLocal(ymax,dy,valQmin,valQmax,dlnlnQ,nloop,order,factscheme)
    implicit none
    real(dp), intent(in) :: ymax   !! highest value of ln1/x user wants to access
    real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
    real(dp), intent(in) :: valQmin, valQmax !! range in Q
    real(dp), intent(in) :: dlnlnQ !! internal table spacing in lnlnQ
    integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=3)
    integer,  intent(in) :: order  !! order of numerical interpolation (+ve v. -ve: see below)
    integer,  intent(in) :: factscheme !! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar
    !-------------------------------------

    ! initialise our grids

    ! the internal interpolation order (with a minus sign allows
    ! interpolation to take fake zero points beyond x=1 -- convolution
    ! times are unchanged, initialisation time is much reduced and
    ! accuracy is slightly reduced)
    !order = -5
    ! Now create a nested grid
    call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
    call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp,  order=order)
    call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp,  order=order)
    call InitGridDef(gdarray(1),dy,       ymax  ,  order=order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)

    ! At N3LO we need to change the grid slightly:
    !  - convolution epsilon is reduced from 10^-7 to 10^-6
    !  - interpolation order is reduced from -6 to -5 for the high x region
    ! Note: The grid points in y have to be the same as for the other grid!
    call SetDefaultConvolutionEps(0.000001_dp) ! anything less breaks integration of N3LO pieces
    call InitGridDef(gdarray_n3lo(4),dy/27.0_dp,0.2_dp, order=-abs(order)+1)!reduce order to -5
    call InitGridDef(gdarray_n3lo(3),dy/9.0_dp,0.5_dp,  order=-abs(order)+1) !reduce order to -5
    call InitGridDef(gdarray_n3lo(2),dy/3.0_dp,2.0_dp,  order=order)
    call InitGridDef(gdarray_n3lo(1),dy,       ymax  ,  order=order)
    call InitGridDef(grid_n3lo,gdarray_n3lo(1:4),locked=.true.)

    ! create the tables that will contain our copy of the user's pdf
    ! as well as the convolutions with the pdf.
    call AllocPdfTable(grid, tables(:), valQmin, valQmax, &
         & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)
    
    ! tables(0) : pdf
    tables(0)%tab = zero
    ! tables(1) : LO structure function : C_LO x f
    tables(1)%tab = zero
    
    ! tables(2) : NLO structure function : C_NLO x f
    tables(2)%tab = zero
    ! tables(3) : NLO contribution : C_LO x P_LO x f
    tables(3)%tab = zero
    
    ! tables(4) : NNLO structure function : C_NNLO x f
    tables(4)%tab = zero
    ! tables(5) : NNLO contribution : C_LO x P_LO^2 x f
    tables(5)%tab = zero
    ! tables(6) : NNLO contribution : C_LO x P_NLO x f
    tables(6)%tab = zero
    ! tables(7) : NNLO contribution : C_NLO x P_LO x f
    tables(7)%tab = zero
    
    ! tables(8)  : N3LO contribution : C_N3LO x f
    tables(8)%tab = zero
    ! tables(9)  : N3LO contribution : C_LO x P_LO^3 x f
    tables(9)%tab = zero
    ! tables(10) : N3LO contribution : C_LO x P_LO x P_NLO x f
    tables(10)%tab = zero
    ! tables(11) : N3LO contribution : C_LO x P_NLO x P_LO x f
    tables(11)%tab = zero
    ! tables(12) : N3LO contribution : C_NLO x P_LO^2 x f
    tables(12)%tab = zero
    ! tables(13) : N3LO contribution : C_NLO x P_NLO x f
    tables(13)%tab = zero
    ! tables(14) : N3LO contribution : C_NNLO x P_LO x f
    tables(14)%tab = zero
    ! tables(15) : N3LO contribution : C_LO x P_NNLO x f
    tables(15)%tab = zero
    
    ! initialise splitting-function holder
    call InitDglapHolder(grid,dh,factscheme=factscheme,&
         &                      nloop=nloop,nflo=3,nfhi=6)
    ! choose a sensible default number of flavours.
    call SetNfDglapHolder(dh,nflcl=nflav)
  end subroutine hoppetStartExtendedLocal


  !----------------------------------------------------------------------
  ! fill the streamlined interface PDF table (possibly using hoppet's
  ! evolution)
  subroutine read_PDF()
    !real(dp) :: muR_Q
    !real(dp) :: Q0pdf, asMZ
    !! define the interfaces for LHA pdf (by default not used)
    !! (NB: unfortunately this conflicts with an internal hoppet name,
    !! so make sure that you "redefine" the internal hoppet name,
    !! as illustrated in the commented "use" line above:
    !! use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, ...)
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    !logical, save :: pre_ev_done = .false.
    !----------------
    real(dp) :: res_lhapdf(-6:6), x, Q
    real(dp) :: res_hoppet(-6:6)
    real(dp) :: toy_pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    integer  :: ix
    real(dp) :: xvals(0:grid%ny)
    real(dp), parameter :: mz = 91.2_dp
    
    ! ! Set parameters of running coupling
    ! asMZ = alphasPDF(MZ)
    ! Q0pdf = 10.0_dp
    ! muR_Q = 1.0_dp
    ! if (hoppet_evolution) then
    !    if (hoppet_pre_evolution) then
    !       if (.not. pre_ev_done) then
    !          call hoppetPreEvolve(asMZ, MZ, nloop, muR_Q, Q0pdf)
    !          pre_ev_done = .true.
    !       end if
    !       call hoppetCachedEvolve(EvolvePDF)
    !    else
    !       call hoppetEvolve(asMZ, MZ, nloop, muR_Q, EvolvePDF, Q0pdf)
    !    end if
    !    write(6,'(a)') "Evolution done!"
    ! else
    !    call hoppetAssign(EvolvePDF)
    !    write(6,'(a)') "PDF assigned to hoppet table!"
    ! end if
    ! !write(0,*) "Size of table = ", size(tables(0)%tab)
    
    if (toy_Q0 > zero) then
       write(6,*) "WARNING: Using toy PDF"
       toy_pdf_at_Q0 = unpolarized_dummy_pdf(xValues(grid))
       call InitRunningCoupling(toy_coupling, alfas=toy_alphas_Q0, &
            &                   nloop = 3, Q = toy_Q0, fixnf=nf_int)
       call EvolvePdfTable(tables(0), toy_Q0, toy_pdf_at_Q0, dh, toy_coupling, nloop=3)
    elseif (dglap_Q0 > zero) then

       write(6,*) "WARNING: Using internal HOPPET DGLAP evolution"
       xvals = xValues(grid)

       do ix = 0, grid%ny
          call EvolvePDF(xvals(ix),dglap_Q0,pdf_at_Q0(ix,:))
       enddo
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, (/  1.275_dp, 4.18_dp, 173.21_dp/),&
            & .true.)
       call EvolvePdfTable(tables(0), dglap_Q0, pdf_at_Q0, dh, coupling, &
            &  muR_Q=xmuR_PDF, nloop=3)
    else
    ! InitRunningCoupling has to be called for the HOPPET coupling to be initialised 
    ! Default is to ask for 4 loop running and threshold corrections at quark masses.  
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, (/  1.275_dp, 4.18_dp, 173.21_dp/),&
            & .true., xmuR_PDF)
       ! InitRunningCoupling(coupling, alfas , Q , nloop,&
       !     & fixnf, quarkmasses,& ! fixnf can be set to a positive number for
                                    ! fixed nf. -1000000045 gives variable nf
                                    ! and threshold corrections at quarkmasses.
       !     & masses_areMSbar)
       call hoppetAssign(EvolvePDF)
    end if
    
    ! quickly test that we have read in the PDFs correctly
    write(6,*) "Quick test that PDFs have been read in correctly"
    x = 0.08_dp
    Q = 17.0_dp
    call EvolvePDF(x, Q, res_lhapdf)
    call EvalPdfTable_xQ(tables(0), x, Q, res_hoppet)
    write(6,*) 'lhapdf: ', res_lhapdf
    write(6,*) 'hoppet: ', res_hoppet
  end subroutine read_PDF

  !----------------------------------------------------------------------
  !! Given a pdf_subroutine with the interface shown below, initialise
  !! our internal pdf table.
  subroutine hoppetAssign(pdf_subroutine)
    implicit none
    interface ! indicate what "interface" pdf_subroutine is expected to have
       subroutine pdf_subroutine(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine pdf_subroutine
    end interface
    !-----------------------------------

    ! set up table(0) by copying the values returned by pdf_subroutine onto 
    ! the x-Q grid in table(0)
    call FillPdfTable_LHAPDF(tables(0), pdf_subroutine)
  end subroutine hoppetAssign

  !----------------------------------------------------------------------
  ! Set up the structure function tables up to a given order
  subroutine set_structure_functions_up_to(order)
    integer, intent(in) :: order

    ! To turn off b quarks completely (only for testing and comparison)
    ! uncomment following two lines:
    ! tables(0)%tab(:,-5,:) = zero
    ! tables(0)%tab(:,+5,:) = zero
    
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

    ! To rescale PDF by the N3LO F2 structure function (evaluated at 10 GeV)
    ! as a check of the size of missing N3LO PDFs, uncomment following line:
    ! call rescale_pdf_nnlo(8.0_dp,order)
    ! call rescale_pdf_n3lo(8.0_dp,order)
  end subroutine set_structure_functions_up_to

  !----------------------------------------------------------------------
  ! Rescale the PDF by F2 NNLO / F2 N3LO
  subroutine rescale_pdf_n3lo (Qval, order)
    real(dp), intent(in) :: Qval
    integer, intent(in)  :: order
    real(dp) :: mF, mR, factor
    real(dp) :: str_fct(-6:6), f2nnlo, f2n3lo, y(0:grid%ny)
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
    real(dp) :: str_fct(-6:6), f2nlo, f2nnlo, y(0:grid%ny) 
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
    real(dp) :: ytest, xval, mR, mF, F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO, res(-6:6)
    integer  :: iy
    !F1 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a,a)') '# x  F1Wp(LO) F1Wm(LO) F1Wp(NLO) F1Wm(NLO) F1Wp(NNLO) F1Wm(NNLO)', &
         & ' F1Wp(N3LO) F1Wm(N3LO) F1Z(LO) F1Z(NLO) F1Z(NNLO) F1Z(N3LO)'
    mF = muF(Qtest)
    mR = muR(Qtest)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
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
       write(idev,*)
    end do
  end subroutine write_f1
  
  !----------------------------------------------------------------------
  ! write the F2 structure function to idev
  subroutine write_f2 (idev, Qtest, ymax, ny)
    real(dp), intent(in) :: Qtest, ymax
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO, res(-6:6)
    integer  :: iy
    !F2 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a,a)') '# x  F2Wp(LO) F2Wm(LO) F2Wp(NLO) F2Wm(NLO) F2Wp(NNLO) F2Wm(NNLO)', &
         & ' F2Wp(N3LO) F2Wm(N3LO) F2Z(LO) F2Z(NLO) F2Z(NNLO) F2Z(N3LO)'
    mF = muF(Qtest)
    mR = muR(Qtest)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
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
       write(idev,*)
    end do
  end subroutine write_f2

  !----------------------------------------------------------------------
  ! write the F3 structure function to idev
  subroutine write_f3 (idev, Qtest, ymax, ny)
    real(dp), intent(in) :: Qtest, ymax
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, mR, mF, F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO, res(-6:6)
    integer  :: iy
    !F3 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a,a)') '# x  F3Wp(LO) F3Wm(LO) F3Wp(NLO) F3Wm(NLO) F3Wp(NNLO) F3Wm(NNLO)', &
         & ' F3Wp(N3LO) F3Wm(N3LO) F3Z(LO) F3Z(NLO) F3Z(NNLO) F3Z(N3LO)'
    mF = muF(Qtest)
    mR = muR(Qtest)
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
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
       write(idev,*)
    end do
  end subroutine write_f3


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
       tables(1)%tab(:,:,iQ) = structure_function_general(C2LO*f, CLLO*f, C3LO*f)
    end do
    
  end subroutine set_LO_structure_functions
  
  !----------------------------------------------------------------------
  ! Set up the LO structure functions for any scale choice
  subroutine set_LO_structure_functions_anyscale()
    integer :: iQ

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       tables(1)%tab(:,:,iQ) = structure_function_general(&
            & tables(0)%tab(:,:,iQ) * C2LO, &
            & tables(0)%tab(:,:,iQ) * CLLO, &
            & tables(0)%tab(:,:,iQ) * C3LO)       
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
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
       PLO_f = dh%P_LO * f
       f2 = CxNLO_with_logs(C2LO, C2NLO, f, PLO_f)
       fL = CxNLO_with_logs(CLLO, CLNLO, f, PLO_f)
       f3 = CxNLO_with_logs(C3LO, C3NLO, f, PLO_f)
       
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
       f2 = (C2NLO * f)
       fL = (CLNLO * f)
       f3 = (C3NLO * f)
       tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now compute the (C_LO x P_LO x f) term
       PLO_f = dh%P_LO * f
       f2 = (C2LO * PLO_f)
       fL = (CLLO * PLO_f)
       f3 = (C3LO * PLO_f)
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
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
       PLO_f   = dh%P_LO  * f
       PNLO_f  = dh%P_NLO * f
       PLO2_f  = dh%P_LO  * PLO_f
       f2 = CxNNLO_with_logs(C2LO, C2NLO, C2NNLO, f, PLO_f, PNLO_f, PLO2_f)
       fL = CxNNLO_with_logs(CLLO, CLNLO, CLNNLO, f, PLO_f, PNLO_f, PLO2_f)
       f3 = CxNNLO_with_logs(C3LO, C3NLO, C3NNLO, f, PLO_f, PNLO_f, PLO2_f)

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
       f2 = (C2NNLO * f)
       fL = (CLNNLO * f)
       f3 = (C3NNLO * f)
       tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO^2 x f) term
       f2 =  (C2LO * PLO2_f)
       fL =  (CLLO * PLO2_f)
       f3 =  (C3LO * PLO2_f)
       tables(5)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO) term
       f2 = (C2LO * PNLO_f)
       fL = (CLLO * PNLO_f)
       f3 = (C3LO * PNLO_f)
       tables(6)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO) term
       f2 = (C2NLO * PLO_f)
       fL = (CLNLO * PLO_f)
       f3 = (C3NLO * PLO_f)
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
    
       f2 = CxN3LO_with_logs(C2LO, C2NLO, C2NNLO, C2N3LO, f, PLO_f, PNLO_f, PLO2_f, &
            &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       fL = CxN3LO_with_logs(CLLO, CLNLO, CLNNLO, CLN3LO, f, PLO_f, PNLO_f, PLO2_f, &
            &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       f3 = CxN3LO_with_logs(C3LO, C3NLO, C3NNLO, C3N3LO, f, PLO_f, PNLO_f, PLO2_f, &
            &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       ! now the fl_11 piece
       f2_fl11 = C2N3LO_fl11 * f
       fL_fl11 = CLN3LO_fl11 * f

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
       f2 = (C2N3LO * f)
       fL = (CLN3LO * f)
       f3 = (C3N3LO * f)
       f2_fl11 = (C2N3LO_fl11 * f)
       fL_fl11 = (CLN3LO_fl11 * f)
       ! tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)

       ! Now calculate the (C_LO x P_LO^3 x f) term
       f2 =  (C2LO * PLO3_f)
       fL =  (CLLO * PLO3_f)
       f3 =  (C3LO * PLO3_f)
       tables(9)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO x P_NLO x f) term
       f2 =  (C2LO * PLONLO_f)
       fL =  (CLLO * PLONLO_f)
       f3 =  (C3LO * PLONLO_f)
       tables(10)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO x P_LO x f) term
       f2 =  (C2LO * PNLOLO_f)
       fL =  (CLLO * PNLOLO_f)
       f3 =  (C3LO * PNLOLO_f)
       tables(11)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO^2 x f) term
       f2 =  (C2NLO * PLO2_f)
       fL =  (CLNLO * PLO2_f)
       f3 =  (C3NLO * PLO2_f)
       tables(12)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_NLO) term
       f2 = (C2NLO * PNLO_f)
       fL = (CLNLO * PNLO_f)
       f3 = (C3NLO * PNLO_f)
       tables(13)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NNLO x P_LO) term
       f2 = (C2NNLO * PLO_f)
       fL = (CLNNLO * PLO_f)
       f3 = (C3NNLO * PLO_f)
       tables(14)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NNLO) term
       f2 = (C2LO * PNNLO_f)
       fL = (CLLO * PNNLO_f)
       f3 = (C3LO * PNNLO_f)
       tables(15)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine set_N3LO_structure_functions_anyscale


 !----------------------------------------------------------------------
 ! Structure function valid up to NNLO
  function structure_function_general(C2_f, CL_f, C3_f) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp) :: C2_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: res(0:ubound(C2_f,dim=1), ncompmin:ncompmax)
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
    real(dp)             :: res(0:ubound(C2_f,dim=1), ncompmin:ncompmax)
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
    
    ! overall factor of two that we still haven't fully looked into as
    ! of 2015-02-24 [GPS TMP]
    res = res * two
    !! GPS+AK TMP: we have just included a factor of 2 but this plainly
    !! should not be present for the electromagnetic part; so here
    !! we eliminate it again...
    res(:,F2EM) = half * res(:,F2EM)
    res(:,F1EM) = half * res(:,F1EM)

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
    if (toy_Q0 > zero) then
       ! we use alphas associated with the toy PDF
       alphasLocal = Value(toy_coupling, muR)
    else
       muR_lcl = max(muR,Qmin)
       ! we use alphas from the LHAPDF PDF
       ! alphasLocal = alphasPDF(muR_lcl)
       ! we use alphas from HOPPET
       alphasLocal = Value(coupling, muR_lcl)
    end if
  end function alphasLocal
  

  !----------------------------------------------------------------------
  ! F_LO
  ! calculate the leading order structure function at x, muF
  !
  function F_LO (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:6)
    real(dp) :: C1f(-6:6), C0P0f(-6:6)
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
    real(dp) :: res(-6:6), as2pi, LFQ2
    real(dp) :: C1f(-6:6), C0P0f(-6:6)
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
    real(dp) :: res(-6:6), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:6), C0P0f(-6:6), C2f(-6:6), C0P0sqf(-6:6)
    real(dp) :: C0P1f(-6:6), C1P0f(-6:6)
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
    real(dp) :: res(-6:6), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:6), C0P0f(-6:6), C3f(-6:6), C0P0sqf(-6:6), C2f(-6:6)
    real(dp) :: C0P1f(-6:6), C1P0f(-6:6), C0P0cbf(-6:6), C0P01f(-6:6), C0P10f(-6:6)
    real(dp) :: C1P0sqf(-6:6), C1P1f(-6:6), C2P0f(-6:6), C0P2f(-6:6)
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
       ! scale_choice = 0 : muF = xmuF * cst_muR
       muR = xmuR * cst_muR
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
       ! scale_choice = 0 : muF = xmuF * cst_muF
       muF = xmuF * cst_muF
    elseif (scale_choice.eq.1) then
       ! scale_choice = 1 : muF = xmuF * Q
       muF = xmuF * Q
    else
       call wae_error('muF(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muF
  
end module structure_functions
