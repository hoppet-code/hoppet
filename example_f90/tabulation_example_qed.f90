!! An example program using a tabulation. It should yield an output identical
!! to the content of the file tabulation_example_qed.default_output
!!
!!
program tabulation_example_qed
  use hoppet
  
  !! if using LHAPDF, rename the hoppet functions EvolvePDF and InitPDF, which
  !! would otherwise conflict with LHAPDF, by replacing the use
  !! statement above with:
  !use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF

  implicit none
  real(dp) :: dy, ymax
  ! the interpolation order (not the perturbative order!)
  integer  :: order
  ! the number of loops in QCD and QED
  integer  :: nloop_qcd, nqcdloop_qed
  ! whether to include the alpha^2 Plq term
  logical  :: with_Plq_nnloqed
  !! holds information about the grid in y=log(1/x)
  type(grid_def) :: grid
  !! holds the splitting functions
  type(dglap_holder) :: dh
  !! holds the PDF tabulation
  type(pdf_table)       :: table
  !! holds the coupling
  real(dp)               :: effective_light_quark_masses, quark_masses(4:6)
  type(running_coupling) :: coupling
  !! holds the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp) :: Q0
  !! holds results at some x, Q
  real(dp) :: Q, pdf_at_xQ(ncompmin:ncompmaxLeptons)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix
  !! Additional things for QED
  type(qed_coupling)  :: coupling_qed
  type(qed_split_mat) :: qed_split


  !! define the interfaces for LHA pdf (by default not used)
  !! (NB: unfortunately this conflicts with an internal hoppet name,
  !! so make sure that you "redefine" the internal hoppet name, 
  !! as illustrated in the commented "use" line above:
  !! use hoppet, EvolvePDF_hoppet => EvolvePDF, ...)
  ! interface
  !    subroutine EvolvePDF(x,Q,res)
  !      use types; implicit none
  !      real(dp), intent(in)  :: x,Q
  !      real(dp), intent(out) :: res(*)
  !    end subroutine EvolvePDF
  ! end interface

  ! set up parameters for grid
  order = -6
  ymax  = 12.0_dp
  dy    = 0.1_dp

  ! set up the grid itself (this call sets up a nested grid composed of 4 subgrids)
  call InitGridDefDefault(grid, dy, ymax, order=order)

  ! set variables to specify the perturbative accuracy
  nloop_qcd = 3
  nqcdloop_qed = 1
  with_Plq_nnloqed = .false.
  
  ! initialise the splitting-function holder
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop_qcd,nflo=3,nfhi=6)

  ! and the QED splitting matrices
  call InitQEDSplitMat(grid, qed_split)

  write(6,'(a)') "Splitting functions initialised!"

  ! allocate a pdf0 array including leptons
  call AllocPDFWithLeptons(grid, pdf0)

  ! initialise a PDF from the function below (must be contained,
  ! in a "used" module, or with an explicitly defined interface)
  ! pdf0(:,-6:6) are the standard pdf for quarks and gluons,
  ! pdf0(:,8) is the photon density,
  ! pdf0(:,9:11) are the sum of lepton-antilepton densities
  ! for the electron, the muon and the tau.
  pdf0(:,:) = unpolarized_dummy_pdf(xValues(grid))  

  Q0 = sqrt(2.0_dp)  ! the initial scale

  ! allocate and initialise the running coupling with a given
  ! set of quark masses (NB: charm mass just above Q0).
  quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  call InitRunningCoupling(coupling,alfas=0.35_dp,Q=Q0,nloop=nloop_qcd,&
       &                   quark_masses = quark_masses)

  ! Initialises coupling_qed object, including thresholds. All
  ! light quarks turn on simultaneously at some specified effective scale;
  ! the choice of 0.109 GeV results in a fairly accurate value of alpha(mZ)
  ! and alpha(m_tau).
  effective_light_quark_masses = 0.109_dp
  call InitQEDCoupling(coupling_qed, effective_light_quark_masses, quark_masses(4:6))

  ! create the table object that will contain the values of the parton
  ! densities in a grid of x and Q interpolation points.
  call AllocPdfTableWithLeptons(grid, table, Qmin=1.0_dp, Qmax=28000.0_dp, & 
       & dlnlnQ = dy/4.0_dp, freeze_at_Qmin=.true.)
  ! add information about the nf transitions to the table (improves
  ! interpolation quality)
  call AddNfInfoToPdfTable(table,coupling)

  ! create the tabulation based on the evolution of pdf0 from scale Q0
  ! up to qmax and down to qmin
  call EvolvePdfTableQED(table, Q0, pdf0, dh, qed_split, &
                         coupling, coupling_qed, nloop_qcd, nqcdloop_qed, with_Plq_nnloqed)

  write(6,'(a)') "Evolution done!"

  ! get the value of the tabulation at some point
  Q = 91.1876_dp;
  write(6,'(a)')
  write(6,'(a,f9.4,a,f8.5)') " At Q = ",Q," GeV, QCD coupling is alpha_s     = ", value(coupling, Q)
  write(6,'(a,f9.4,a,f8.3)') " At Q = ",Q," GeV, QED coupling is 1/alpha_qed = ", one/value(coupling_qed, Q)
  write(6,'(a,f6.4,a,f8.3)') " At mtau = ",m_tau," GeV, QED coupling is 1/alpha_qed = ", one/value(coupling_qed, m_tau/1.0000001_dp)
  write(6,'(a)')

  ! Get values at some points and print them 
  Q = 100.0_dp
  write(6,'(a,f9.4,a)') "           Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,2a12,a14,a10,a11,4a13)') "x",&
       & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon",&
       & "photon","e+ + e-","mu+ + mu-","tau+ + tau-"
  do ix = 1, size(heralhc_xvals)
     call EvalPdfTable_xQ(table,heralhc_xvals(ix),Q,pdf_at_xQ)
     write(6,'(es7.1,9es12.4)') heralhc_xvals(ix), &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  pdf_at_xQ(0), &
          &  pdf_at_xQ(iflv_photon), &
          &  pdf_at_xQ(iflv_electron), &
          &  pdf_at_xQ(iflv_muon),&
          &  pdf_at_xQ(iflv_tau)          
  end do
  
  ! some cleaning up (not strictly speaking needed, but illustrates
  ! how it's done)
  call Delete(table)
  call Delete(pdf0)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)
  call Delete(qed_split)
  call Delete(coupling_qed)

contains 
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmaxLeptons)
    real(dp) :: gluon(size(xvals)), uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    ! (not actually needed after setting pdf=0
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    gluon = N_g * xvals**(-0.1_dp) * (1-xvals)**5
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)

    ! labels iflv_g, etc., come from the hoppet module, inherited
    ! from the main program
    pdf(:, iflv_g) = 0.99_dp * gluon
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + 0.9_dp * ubar
    pdf(:,-iflv_u) = 0.9_dp * ubar
    pdf(:, iflv_d) = dv + 0.9_dp * dbar
    pdf(:,-iflv_d) = 0.9_dp * dbar
    pdf(:,iflv_photon) = gluon/100._dp + ubar/10._dp 
    pdf(:,iflv_electron) = dbar/10._dp
    pdf(:,iflv_muon) = ubar/10._dp
    pdf(:,iflv_tau) = dbar/10._dp         

  end function unpolarized_dummy_pdf

end program tabulation_example_qed


