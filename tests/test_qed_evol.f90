program test_qed_evol
  use hoppet_v1
  use qed_evolution; use qed_objects; use qed_coupling_module
  implicit none
  type(grid_def)      :: grid
  type(qed_split_mat) :: qed_split
  type(dglap_holder)  :: dh
  real(dp)               :: quark_masses(4:6)
  type(running_coupling) :: coupling
  type(qed_coupling)     :: coupling_qed
  real(dp)               :: ymax, dy
  real(dp), pointer :: pdf_in(:,:), pdf_out(:,:)
  real(dp)         :: Q0, Qfinal
  integer :: nloop_qcd, nqcdloop_qed
  
  ymax  = 12.0_dp
  dy    = 0.1_dp
  nloop_qcd = 3
  nqcdloop_qed = 0
  !Q0 = sqrt(2.0_dp)  ! the initial scale
  Q0 = 1.7_dp  ! the initial scale
  Qfinal = 40.0_dp
  
  call InitGridDefDefault(grid, dy, ymax)
  call SetDefaultEvolutionDu(dy/3.0_dp)  ! generally a good choice
  call InitQEDSplitMat(grid, qed_split)
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop_qcd,nflo=3,nfhi=6)

  quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  call InitRunningCoupling(coupling,alfas=0.35_dp,Q=sqrt(two),nloop=nloop_qcd,&
       &                   quark_masses = quark_masses)
  call InitQEDCoupling(coupling_qed, quark_masses(4), quark_masses(5), quark_masses(6))

  call AllocPDFWithLeptons(grid, pdf_in)
  call AllocPDFWithLeptons(grid, pdf_out)
  pdf_in = zero
  ! fill in the QCD part
  pdf_in(:,:ncompmax) = unpolarized_dummy_pdf(xValues(grid))
  ! then add in a photon of some kind
  pdf_in(:,iflv_photon) = Value(coupling_qed, Q0) * (one - xValues(grid))**4

  call QEDQCDEvolvePDF(dh, qed_split, pdf_in, coupling, coupling_qed,&
       &               Q0, Qfinal, nloop_qcd, nqcdloop_qed)
  !! NEED TO
  !! - TEST AGAINST PLAIN QCD EVOLUTION (with QED alpha turned low)
  !! - test momentum sum rule, with both nqcdloop_qed values
  !! - compare to photon evolution from another source?
contains

  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
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
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)

    ! labels iflv_g, etc., come from the hoppet_v1 module, inherited
    ! from the main program
    pdf(:, iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf


end program test_qed_evol
