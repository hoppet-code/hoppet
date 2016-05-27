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
  real(dp)         :: Q0, Qfinal, moments(ncompmin:ncompmaxLeptons)
  integer :: nloop_qcd, nqcdloop_qed
  
  ymax  = 25.0_dp
  dy    = 0.1_dp
  nloop_qcd = 3
  nqcdloop_qed = 1
  Q0 = sqrt(2.0_dp)  ! the initial scale
  !Q0 = 1.7_dp  ! the initial scale
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
  write(0,*) "# QED coupling at 1 GeV = ", Value(coupling_qed,one)
  
  call AllocPDFWithLeptons(grid, pdf_in)
  call AllocPDFWithLeptons(grid, pdf_out)
  pdf_in = zero
  ! fill in the QCD part
  pdf_in(:,:ncompmax) = unpolarized_dummy_pdf(xValues(grid))
  ! then add in a photon of some kind
  pdf_in(:,iflv_photon) = alpha_qed_scale_0 * (one - xValues(grid))**4
  pdf_out = pdf_in

  
  call QEDQCDEvolvePDF(dh, qed_split, pdf_out, coupling, coupling_qed,&
       &               Q0, Qfinal, nloop_qcd, nqcdloop_qed)
  !call EvolvePDF(dh, pdf_out(:,ncompmin:ncompmax), coupling, Q0, Qfinal, nloop=nloop_qcd)


  call write_moments(pdf_in, ' in')
  call write_moments(pdf_out,'out')
  call write_pdf(pdf_in,'in')
  call write_pdf(pdf_out,'out')
  
  !! TESTS
  !! * Test against plain qcd evolution, with QED alpha turned to zero
  !!   -> works OK, but watch out for the fact that the tau introduces
  !!      an extra threshold, so that numerically there are differences
  !!      at the evolution precision; reducing du by a factor of 10
  !!      significantly reduces these differences
  !!
  !! * test momentum sum rule, with both nqcdloop_qed values
  !!   -> this seems to work fine
  !!
  !! * compare to photon evolution from another source?
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

  !----------------------------------------------------------------------
  subroutine write_moments(pdf,label)
    real(dp),         intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !-----------------
    real(dp) :: moments(ncompmin:ubound(pdf,dim=2))
    
    write(6,"(a)", advance="no"), "# total momentum "//trim(label)//" (& components)"
    moments = TruncatedMoment(grid, pdf, one)
    write(6,"(40f10.7)") sum(moments), moments

    write(6,"(a)", advance="no"), "# total number "//trim(label)//" (& components)"
    moments(1:6) = TruncatedMoment(grid, pdf(:,1:6)-pdf(:,-1:-6:-1), zero)
    write(6,"(40f11.7)") sum(moments(1:6)), moments(1:6)


  end subroutine write_moments
  
  
  !----------------------------------------------------------------------
  subroutine write_pdf(pdf,label)
    real(dp), intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !----------------------
    real(dp) :: y, yVals(0:grid%ny), last_y
    integer  :: i

    write(6,"(a)"), "# pdf "//trim(label)
    !do i = 0, 100
    !   y = 0.1_dp * i
    !   write(6,*) y, pdf.aty.(y.with.grid)
    !end do
    yVals = yValues(grid)
    do i = 0, grid%ny
       y = yVals(i)
       if (y < 1.000000001_dp * last_y) cycle
       write(6,*) y, pdf(i,:)
       last_y = y
    end do
    write(6,*)
    write(6,*)
    
  end subroutine write_pdf
  
end program test_qed_evol
