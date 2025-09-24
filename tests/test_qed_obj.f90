program test_qed_obj
  use hoppet; use convolution
  use qed_objects
  use qed_coupling_module
  use io_utils
  implicit none
  type(grid_def)     :: grid
  type(qed_split_mat) :: qed_split
  real(dp) :: ymax = 25.0_dp, dy = 0.10_dp
  real(dp), pointer :: pdf(:,:), dpdf(:,:), moments(:)
  integer :: ncompmaxQED = 11 ! 8 without leptons, 11 with leptons
  type(qed_coupling) :: qed_alpha
  
  call InitGridDefDefault(grid, dy, ymax)
  call InitQEDSplitMat(grid, qed_split)
  !call qed_split%SetNf(3,3,2)
  call QEDSplitMatSetNf(qed_split, 3,3,2)
  
  call AllocGridQuant(grid,  pdf, -6, ncompmaxQED)
  call AllocGridQuant(grid, dpdf, -6, ncompmaxQED)
  allocate(moments(ncompmin:ncompmaxQED))
  
  pdf = zero
  pdf(:,ncompmin:ncompmax) = unpolarized_dummy_pdf(xValues(grid))
  pdf(:,8) = (one - xValues(grid))**4
  if (ncompmaxQED >= 9) pdf(:,9) = (one - xValues(grid))**5
  
  dpdf = qed_split%lo * pdf
  moments = TruncatedMoment(grid, dpdf, one)
  write(6,*) '*********** LO QED splitting checks **************'
  write(6,*) '1st moment of dflavours = ', moments
  write(6,*) 'sum of dflavours = ', sum(moments)
  write(6,*)
  moments(1:6) = TruncatedMoment(grid, dpdf(:,1:6)-dpdf(:,-1:-6:-1), zero)
  write(6,*) 'zeroth moment of dflav-dantiflav = ', moments(1:6)
  write(6,*)

  dpdf = qed_split%nlo * pdf
  moments = TruncatedMoment(grid, dpdf, one)
  write(6,*) '*********** NLO QED splitting checks **************'
  write(6,*) '1st moment of dflavours = ', moments
  write(6,*) 'sum of dflavours = ', sum(moments)
  write(6,*)
  moments(1:6) = TruncatedMoment(grid, dpdf(:,1:6)-dpdf(:,-1:-6:-1), zero)
  write(6,*) 'zeroth moment of dflav-dantiflav = ', moments(1:6)
  write(6,*)


  ! now check the pure photon-branching parts
  call InitQEDCoupling(qed_alpha, 1.5_dp, 4.5_dp, 173.2_dp)
  pdf(:,-6:6) = zero
  if (ncompmaxQED == ncompmaxLeptons) pdf(:,9:11) = zero
  dpdf = qed_split%lo * pdf
  moments(8) = TruncatedMoment(grid, dpdf(:,8), two)
  moments(7) = TruncatedMoment(grid,  pdf(:,8), two)
  write(6,*) 'LO fractional dpdf/pdf for photon is', moments(8)/moments(7)
  write(6,*) 'LO qed beta function is ', qed_alpha%b0_values(6) * twopi
  dpdf = qed_split%nlo * pdf
  moments(8) = TruncatedMoment(grid, dpdf(:,8), two)
  moments(7) = TruncatedMoment(grid,  pdf(:,8), two)
  write(6,*) 'NLO fractional dpdf/pdf for photon is', moments(8)/moments(7)
  !write(6,*) 'NLO qed beta function is ', qed_alpha%b0_values(6) * twopi
  write(6,*)
contains

  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolazrized evolution
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
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)
    pdf(:,iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
        
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf

end program test_qed_obj

