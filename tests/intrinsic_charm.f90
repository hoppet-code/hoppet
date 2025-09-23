!! program to test behaviour with intrinsic charm
!!
!!
program msbar_tests
  use hoppet; use io_utils
  implicit none
  !-- for alpha_s
  type(running_coupling) :: coupling_pole, coupling_ms
  real(dp) :: asMZ, MZ, Q0, Q, dt
  real(dp) :: pole_masses(3)
  real(dp) :: test_scales(6) = (/1.0_dp, 1.3_dp, 1.5_dp, 3.0_dp, 10.0_dp, 200.0_dp /)
  integer  :: i
  !-- for PDFs
  type(grid_def) :: grid
  type(dglap_holder) :: dh
  type(pdf_table)    :: table_pole, table_ms
  real(dp), pointer  :: pdf_init(:,:), pdf_final(:,:)
  real(dp) :: x, pdf_pole(-6:6), pdf_ms(-6:6)
  integer  :: iflv


  dt = dble_val_opt("-dt", 0.2_dp)
  call SetDefaultCouplingDt(dt)

  Q0 = 1.3_dp
  
  asMZ=0.118_dp
  MZ=91.2_dp
  pole_masses(1:3) = (/1.414213563_dp, 4.5_dp, 175.0_dp /)
  
  ! get a coupling with pole masses
  call InitRunningCoupling(coupling_pole, asMZ, MZ, nloop=3, quark_masses = pole_masses)

  

  !======================================================================
  ! now set up 
  call SetDefaultEvolutionDu(0.05_dp)
  call InitGridDef(grid,ymax=10.0_dp, dy=0.1_dp,order=-5)
  call InitDglapHolder(grid,dh,nloop=3,nflo=3,nfhi=6)
  call AllocPDF(grid, pdf_init)
  call AllocPDF(grid, pdf_final)
  pdf_init = unpolarized_dummy_pdf(xValues(grid))

  call AllocPDFTable(grid, table_pole, 1.0_dp, 1e4_dp)
  call AddNfInfoToPDFTable(table_pole, coupling_pole)

  !call PreEvolvePDFTable(table_pole, Q0, dh, coupling_pole)
  !call EvolvePDFTable(table_pole, pdf_init)

  call EvolvePDFTable(table_pole, Q0, pdf_init, dh, coupling_pole)

  pdf_final = pdf_init
  call EvolvePDF(dh, pdf_final, coupling_pole, Q0, 1.5_dp)

  ! side test
  x = 1e-1_dp
  iflv = int_val_opt("-iflv",4)
  
  do i = 1, size(test_scales)
    Q = test_scales(i)
    call EvalPDFTable_xQ(table_pole, x, Q, pdf_pole)
    write(6,*) 'Q, pdf_pole(iflv), pdf_ms(iflv):', &
         & Q, pdf_pole(iflv), pdf_final(:,iflv).atx.(x.with.grid)
  end do



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

    pdf(:, iflv_c) = dbar

  end function unpolarized_dummy_pdf

end program msbar_tests


