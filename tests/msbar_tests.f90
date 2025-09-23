!! program to test implementation of MSbar-mass thresholds
!!
!!
program msbar_tests
  use hoppet; use io_utils
  implicit none
  !-- for alpha_s
  type(running_coupling) :: coupling_pole, coupling_ms
  real(dp) :: asQ0, Q0, asQ_ms, asQ_pole, Q, dt
  real(dp) :: pole_masses(3), msbar_masses(3)
  real(dp) :: test_scales(4) = (/1.0_dp, 3.0_dp, 10.0_dp, 200.0_dp /)
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

  Q0 = 3.0_dp
  asQ0 = dble_val_opt('-asQ0',0.255786_dp)

  pole_masses(1:3) = (/1.414213563_dp, 4.5_dp, 175.0_dp /)
  
  ! get a coupling with pole masses
  call InitRunningCoupling(coupling_pole, asQ0, Q0, nloop=3, quark_masses = pole_masses,&
       &                   masses_are_MSbar = .false.)
  
  ! get the conversion of pole masses to MSbar masses (requires an initialised coupling)
  do i = 1, 3
    msbar_masses(i) = pole_masses(i)*(1 - four/three/pi*value(coupling_pole, pole_masses(i)))
    write(6,*) 'pole, MSbar masses are:', pole_masses(i), msbar_masses(i)
  end do

  ! get the coupling with the MSbar masses
  call InitRunningCoupling(coupling_ms, asQ0, Q0, nloop=3, quark_masses = msbar_masses,&
       &                   masses_are_MSbar = .true.)

  ! verify that the difference between the two couplings scales as as^4
  ! (since matching is provided at relative O(as^2), i.e. absolute as^3)
  do i = 1, size(test_scales)
    Q = test_scales(i)
    asQ_pole = value(coupling_pole, Q)
    asQ_ms = value(coupling_ms, Q)
    write(6,*) 'Q, as(Q,pole), as(Q,msbar):', &
         & Q, asQ_pole, asQ_ms, (asQ_pole-asQ_ms)/asQ0**4
  end do
  

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
  !call EvolvePDFTable(table_pole, Q0, pdf_init, dh, coupling_pole)
  call PreEvolvePDFTable(table_pole, Q0, dh, coupling_pole)

  call AllocPDFTable(grid, table_ms, 1.0_dp, 1e4_dp)
  call AddNfInfoToPDFTable(table_ms, coupling_ms)
  !call EvolvePDFTable(table_ms, Q0, pdf_init, dh, coupling_ms)
  call PreEvolvePDFTable(table_ms, Q0, dh, coupling_ms)


  call EvolvePDFTable(table_pole, pdf_init)
  call EvolvePDFTable(table_ms, pdf_init)
  
  x = 1e-1_dp
  iflv = int_val_opt("-iflv",0)
  
  do i = 1, size(test_scales)
    Q = test_scales(i)
    call EvalPDFTable_xQ(table_pole, x, Q, pdf_pole)
    call EvalPDFTable_xQ(table_ms  , x, Q, pdf_ms)
    write(6,*) 'Q, pdf_pole(iflv), pdf_ms(iflv):', &
         & Q, pdf_pole(iflv), pdf_ms(iflv), (pdf_pole(iflv)-pdf_ms(iflv))/asQ0**3
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
  end function unpolarized_dummy_pdf

end program msbar_tests


