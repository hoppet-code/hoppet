program test_qed
  use types; use consts_dp
  use convolution
  use qed_coupling_module
  use qed_splitting_functions
  use splitting_functions
  use qcd_qed_coupling
  use qcd_coupling
  implicit none
  type(grid_def)     :: grid
  type(qed_coupling) :: qed_alpha
  type(grid_conv)    :: Pqq_01, Pqy_01, Pyq_01, Pyy_01
  type(grid_conv)    :: Pqg_11, Pyg_11, Pgg_11
  type(grid_conv)    :: PqqV_11, PqqbarV_11, Pgq_11, Pyq_11
  real(dp), pointer  :: res(:), source(:)
  real(dp) :: moment

  real(dp) :: scales(9) = (/m_electron/two, 0.01_dp, 0.1_dp, 1.0_dp, 1.68_dp, 10.0_dp, 91.2_dp, 100.0_dp, 1000.0_dp/)
  real(dp) :: quark_masses(4:6) = (/ 1.67_dp, 4.78_dp, 173.21_dp/)
  real(dp) :: cpls(2), alphas, mu
  integer  :: i, n
  type(QCDQEDCoupling) :: qcd_qed
  type(running_coupling) :: qcd_cpl
  
  write(6,'(a)') "========== testing QED coupling ========================="
  call InitQEDCoupling(qed_alpha, quark_masses(4), quark_masses(5), quark_masses(6))
  write(6,*) IQEDThresholdAtQ(qed_alpha, 0.2_dp)

  qcd_qed = new_QCDQEDCoupling(0.118_dp, 1/128.6496065_dp, 91.2_dp, &
       & '01 02 03 10 11',&
       & quark_masses, &
       & (/ qcd_mass_scheme_pole, qcd_mass_scheme_pole, qcd_mass_scheme_pole /) )
  write(6,*) coupling_array(qcd_qed%Coupling, 91.2_dp)
  call InitRunningCoupling(0.118_dp, 91.2_dp, nloop=3, quark_masses = quark_masses)
  
  write(6,'(10a13)') 'mu','1/alpha','1/al(new)','diff','alphas','alphas(new)','diff'
  do i = 1, size(scales)
     mu = scales(i)
  !n=200; do i = 0, n
  !   mu = exp(i*log(300.0_dp)/n)
     cpls = coupling_array(qcd_qed%Coupling, mu)
     alphas = zero
     if (mu > 0.5_dp) alphas = Value(mu)
     write(6,'(f13.6,3f13.7,3f13.9)') mu, &
          & one/Value(qed_alpha, mu), 1/cpls(2), one/Value(qed_alpha, mu)-1/cpls(2), &
          & alphas, cpls(1), alphas-cpls(1)
  end do
  !stop !************************************
  
  write(6,'(a)') "========== testing sum rules ========================="
  call setup_grid_def()

  call InitGridConv(grid, Pqq_01, sf_Pqq_01)
  call InitGridConv(grid, Pyq_01, sf_Pyq_01)
  call InitGridConv(grid, Pqy_01, sf_Pqy_01)
  call InitGridConv(grid, Pyy_01, sf_Pyy_01)

  call InitGridConv(grid, Pqg_11, sf_Pqg_11)
  call InitGridConv(grid, Pyg_11, sf_Pyg_11)
  call InitGridConv(grid, Pgg_11, sf_Pgg_11)

  call InitGridConv(grid, PqqV_11, sf_PqqV_11)
  call InitGridConv(grid, PqqbarV_11, sf_PqqbarV_11)
  call InitGridConv(grid, Pgq_11, sf_Pgq_11)
  call InitGridConv(grid, Pyq_11, sf_Pgq_11)
  
  call AllocGridQuant(grid, source)
  call AllocGridQuant(grid, res)

  moment = one
  write(6,*) "looking at moment=",moment
  source = exp(moment * yValues(grid))
  res = Pqq_01 * source
  write(6,*) "Pqq_01 = ", res(grid%ny)/source(grid%ny)
  res = Pyq_01 * source
  write(6,*) "Pyq_01 = ", res(grid%ny)/source(grid%ny)
  res = Pqy_01 * source
  write(6,*) "Pqy_01 = ", res(grid%ny)/source(grid%ny)
  res = Pyy_01 * source
  write(6,*) "Pyy_01 = ", res(grid%ny)/source(grid%ny)
  write(6,*)

  
  res = Pqg_11 * source
  write(6,*) "Pqg_11 = ", res(grid%ny)/source(grid%ny)
  res = Pyg_11 * source
  write(6,*) "Pyg_11 = ", res(grid%ny)/source(grid%ny)
  res = Pgg_11 * source
  write(6,*) "Pgg_11 = ", res(grid%ny)/source(grid%ny)
  res = (two*(Pqg_11*source) + Pyg_11 * source + Pgg_11 * source)
  write(6,*) "(2*Pqg + Pyg + Pgg)_11 = ", res(grid%ny)/source(grid%ny)
  write(6,*)
  
  res = PqqV_11 * source
  write(6,*) "PqqV_11 = ", res(grid%ny)/source(grid%ny)
  res = PqqbarV_11 * source
  write(6,*) "PqqbarV_11 = ", res(grid%ny)/source(grid%ny)
  res = Pgq_11 * source
  write(6,*) "Pgq_11 = ", res(grid%ny)/source(grid%ny)
  res = Pyq_11 * source
  write(6,*) "Pyq_11 = ", res(grid%ny)/source(grid%ny)
  res = (PqqV_11*source + PqqbarV_11*source + Pgq_11*source + Pyq_11*source)
  write(6,*) "(PqqV+PqqbarV+Pgq+Pyq)_11 = ", res(grid%ny)/source(grid%ny)

  
contains
  subroutine setup_grid_def()
    real(dp) :: ymax, dy
    integer  :: order
    type(grid_def) :: gdarray(4)
    ! set up parameters for grid
    order = -6
    ymax  = 20.0_dp  ! you may see significant violations with too small a ymax
    dy    = 0.1_dp
    
    !call SetDefaultEvolutionDu(dy/3.0_dp)  ! generally a good choice
    
    ! set up the grid itself -- we use 4 nested subgrids
    call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
    call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
    call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
    call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)
    
  end subroutine setup_grid_def

  
end program test_qed
