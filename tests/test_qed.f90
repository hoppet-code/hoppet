program test_qed
  use types; use consts_dp
  use convolution
  use qed_coupling_module
  use qed_splitting_functions
  use splitting_functions
  implicit none
  type(grid_def)     :: grid
  type(qed_coupling) :: qed_alpha
  type(grid_conv)    :: Pqq_01, Pqy_01, Pyq_01, Pyy_01
  type(grid_conv)    :: Pqg_11, Pyg_11, Pgg_11
  type(grid_conv)    :: PqqV_11, PqqbarV_11, Pgq_11, Pyq_11
  real(dp), pointer  :: res(:), source(:)
  real(dp) :: moment

  write(6,'(a)') "========== testing QED coupling ========================="
  call InitQEDCoupling(qed_alpha, 1.5_dp, 4.5_dp, 173.0_dp)
  write(6,*) IQEDThresholdAtQ(qed_alpha, 0.2_dp)

  
  write(6,*) 'alpha QED at 0       = ', one/Value(qed_alpha, 0.0_dp)
  write(6,*) 'alpha QED at MZ      = ', one/Value(qed_alpha, 91.2_dp)
  write(6,*) 'alpha QED at 0.1 GeV = ', one/Value(qed_alpha, 0.1_dp)
  write(6,*) 'alpha QED at 1 GeV   = ', one/Value(qed_alpha, 1.0_dp)
  write(6,*) 'alpha QED at 10 GeV  = ', one/Value(qed_alpha, 10.0_dp)
  write(6,*) 'alpha QED at 100 GeV = ', one/Value(qed_alpha, 100.0_dp)
  write(6,*) 'alpha QED at 1 TeV   = ', one/Value(qed_alpha, 1000.0_dp)
  
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
