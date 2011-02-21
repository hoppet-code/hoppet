!! program to test implementation of MSbar masses
!!
!!
program msbar_tests
  use hoppet_v1; use sub_defs_io
  implicit none
  type(running_coupling) :: coupling_pole, coupling_ms
  real(dp) :: asQ0, Q0, asQ_ms, asQ_pole, Q, dt
  real(dp) :: pole_masses(3), msbar_masses(3)
  real(dp) :: test_scales(4) = (/1.0_dp, 3.0_dp, 10.0_dp, 200.0_dp /)
  integer  :: i

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
  
  
end program msbar_tests


