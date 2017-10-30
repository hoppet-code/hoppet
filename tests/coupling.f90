program coupling
  use types; use consts_dp
  use qcd_coupling
  implicit none

  real(dp), parameter :: mc_pole = 1.67_dp
  real(dp), parameter :: mb_pole = 4.78_dp
  real(dp), parameter :: mt_pole = 173.21_dp

  real(dp), parameter :: mc_msbar = 1.275_dp
  real(dp), parameter :: mb_msbar = 4.18_dp

  type(running_coupling) :: this_coupling

  real(dp), parameter :: mz = 91.1876_dp
  real(dp) :: alphas_mz = 0.118_dp
  real(dp) :: mu_vals(5) = (/ 1.1_dp, 2.0_dp, 6.0_dp, mz, 1000.0_dp /)
  real(dp) :: mu
 
  integer  :: i_mu

  call InitRunningCoupling(this_coupling, alphas_mz, Q=mz, nloop = 4, &
       &                   quark_masses = (/mc_pole, mb_pole, mt_pole /),&
       &                   masses_are_MSbar = .false.)

  do i_mu = 1, size(mu_vals)
     mu = mu_vals(i_mu)
     write(6,*) mu, Value(this_coupling, mu)
  end do

  write(6,*)
  call delete(this_coupling)
  call InitRunningCoupling(this_coupling, alphas_mz, Q=mz, nloop = 4, &
       &                   quark_masses = (/mc_msbar, mb_msbar, 161.0_dp /),&
       &                   masses_are_MSbar = .true.)
  do i_mu = 1, size(mu_vals)
     mu = mu_vals(i_mu)
     write(6,*) mu, Value(this_coupling, mu)
  end do
  
  
end program coupling
