!---------------------------------------------------------------
!! provides facilities related to the calculation of the QED coupling
!! with plain 1-loop running (no QCD corrections)
!!  
!-------------------------------------------------------------------
module qed_coupling_module
  use types; use consts_dp
  use warnings_and_errors
  use qcd
  implicit none
  private

  ! use 0.109 to get alpha(mZ) = 1/127.95
  !real(dp), parameter :: m_light_quarks = 0.109_dp ! DECIDE ON THIS
  real(dp), parameter :: m_light_quarks_default = 1.000_dp ! DECIDE ON THIS

  ! lepton masses from from 2022 PDG (plus 2023 update)
  real(dp), parameter, public :: m_electron = 0.51099895000e-3_dp ! +-(15) MeV
  real(dp), parameter, public :: m_muon     = 105.6583755e-3_dp   ! +-(23)
  real(dp), parameter, public :: m_tau      = 1776.86e-3_dp ! +-0.12

  real(dp), parameter, public :: e_dn2 = (one/three)**2
  real(dp), parameter, public :: e_up2 = (two/three)**2

  ! 2022 PDG: summary table + EW review
  ! - 1/alpha(0) = 137.035999180(10)
  ! - 1/alpha(mtau) = 1/133.471 ± 0.007
  ! - 1/alpha(mZ) = 127.951±0.009
  real(dp), parameter :: alpha_qed_scale_0 = one/137.035999180_dp ! +-(10) is uncertainty
  public :: alpha_qed_scale_0
  
  integer, parameter :: n_thresholds = 7

  public :: IQEDThresholdAtQ
  public :: qed_squared_charge
  
  type qed_coupling 
     !private
     real(dp) :: m_light_quarks
     real(dp) :: mc, mb, mt
     integer  :: n_thresholds
     integer  :: nflav(3, 0:n_thresholds) ! first index: 1 = nleptons, 2=ndown, 3=nup
     real(dp) :: thresholds(0:n_thresholds+1)  !
     real(dp) :: b0_values(n_thresholds)
     real(dp) :: alpha_values(n_thresholds)
  end type qed_coupling
  public :: qed_coupling


  ! bring these under one interface to maintain backwards compatibility
  ! but also incorporate a more sensible new interface
  interface InitQEDCoupling
    module procedure  InitQEDCoupling_old, InitQEDCoupling_new
  end interface InitQEDCoupling
  public :: InitQEDCoupling

  interface Value
     module procedure  QEDCouplingValue
  end interface Value
  public :: Value

  interface Delete
     module procedure  DeleteQEDCoupling 
  end interface
  public :: Delete

contains

!----------------------------------------------------------------------
!! Initialises a QED coupling. New interface, which explicitly includes
!! the effective light quark mass, which can be tuned to reproduce the
!! correct high-scale QED behaviour
subroutine InitQEDCoupling_new(coupling, m_light_quarks, m_heavy_quarks, value_at_scale_0)
    use assertions
    type(qed_coupling), intent(out) :: coupling
    real(dp),           intent(in)  :: m_light_quarks, m_heavy_quarks(4:6)
    real(dp), optional, intent(in)  :: value_at_scale_0 ! defaults to alpha_qed_scale_0
    !--------------------------------------------
    real(dp) :: mc, mb, mt
    integer :: i
    
    mc = m_heavy_quarks(4)
    mb = m_heavy_quarks(5)
    mt = m_heavy_quarks(6)
    
    if (mc > m_tau) call wae_error("InitQEDCoupling", "mc > m_tau, but that is not currently supported")
    if (m_light_quarks < m_muon) call wae_error("InitQEDCoupling", "m_light_quarks < m_muon, but that is not currently supported")
    if (m_light_quarks > mc) call wae_error("InitQEDCoupling", "m_light_quarks > mc, but that is not physical (or supported)")

    coupling%n_thresholds = n_thresholds
    
    ! set up the thresholds
    coupling%thresholds(0) = zero
    coupling%thresholds(1) = m_electron
    coupling%thresholds(2) = m_muon
    coupling%thresholds(3) = m_light_quarks
    coupling%thresholds(4) = mc
    coupling%thresholds(5) = m_tau
    coupling%thresholds(6) = mb
    coupling%thresholds(7) = mt
    coupling%thresholds(8) = 1e200_dp

    ! set up the numbers of flavours just above each threshold
    !                       nlept  ndown  nup
    coupling%nflav(:,0) = (/   0,    0,    0 /)
    coupling%nflav(:,1) = (/   1,    0,    0 /)
    coupling%nflav(:,2) = (/   2,    0,    0 /)
    coupling%nflav(:,3) = (/   2,    2,    1 /)
    coupling%nflav(:,4) = (/   2,    2,    2 /)
    coupling%nflav(:,5) = (/   3,    2,    2 /)
    coupling%nflav(:,6) = (/   3,    3,    2 /)
    coupling%nflav(:,7) = (/   3,    3,    3 /)

    ! set up the b0 values by first setting the squared charges
    do i = 1, n_thresholds
       coupling%b0_values(i) = qed_squared_charge(coupling%nflav(1,i),&
            &                                     coupling%nflav(2,i),&
            &                                     coupling%nflav(3,i))
    end do
    ! then include the normalisation
    coupling%b0_values = -two/(6*pi) * coupling%b0_values

    ! finally set up the alpha values at the threshold
    coupling%alpha_values(1) = default_or_opt(alpha_qed_scale_0, value_at_scale_0)
    do i = 2, n_thresholds
       coupling%alpha_values(i) = coupling%alpha_values(i-1)/&
            &  (1 + coupling%b0_values(i-1) * coupling%alpha_values(i-1) &
            &        * two*log(coupling%thresholds(i)/coupling%thresholds(i-1)))
    end do
  end subroutine InitQEDCoupling_new

  !----------------------------------------------------------------------
  !! Initialises a QED coupling. This interface is deprecated because it
  !! uses a default fairly large value for the effective light quark masses,
  !! resulting in a QED coupling that is about 1% too large at mZ
  subroutine InitQEDCoupling_old(coupling, mc, mb, mt, value_at_scale_0)
    type(qed_coupling), intent(out) :: coupling
    real(dp),           intent(in)  :: mc, mb, mt
    real(dp), optional, intent(in)  :: value_at_scale_0 ! defaults to alpha_qed_scale_0
    call InitQEDCoupling_new(coupling, m_light_quarks_default, (/mc,mb,mt/), value_at_scale_0)
  end subroutine InitQEDCoupling_old
      

  !----------------------------------------------------------------------
  !! returns the qed squared charge, including colour factors
  function qed_squared_charge(nleptons, ndown, nup) result(e2sum)
    real(dp)            :: e2sum
    integer, intent(in) :: nleptons, ndown, nup
    e2sum = nleptons + ca * (nup*e_up2 + ndown * e_dn2)
  end function qed_squared_charge
  
  !----------------------------------------------------------------------
  !! returns which threshold we are dealing with at a given scale
  !! (linear time in the number of thresholds, i.e. not super optimal...)
  function IQEDThresholdAtQ(coupling, mu)
    integer :: IQEDThresholdAtQ
    type(qed_coupling), intent(in) :: coupling
    real(dp)          , intent(in) :: mu
    !------------------------------
    integer i

    IQEDThresholdAtQ = 0
    do i = 1, n_thresholds
       if (mu < coupling%thresholds(i)) exit
       IQEDThresholdAtQ = i
    end do
  end function IQEDThresholdAtQ

  !----------------------------------------------------------------------
  !! 
  function QEDCouplingValue(coupling, mu)
    real(dp)                       :: QEDCouplingValue
    type(qed_coupling), intent(in) :: coupling
    real(dp)          , intent(in) :: mu
    !----------------------------------
    integer :: i

    i = IQEDThresholdAtQ(coupling,mu)
    if (i == 0) then
       ! coupling is frozen below the electron mass
       QEDCouplingValue = coupling%alpha_values(1)
    else
       ! run it from the threshold below
       QEDCouplingValue = coupling%alpha_values(i)/&
            &  (1 + coupling%b0_values(i) * coupling%alpha_values(i) &
            &        * two*log(mu/coupling%thresholds(i)))
    end if
  end function QEDCouplingValue

  !----------------------------------------------------------------------
  subroutine DeleteQEDCoupling(coupling)
    type(qed_coupling), intent(in) :: coupling
    ! this routine does nothing but is there for consistency with
    ! interfaces to other objects
  end subroutine DeleteQEDCoupling
  
end module qed_coupling_module

!!$program astest
!!$  use types; use consts_dp
!!$  use qcd_coupling
!!$  implicit none
!!$  integer  :: i
!!$  real(dp) :: Q
!!$
!!$  call InitRunningCoupling(alfas=0.118_dp,nloop=2)
!!$  do i = 0,100
!!$     Q = 91.2_dp**(i/100.0_dp)
!!$     write(6,*) Q, Value(Q)
!!$  end do
!!$  
!!$end program astest
