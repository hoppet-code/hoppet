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

  real(dp), parameter :: m_light_quarks = 1.000_dp ! DECIDE ON THIS

  ! lepton masses from from 2014 PDG
  real(dp), parameter, public :: m_electron = 0.510998928e-3_dp ! +-(11) on last digits
  real(dp), parameter, public :: m_muon     = 105.6583715e-3_dp ! +-(35) on last digits
  real(dp), parameter, public :: m_tau      = 1776.82e-3_dp ! +-(16) on last digits

  real(dp), parameter, public :: e_dn2 = (one/three)**2
  real(dp), parameter, public :: e_up2 = (two/three)**2

  ! from the PDG rpp2014-rev-phys-constants
  real(dp), parameter :: alpha_qed_scale_0 = one/137.035999074_dp ! +-(44) is uncertainty
  public :: alpha_qed_scale_0
  
  integer, parameter :: n_thresholds = 7

  public :: IQEDThresholdAtQ
  public :: qed_squared_charge
  
  type qed_coupling 
     !private
     real(dp) :: mc, mb, mt
     integer  :: n_thresholds
     integer  :: nflav(3, 0:n_thresholds) ! first index: 1 = nleptons, 2=ndown, 3=nup
     real(dp) :: thresholds(0:n_thresholds+1)  !
     real(dp) :: b0_values(n_thresholds)
     real(dp) :: alpha_values(n_thresholds)
  end type qed_coupling
  public :: qed_coupling

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
  !! Initialises a QED coupling
  subroutine InitQEDCoupling(coupling, mc, mb, mt, value_at_scale_0)
    use assertions
    type(qed_coupling), intent(out) :: coupling
    real(dp),           intent(in)  :: mc, mb, mt
    real(dp), optional, intent(in)  :: value_at_scale_0 ! defaults to alpha_qed_scale_0
    !--------------------------------------------
    integer :: i
    
    if (mc > m_tau) call wae_error("InitQEDCoupling", "mc > m_tau")

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
  end subroutine InitQEDCoupling
  

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
