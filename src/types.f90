!======================================================================
! Basic constants needed everywhere...
!
! GPS 16/10/96
! $Id: types.f90,v 1.2 2002/01/25 16:20:48 gsalam Exp $
!======================================================================



!----------------------------------------------------------------------
! The next two are maybe to be used only when it's time for a general
! rewrite?
module types
  implicit none
  integer, parameter  :: dp = kind(1.0d0), sp = kind(1.0)
end module types
!----------------------------------------------------------------------
module consts_dp
  use types
  implicit none
  private

  real(dp), public, parameter :: pi =&
       & 3.141592653589793238462643383279502884197_dp
  real(dp), public, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_dp
  real(dp), public, parameter :: zeta2 =&
       & 1.644934066848226436472415166646025189219_dp
  real(dp), public, parameter :: zeta3 =&
       & 1.202056903159594285399738161511449990765_dp
  real(dp), parameter, public :: pisq =&
       & 9.869604401089358618834490999876151135314_dp
  real(dp), parameter, public :: eulergamma =&
       & 0.577215664901532860606512090082402431042_dp
  real(dp), parameter, public :: ln2 =&
       & 0.693147180559945309417232121458176568076_dp
  real(dp), public, parameter :: half = 0.5_dp, two = 2.0_dp
  real(dp), public, parameter :: zero = 0.0_dp, one = 1.0_dp
  real(dp), public, parameter :: three = 3.0_dp, four = 4.0_dp

end module consts_dp

