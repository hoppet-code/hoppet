! $Id: interpolation.f90,v 1.3 2002/12/30 17:20:26 gsalam Exp $
module interpolation
  use types
  implicit none
  private

  public :: uniform_interpolation_weights
contains


  !--------------------------------------------------------------------
  ! returns the weights for the uniform interpolation,
  ! where the first entry of weights corresponds to x=1
  ! and the spacing is 1.
  !
  ! formula should be weight(i) = [Prod_{j/=i} (x-j)] (-1)^(n-1) / (i! (n-i)!)
  !
  ! Algorithm uses instead [Prod_{j} (x-j)] / (x-i) unless one
  ! of the x-i==0 in which case the result is just 0,...,0,,1,0...,0
  !
  ! For simplicity of caching, n == ubound(weights), is limited .le. nmax
  subroutine uniform_interpolation_weights(x,weights)
    use warnings_and_errors; use consts_dp
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: weights(0:)
    !-----------------------------------------
    integer,  parameter :: nmax = 9
    !                                           order=n
    real(dp), save      :: normalisation(0:nmax,0:nmax) = zero
    real(dp)            :: dists(0:nmax)
    real(dp)            :: prod
    integer             :: n, i
    
    n = ubound(weights,dim=1)
    if (n > nmax) call wae_error('uniform_interpolation_weights',&
         &'ubound of weights is too large:',intval=n)

    !-- intialise once for each n
    if (normalisation(0,n) == zero) then
       !-- calculate factorial
       normalisation(0,n) = one
       do i = 1, n
         normalisation(0,n) = normalisation(0,n) * (-i)
       end do
       !-- calculate inverse weight "normalisations"
       do i = 1, n
         normalisation(i,n) = (normalisation(i-1,n) * i)/(i-1-n)
       end do
       normalisation(:n,n) = one / normalisation(:n,n)
    end if
    
    prod = one
    do i = 0, n
      dists(i) = x - i
      if (dists(i) == zero) then
        weights(:) = zero
        weights(i) = one
        return
      end if
      prod = prod * dists(i)
    end do
    !prod = product(dists)
    do i = 0, n
      weights(i) = prod * normalisation(i,n) / dists(i)
    end do
    
  end subroutine uniform_interpolation_weights
  
end module interpolation

!! Module with coefficients for interpolation of uniformly spaced points. 
!! 
!! With Nth order interpolation, points at x=0,...,N, 
!! interp_coeffsN(i) should multiply $f(x=i) \prod_{j=0,j\neq i}^N (x-j)$
!!
!! When speed is critical, using these hard-coded coefficients will be
!! faster than calling uniform_interpolation_weights. Faster still
!! is to use the precomputed values directly in your code, avoiding
!! any overhead from function calls.
module interpolation_coeffs
  use types
  implicit none
  private

  !! Interp_coeffsN(i) should multiply $f(x=i) \prod_{j=0,j\neq i}^N (x-j)$
  real(dp), parameter :: interp_coeffs1(0:1) = [-1.0_dp, 1.0_dp]
  real(dp), parameter :: interp_coeffs2(0:2) = [ 0.5_dp, -1.0_dp, 0.5_dp ]
  real(dp), parameter :: interp_coeffs3(0:3) = [ -1.0_dp/6.0_dp, 1.0_dp/2.0_dp, -1.0_dp/2.0_dp, 1.0_dp/6.0_dp ]
  real(dp), parameter :: interp_coeffs4(0:4) = [ 1.0_dp/24.0_dp, -1.0_dp/6.0_dp, 1.0_dp/4.0_dp,  1.0_dp/6.0_dp,  -1.0_dp/24.0_dp ]

  public :: interp_coeffs1, interp_coeffs2, interp_coeffs3, interp_coeffs4
  public :: fill_interp_weights1, fill_interp_weights2, fill_interp_weights3, fill_interp_weights4

contains

  pure subroutine fill_interp_weights1(x, weights)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: weights(0:1)
    !-----------------------------------------
    weights(0) = (x-1.0_dp) * interp_coeffs1(0)
    weights(1) = (x       ) * interp_coeffs1(1)
  end subroutine fill_interp_weights1

  pure subroutine fill_interp_weights2(x, weights)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: weights(0:2)
    !-----------------------------------------
    real(dp) :: xm1,xm2
    xm1 = x-1.0_dp
    xm2 = x-2.0_dp
    weights(0) = xm1*xm2 * interp_coeffs2(0)
    weights(1) = x  *xm2 * interp_coeffs2(1)
    weights(2) = x  *xm1 * interp_coeffs2(2)
  end subroutine fill_interp_weights2

  pure subroutine fill_interp_weights3(x, weights)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: weights(0:3)
    !-----------------------------------------
    real(dp) :: xm1,xm2,xm3
    xm1 = x-1.0_dp
    xm2 = x-2.0_dp
    xm3 = x-3.0_dp
    weights(0) = xm1*xm2*xm3 * interp_coeffs3(0)
    weights(1) = x  *xm2*xm3 * interp_coeffs3(1)
    weights(2) = x  *xm1*xm3 * interp_coeffs3(2)
    weights(3) = x  *xm1*xm2 * interp_coeffs3(3)
  end subroutine fill_interp_weights3

  pure subroutine fill_interp_weights4(x, weights)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: weights(0:4)
    !-----------------------------------------
    real(dp) :: xm1,xm2,xm3,xm4, xm34, xm01
    xm1 = x-1.0_dp
    xm2 = x-2.0_dp
    xm3 = x-3.0_dp
    xm4 = x-4.0_dp
    ! on M2Pro-gfortran15-O3, pre-computing just these two seems to give
    ! the best speed
    xm34 = xm3 * xm4
    xm01 = x * xm1
    weights(0) = xm1*xm2*xm34 * interp_coeffs4(0)
    weights(1) = x  *xm2*xm34 * interp_coeffs4(1)
    weights(2) = xm01*xm34    * interp_coeffs4(2)
    weights(3) = xm01*xm2*xm4 * interp_coeffs4(3)
    weights(4) = xm01*xm2*xm3 * interp_coeffs4(4)
  end subroutine fill_interp_weights4

  ! this variant is a bit slower, by about 1ns on M2Pro-gfortran15-O3
  !pure subroutine fill_interp_weights4(x, weights)
  !  real(dp), intent(in)  :: x
  !  real(dp), intent(out) :: weights(0:4)
  !  !-----------------------------------------
  !  real(dp) :: xm(0:4), prodsl(0:4), prodsr(0:4)
  !  integer  :: i
  !  xm = x - [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  !  prodsl(0) = 1.0_dp
  !  prodsr(4) = 1.0_dp
  !  do i = 1,4
  !    prodsl(i) = prodsl(i-1) * xm(i-1)
  !    prodsr(4-i) = prodsr(5-i) * xm(5-i)
  !  end do
  !  weights = prodsl * prodsr * interp_coeffs4
  !end subroutine fill_interp_weights4

end module interpolation_coeffs
