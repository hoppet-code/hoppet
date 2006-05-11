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
    real(dp)            :: dists(0:ubound(weights,dim=1)), prod
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
    
    do i = 0, n
       dists(i) = x - i
       if (dists(i) == zero) then
          weights(:) = zero
          weights(i) = one
          return
       end if
    end do
    prod = product(dists)
    do i = 0, n
       weights(i) = prod * normalisation(i,n) / dists(i)
    end do
    
  end subroutine uniform_interpolation_weights
  

end module interpolation

