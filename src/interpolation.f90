! $Id: interpolation.f90,v 1.3 2002/12/30 17:20:26 gsalam Exp $
module interpolation
  use types
  implicit none
  private

  public :: uniform_interpolation_weights
  public :: general_interpolation_weights
  public :: interpolation_grid

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
    real(dp), intent(in)   :: x
    real(dp), intent(out)  :: weights(0:)
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


!--------------------------------------------------------------------
  ! returns the weights for a general, non-uniform interpolation
  ! the x grid is now another argument of the procedure
  !
  ! formula is weight(i) = [Prod_{j/=i} (x-xgrid(j))] 
  !    / [Prod_{j=0}^{i-1} (xgrid(i)-xgrid(j))] 
  !    / [Prod_{j=i+1}^n (xgrid(i)-xgrid(j))]
  !
  ! Algorithm uses instead if one
  ! of the x-i==0 in the result is just 0,...,0,,1,0...,0
  !
  !
  ! For simplicity of caching, n == ubound(weights), is limited .le. nmax
  subroutine general_interpolation_weights(x,xgrid,weights)
    use warnings_and_errors; use consts_dp
    real(dp), intent(in)  :: x,xgrid(0:)
    real(dp), intent(out) :: weights(0:)
    !-----------------------------------------
    !integer,  parameter :: nmax = 9
    integer             :: n, i, ngrid, j
    real(dp)            :: normalisation_general(0:ubound(weights,dim=1))
    real(dp)            :: dists(0:ubound(weights,dim=1)), prod
    real(dp)            :: dist

    !Check dimensions    
    n = ubound(weights,dim=1)
    ngrid= ubound(xgrid, dim=1)
    if (n.ne.ngrid) call wae_error('general_interpolation_weights',&
         &'ubound of weights different ubound xgrid:',intval=n)
    
    !-- calculate normalization
    do i = 0, n
       normalisation_general(i) = one
       do j = 0, i-1
          dist  = xgrid(i) - xgrid(j)
          normalisation_general(i) = normalisation_general(i) * dist
       end do
       do j = i+1, n
          dist  = xgrid(i) - xgrid(j)
          normalisation_general(i) = normalisation_general(i) * dist
       end do
       ! print *, i,normalisation_general(i)
    enddo
    !-- calculate inverse weight "normalisation_generals"
    normalisation_general(:n) = one / normalisation_general(:n)

    ! Compute interpolation weights
    do i = 0, n
       dists(i) = x - xgrid(i)
       if (dists(i) == zero) then
          weights(:) = zero
          weights(i) = one
          return
       end if
    end do
    prod = product(dists)
    do i = 0, n
       weights(i) = prod * normalisation_general(i) / dists(i)
    end do
    
  end subroutine general_interpolation_weights

!-----------------------------------------


  !! ! Alternative form
  !! function interpolation_ixlo(x, xgrid, n_interp) result(ixlo)
  !!   real(dp), intent(in) :: x, xgrid(0:)
  !!   integer              :: ixlo
  !!   !...
  !!   ! set ixlo
  !!   stop
  !! end function interpolation_ixlo
  

  ! Select the suitable grid for the general interpolation
  ! starting from the initial large grid
  subroutine interpolation_grid(n_interp,x, xgrid, ix_interp, x_interp)
    use consts_dp
    implicit none
    
    integer, intent(in) :: n_interp
    real(dp), intent(in)  :: x,xgrid(:)
    
    real(dp), intent(out) :: x_interp(0:)
    integer, intent(out) :: ix_interp(0:)
    
    integer ix,jx, nxgrid,i_interp

    nxgrid = ubound(xgrid,dim=1)

    ! Find the associated interpolation grid
    ! [ Either binary search in general, or in particular
    ! make use of knowledge of how grid has been generated ]
    ! Binary grid is the simplest technique
    ! Only assumption: ordering of grid points

    ! Binary search
    jx=1
    do ix=1,nxgrid-1 ! Exclude last point so that jx <= nxgrid
       if( x.lt. xgrid(ix) ) then
          jx = jx + 1
       else
          exit ! Exit the loop
       endif
    enddo
    !write(6,*) "x, jx, nxgrid = ",x,jx,nxgrid

 
    ! Now define grid for interpolation
    ! Different according to the value of x
   
    if( jx.le.n_interp ) then
       ! Left region
       ! GPS illustration
       !x_interp(:) = xgrid(i_interp+1: )
       do i_interp = 0, n_interp
          ix_interp(i_interp) = i_interp + 1
          x_interp(i_interp) = xgrid( ix_interp(i_interp) )
       enddo
      
    else if( jx.lt.( nxgrid - n_interp/2 ) ) then
        ! Central region
       do i_interp = 0, n_interp
          ix_interp(i_interp) =  jx + i_interp - n_interp/2
          x_interp(i_interp) = xgrid( ix_interp(i_interp) )
       enddo
       
    else
       ! Right region
       do i_interp = 0, n_interp
          ix_interp(i_interp) =  nxgrid + i_interp - n_interp
          x_interp(i_interp) = xgrid( ix_interp(i_interp) )
       enddo

    endif

!    write(6,*) "i_interp, ix_interp(i_interp),  x_interp(i_interp)"
!    do i_interp = 0,  n_interp
!       write(6,*) i_interp, ix_interp(i_interp),  x_interp(i_interp)
!    enddo
    
    
  end subroutine interpolation_grid
  

end module interpolation

