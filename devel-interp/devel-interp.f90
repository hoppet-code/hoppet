!-----------------------------------------------------------
!
! Devel program to test the new general interpolation module
! of hoppet designed to work with interpolated splitting functions
!
! JRC, 06/08
!
!------------------------------------------------------------


program test_hoppet

! Use the main module
use hoppet_v1
use interpolation


implicit none

integer, parameter :: n=3
integer :: i,k

real(dp) :: x, interp
real(dp) :: weights(0:n),weights2(0:n)
real(dp) :: xgrid(0:n)
real(dp) :: ygrid(0:n)

! Check same results as uniform interpolation in
! the uniform limit

! Set uniform grid
! The first point corresponds to x=1
print *,"Interpolation grid"
open(unit=10,status="unknown",file="grid.res")
do i = 0 , n
   xgrid(i) = i
   ygrid(i) = exp(dfloat(-i))
   print *,i,xgrid(i),ygrid(i)
   write(10,*) xgrid(i),ygrid(i)
enddo
close(10)

! Compare the two methods
x=1.5
print *,x
call uniform_interpolation_weights(x,weights)
call general_interpolation_weights(x,xgrid,weights2)

print *,"Interpolation weights"
do i = 0 , n
   print *,i,weights(i),weights2(i)
enddo

! Check general interpolation in the uniform case
print *,"Save interpolated function"
open(unit=10,status="unknown",file="interp.res")
do k=0,100
   x=-1+5*dfloat(k)/100

   interp=0
   call general_interpolation_weights(x,xgrid,weights2)
   do i = 0 , n
      interp = interp + ygrid(i)*weights2(i)
   enddo
   
   write(10,*) x,interp

enddo
close(10)

! Check interpolating polynomial in the general,
! non uniform case
print *,"Interpolation grid"
open(unit=10,status="unknown",file="grid2.res")
do i = 0 , n
   xgrid(i) = -log( 1.0/(i+1) )
   ygrid(i) = exp( dfloat(-i) )
   print *,i,xgrid(i),ygrid(i)
   write(10,*) xgrid(i),ygrid(i)
enddo
close(10)

! Check general interpolation in the uniform case
print *,"Save interpolated function"
open(unit=10,status="unknown",file="interp2.res")
do k=0,1000
   x=-1+3*dfloat(k)/1000

   interp=0
   call general_interpolation_weights(x,xgrid,weights)
   do i = 0 , n
      interp = interp + ygrid(i)*weights(i)
   enddo
   
   write(10,*) x,interp

enddo
close(10)



end program test_hoppet

!-----------------------------------------------------------


