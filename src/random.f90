!--------------------------------------------------------------------
!! contains the random number generator; use l'Ecuyer's method which
!! is sufficient for simple uses (want something that above all is
!! quick) . Original code comes via herwig.
!!
module random
  use types
  implicit none
  private

  public :: ran,rangen, SetSeed, GetSeed

  integer, save :: iseed(2) = (/123450,678900/)
contains
  function ran()
    real(dp) :: ran
    real(dp) :: rann(1)
    !call random_number(ran)
    call rangen(rann)
    ran = rann(1)
  end function ran
 
  subroutine rangen(r)
    implicit none
    !---random number generator
    !   uses method of l'Ecuyer, (via F.James, Comp Phys Comm 60(1990)329)
    !   returns the vector r(:) of random values
    real(dp), intent(out) :: r(:)
    integer i,k,iz
    do i=1,size(r)
       k=iseed(1)/53668
       iseed(1)=40014*(iseed(1)-k*53668)-k*12211
       if (iseed(1).lt.0) iseed(1)=iseed(1)+2147483563
       k=iseed(2)/52774
       iseed(2)=40692*(iseed(2)-k*52774)-k*3791
       if (iseed(2).lt.0) iseed(2)=iseed(2)+2147483399
       iz=iseed(1)-iseed(2)
       if (iz.lt.1) iz=iz+2147483562
       r(i)=real(iz,kind=dp)*4.656613d-10
    enddo
  end subroutine rangen

  
  subroutine SetSeed(iseed_in)
    integer, intent(in) :: iseed_in(2)
    iseed = iseed_in
  end subroutine SetSeed
  

  subroutine GetSeed(iseed_out)
    integer, intent(out) :: iseed_out(2)
    iseed_out = iseed
  end subroutine GetSeed

end module random
