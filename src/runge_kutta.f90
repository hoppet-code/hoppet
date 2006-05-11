! $Id: runge_kutta.f90,v 1.1 2001/06/27 13:40:17 gsalam Exp $
module runge_kutta
  use types; use consts_dp
  implicit none
  private


  real(dp),  parameter :: third = one/three

  interface rkstp
     module procedure rkstp_0d, rkstp_1d, rkstp_2d
  end interface
  public :: rkstp


contains
  !----------------------------------------------------------------------
  ! write the Runge Kutta routine; try and be efficient with work-space
  ! This is a scalar version...
  subroutine rkstp_0d(h,x,y,conv)
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: x, y
    interface
       subroutine conv(x,y,dy)
         use types
         real(dp), intent(in)  :: x, y
         real(dp), intent(out) :: dy
       end subroutine conv
    end interface
    !------------------------------------------------------------
    real(dp) :: w1, w2, w3
    real(dp) :: hh
    hh = half * h
    call conv(x,y,w1); w1 = w1 * hh            ! w1 = k1/2
    call conv(x+hh, y + w1, w2); w2 = w2 * hh  ! w2 = k2/2
    call conv(x+hh, y + w2, w3); w3 = w3 * h   ! w3 = k3
    w2 = w1 + two*w2                           ! w2 = half*k1 + k2
    call conv(x+h , y + w3, w1); w1 = w1 * hh  ! w1 = k4/2
    
    x = x + h
    !                k1/2 + k2  + k3 + k4/2
    y = y + third * (w2         + w3 + w1)
  end subroutine rkstp_0d

  !----------------------------------------------------------------------
  subroutine rkstp_1d(h,x,y,conv)
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: x, y(:)
    interface
       subroutine conv(x,y,dy)
         use types
         real(dp), intent(in)  :: x, y(:)
         real(dp), intent(out) :: dy(:)
       end subroutine conv
    end interface
    !------------------------------------------------------------
    real(dp) :: w1(size(y)), w2(size(y)), w3(size(y))
    real(dp) :: hh
    hh = half * h
    call conv(x,y,w1); w1 = w1 * hh            ! w1 = k1/2
    call conv(x+hh, y + w1, w2); w2 = w2 * hh  ! w2 = k2/2
    call conv(x+hh, y + w2, w3); w3 = w3 * h   ! w3 = k3
    w2 = w1 + two*w2                           ! w2 = half*k1 + k2
    call conv(x+h , y + w3, w1); w1 = w1 * hh  ! w1 = k4/2
    
    x = x + h
    !                k1/2 + k2  + k3 + k4/2
    y = y + third * (w2         + w3 + w1)
  end subroutine rkstp_1d

  !----------------------------------------------------------------------
  subroutine rkstp_2d(h,x,y,conv)
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: x, y(:,:)
    interface
       subroutine conv(x,y,dy)
         use types
         real(dp), intent(in)  :: x, y(:,:)
         real(dp), intent(out) :: dy(:,:)
       end subroutine conv
    end interface
    !------------------------------------------------------------
    real(dp) :: w1(size(y,dim=1),size(y,dim=2)), &
         &      w2(size(y,dim=1),size(y,dim=2)), &
         &      w3(size(y,dim=1),size(y,dim=2))
    real(dp) :: hh
    hh = half * h
    call conv(x,y,w1); w1 = w1 * hh            ! w1 = k1/2
    call conv(x+hh, y + w1, w2); w2 = w2 * hh  ! w2 = k2/2
    call conv(x+hh, y + w2, w3); w3 = w3 * h   ! w3 = k3
    w2 = w1 + two*w2                           ! w2 = half*k1 + k2
    call conv(x+h , y + w3, w1); w1 = w1 * hh  ! w1 = k4/2
    
    x = x + h
    !                k1/2 + k2  + k3 + k4/2
    y = y + third * (w2         + w3 + w1)
  end subroutine rkstp_2d

end module runge_kutta
