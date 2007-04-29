!----------------------------------------------------------------------
! A variety of routines that need to be implementation specific
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Reading and writing very large 2-dim arrays (with unformatted i/o) causes
! problems on certain systems -- so provide routines which split up
! the reading and writing if that is necessary
! Splitting not required for dec.
!----------------------------------------------------------------------
subroutine lcl_write_2d(idev, array)
      implicit none
      integer,         intent(in) :: idev
      real(kind(1d0)), intent(in) :: array(1:, 1:)
      !------------------------------------------------------------
      integer  ulim(2), i

      ulim = ubound(array)
      do i = 1, ulim(2)
         write(idev) real(array(:,i))
      end do
      
end subroutine lcl_write_2d


subroutine lcl_read_2d(idev, array, ifail)
      implicit none
      integer,           intent(in)  :: idev
      real(kind(1d0)),   intent(out) :: array(1:, 1:)
      integer, optional, intent(out) :: ifail
      !------------------------------------------------------------
      real, allocatable :: sngl_prec_buffer(:)
      integer  ulim(2), i!, j

      ulim = ubound(array)
      allocate(sngl_prec_buffer(1:ulim(1)))

      do i = 1, ulim(2)
         !-- data stored in sngl precision
         read(idev, end=10) sngl_prec_buffer(:)
         !-- convert it back into double precision
         array(:,i) = dble(sngl_prec_buffer(:))
      end do

      !-- correct return code ---
      if (present(ifail)) then
         ifail = 0
         return
      end if

10    if (present(ifail)) then
         ifail = -1
         return
      end if
end subroutine lcl_read_2d

!!$subroutine lcl_write_2d(idev, array)
!!$      implicit none
!!$      integer, intent(in) :: idev
!!$      real,    intent(in) :: array(1:, 1:)
!!$      !------------------------------------------------------------
!!$
!!$      write(idev) array(:,:)
!!$end subroutine lcl_write_2d
!!$
!!$
!!$subroutine lcl_read_2d(idev, array)
!!$      implicit none
!!$      integer, intent(in)  :: idev
!!$      real,    intent(out) :: array(1:, 1:)
!!$      !------------------------------------------------------------
!!$
!!$      read(idev) array(:,:)
!!$end subroutine lcl_read_2d



!----------------------------------------------------------------------
! Interfaces to the dec f90 iargc and getarg routines
!
! GPS 4/11/95 (CCN8 9)
!----------------------------------------------------------------------
integer function lcl_iargc()
      implicit none
      integer iargc

      lcl_iargc = iargc()
end function lcl_iargc

subroutine lcl_getarg(k, argk)
      implicit none
      integer,      intent(in)  :: k
      character(*), intent(out) :: argk

      call getarg(k, argk)
end subroutine lcl_getarg


subroutine lcl_flush(idev)
  implicit none
  integer, intent(in) :: idev
  
  call flush(idev)
end subroutine lcl_flush

subroutine lcl_system(string)
      implicit none
      character(*), intent(in) :: string
      !------------------------------------------------------------
      !integer return_val, system

      call  system(string)
end subroutine lcl_system

      
